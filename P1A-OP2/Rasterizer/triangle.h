#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// Simple support class for a 2D vector
class vec2D {
public:
    float x, y;

    // Default constructor initializes both components to 0
    vec2D() { x = y = 0.f; };

    // Constructor initializes components with given values
    vec2D(float _x, float _y) : x(_x), y(_y) {}

    // Constructor initializes components from a vec4
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }

    // Display the vector components
    void display() { std::cout << x << '\t' << y << std::endl; }

    // Overloaded subtraction operator for vector subtraction
    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// Class representing a triangle for rendering purposes
class triangle {
    Vertex v[3];       // Vertices of the triangle
    float area;        // Area of the triangle
    colour col[3];     // Colors for each vertex of the triangle

public:
    // Constructor initializes the triangle with three vertices
    // Input Variables:
    // - v1, v2, v3: Vertices defining the triangle
    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        // Calculate the 2D area of the triangle
        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = std::fabs(e1.x * e2.y - e1.y * e2.x);
    }

    // Helper function to compute the cross product for barycentric coordinates
    // Input Variables:
    // - v1, v2: Edges defining the vector
    // - p: Point for which coordinates are being calculated
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    // Compute barycentric coordinates for a given point
    // Input Variables:
    // - p: Point to check within the triangle
    // Output Variables:
    // - alpha, beta, gamma: Barycentric coordinates of the point
    // Returns true if the point is inside the triangle, false otherwise
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        return true;
    }

    // Template function to interpolate values using barycentric coordinates
    // Input Variables:
    // - alpha, beta, gamma: Barycentric coordinates
    // - a1, a2, a3: Values to interpolate
    // Returns the interpolated value
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    // Draw the triangle on the canvas
    // Input Variables:
    // - renderer: Renderer object for drawing
    // - L: Light object for shading calculations
    // - ka, kd: Ambient and diffuse lighting coefficients
    void draw(Renderer& renderer, Light& L, float ka, float kd) {
        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);

        if (area < 1.f) return;

        // Optimisation 1: precompute loop bounds
        const int yStart = (int)(minV.y);
        const int yEnd = (int)ceil(maxV.y);
        const int xStart = (int)(minV.x);
        const int xEnd = (int)ceil(maxV.x);

        // Optimisation 1: normalise light direction once per draw()
        L.omega_i.normalise();

        // Optimisation 2: edge functions + incremental update
        const float invArea = 1.0f / area;

        const vec2D p0(v[0].p), p1(v[1].p), p2(v[2].p);

        const float e01x = p1.x - p0.x, e01y = p1.y - p0.y;
        const float e12x = p2.x - p1.x, e12y = p2.y - p1.y;
        const float e20x = p0.x - p2.x, e20y = p0.y - p2.y;

        const float e01_dx = -e01y, e01_dy = e01x;
        const float e12_dx = -e12y, e12_dy = e12x;
        const float e20_dx = -e20y, e20_dy = e20x;

        auto edgeAt = [](float v1x, float v1y, float ex, float ey, float px, float py) {
            return ex * (py - v1y) - ey * (px - v1x);
            };

        const float fx0 = (float)xStart;
        const float fy0 = (float)yStart;

        float rowE01 = edgeAt(p0.x, p0.y, e01x, e01y, fx0, fy0);
        float rowE12 = edgeAt(p1.x, p1.y, e12x, e12y, fx0, fy0);
        float rowE20 = edgeAt(p2.x, p2.y, e20x, e20y, fx0, fy0);

        for (int y = yStart; y < yEnd; y++) {
            float E01 = rowE01;
            float E12 = rowE12;
            float E20 = rowE20;

            for (int x = xStart; x < xEnd; x++) {

                if (E01 >= 0.f && E12 >= 0.f && E20 >= 0.f) {

                    const float alpha = E01 * invArea;
                    const float beta = E12 * invArea;
                    const float gamma = E20 * invArea;

                    // Optimisation 1: depth first
                    const float depth =
                        v[0].p[2] * beta +
                        v[1].p[2] * gamma +
                        v[2].p[2] * alpha;

                    // Optimisation 1: single zbuffer access
                    float& z = renderer.zbuffer(x, y);
                    if (z > depth && depth > 0.001f) {

                        vec4 normal =
                            v[0].normal * beta +
                            v[1].normal * gamma +
                            v[2].normal * alpha;
                        normal.normalise();

                        colour c =
                            v[0].rgb * beta +
                            v[1].rgb * gamma +
                            v[2].rgb * alpha;
                        c.clampColour();

                        float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                        colour a = (c * kd) * (L.L * dot) + (L.ambient * ka);

                        unsigned char r, g, b;
                        a.toRGB(r, g, b);
                        renderer.canvas.draw(x, y, r, g, b);

                        z = depth; // write back once
                    }
                }

                E01 += e01_dx;
                E12 += e12_dx;
                E20 += e20_dx;
            }

            rowE01 += e01_dy;
            rowE12 += e12_dy;
            rowE20 += e20_dy;
        }
    }





    // Compute the 2D bounds of the triangle
    // Output Variables:
    // - minV, maxV: Minimum and maximum bounds in 2D space
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = std::min(minV.x, v[i].p[0]);
            minV.y = std::min(minV.y, v[i].p[1]);
            maxV.x = std::max(maxV.x, v[i].p[0]);
            maxV.y = std::max(maxV.y, v[i].p[1]);
        }
    }

    // Compute the 2D bounds of the triangle, clipped to the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    // Output Variables:
    // - minV, maxV: Clipped minimum and maximum bounds
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = std::max(minV.x, static_cast<float>(0));
        minV.y = std::max(minV.y, static_cast<float>(0));
        maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
        maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
    }

    // Debugging utility to display the triangle bounds on the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};
