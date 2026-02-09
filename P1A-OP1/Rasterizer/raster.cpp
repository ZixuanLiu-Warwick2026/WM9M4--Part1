#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"

// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;

    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal; 
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
            t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        // Create a triangle object and render it
        triangle tri(t[0], t[1], t[2]);
        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
    Renderer renderer;
    // create light source {direction, diffuse intensity, ambient intensity}
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };
    // camera is just a matrix
    matrix camera = matrix::makeIdentity(); // Initialize the camera with identity matrix

    bool running = true; // Main loop control variable

    std::vector<Mesh*> scene; // Vector to store scene objects

    // Create a sphere and a rectangle mesh
    Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
    //Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

    // add meshes to scene
    scene.push_back(&mesh);
   // scene.push_back(&mesh2); 

    float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
    mesh.world = matrix::makeTranslation(x, y, z);
    //mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput(); // Handle user input
        renderer.clear(); // Clear the canvas for the next frame

        // Apply transformations to the meshes
     //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
        mesh.world = matrix::makeTranslation(x, y, z);

        // Handle user inputs for transformations
        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
        if (renderer.canvas.keyPressed('A')) x += -0.1f;
        if (renderer.canvas.keyPressed('D')) x += 0.1f;
        if (renderer.canvas.keyPressed('W')) y += 0.1f;
        if (renderer.canvas.keyPressed('S')) y += -0.1f;
        if (renderer.canvas.keyPressed('Q')) z += 0.1f;
        if (renderer.canvas.keyPressed('E')) z += -0.1f;

        // Render each object in the scene
        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present(); // Display the rendered frame
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix::makeIdentity();
    }
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    std::cout << "######### scene 1 / Mode: Release #########" << std::endl;
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    bool running = true;

    std::vector<Mesh*> scene;

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        for (auto& m : scene)
            render(renderer, m, camera, L);
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    std::cout << "######### scene 2 / Mode: Release #########" << std::endl;
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    std::vector<Mesh*> scene;

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.push_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        for (auto& m : scene)
            render(renderer, m, camera, L);
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene3: nebula-like particle cloud (only cubes + spheres), dense near the center,
// orbiting the center with uniform angular speed.
void scene3() {
    std::cout << "######### scene 3 / Mode: Release #########" << std::endl;
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    std::vector<Mesh*> scene;

    // --- Cloud settings ---
    const unsigned int N = 1800;            // thousands of small objects
    const float cx = 0.f, cy = 0.f, cz = -12.f;   // cloud center (in front of the camera)
    const float sigma = 3.0f;               // controls the size of the cloud (smaller = more concentrated)
    float angleY = 0.f;
    float angleX = 0.f;
    const float wY = 0.010f;                // constant angular velocity of rotation
    const float wX = 0.006f;

    // Use "sum of 12 uniform random numbers minus 6" to approximate a Gaussian distribution
    // (denser at the center, like a nebula)
    auto gaussian01 = [&]() -> float {
        float s = 0.f;
        for (int i = 0; i < 12; i++) s += rng.getRandomFloat(0.f, 1.f);
        return s - 6.f; // ~ N(0,1)
        };

    struct Obj {
        Mesh* m;
        float ox, oy, oz;     // offset relative to the center (in the cloud's coordinate system)
        float selfSpin;       // slight self-rotation (optional)
    };
    std::vector<Obj> objs;
    objs.reserve(N);

    // --- Create particle cloud ---
    for (unsigned int i = 0; i < N; i++) {
        Mesh* m = new Mesh();

        // Most use small spheres to look more like "particles";
        // a few small cubes add some visual texture
        if (rng.getRandomFloat(0.f, 1.f) < 0.80f) {
            *m = Mesh::makeSphere(0.15f, 5, 9);     // low tessellation: faster
        }
        else {
            *m = Mesh::makeCube(0.25f);
        }

        // 3D Gaussian distribution: dense at the center, sparse towards the edges
        float ox = gaussian01() * sigma;
        float oy = gaussian01() * sigma;
        float oz = gaussian01() * (0.8f * sigma);   // slightly flattened in depth to look more cloud-like

        // Prevent a few extreme outliers (clip points that are too far away)
        const float maxR = 7.0f;
        if (ox * ox + oy * oy + oz * oz > maxR * maxR) {
            // simple resampling once (good enough)
            ox = gaussian01() * sigma;
            oy = gaussian01() * sigma;
            oz = gaussian01() * (0.8f * sigma);
        }

        Obj o;
        o.m = m;
        o.ox = ox; o.oy = oy; o.oz = oz;
        o.selfSpin = rng.getRandomFloat(-0.03f, 0.03f);

        objs.push_back(o);
        scene.push_back(m);

        // initial position
        m->world = matrix::makeTranslation(cx, cy, cz) * matrix::makeTranslation(ox, oy, oz);
    }

    bool running = true;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0; // used as a counter (here represents print count / frame segment)


    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        // constant-speed rotation (the whole cloud orbits around the center)
        angleY += wY;
        angleX += wX;

        matrix cloudOrbit =
            matrix::makeTranslation(cx, cy, cz) *
            matrix::makeRotateY(angleY) *
            matrix::makeRotateX(angleX);

        // Update each particle: global rotation around the center
        // + a bit of self-rotation to enhance the "liveliness"
        for (unsigned int i = 0; i < objs.size(); i++) {
            Obj& o = objs[i];
            matrix local = matrix::makeTranslation(o.ox, o.oy, o.oz);
            matrix spin = matrix::makeRotateY(o.selfSpin);
            o.m->world = cloudOrbit * local * spin;
        }

        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present();

        // Print the total time cost for every 120 frames
        // (larger value means slower performance)
        if (++cycle % 120 == 0) {
            end = std::chrono::high_resolution_clock::now();
            std::cout << (cycle / 120) << " :"
                << std::chrono::duration<double, std::milli>(end - start).count()
                << "ms\n";
            start = std::chrono::high_resolution_clock::now();
        }

    }

    for (auto& m : scene)
        delete m;
}




// Entry point of the application
// No input variables
int main() {
    // Uncomment the desired scene function to run
    //scene1();
    //scene2();
    scene3();
    //sceneTest(); 
    

    return 0;
}