#include <SFML/Graphics.h>
#include <SFML/Window.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IX(x, y, z, N) ((x) + (y) * (N) + (z) * (N) * (N))
#define SWAP(x0, x) {float *tmp = x0; x0 = x; x = tmp;}
#define MAX_SOLVER_ITERATIONS 100

typedef struct {
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;
    float *Vz;

    float *Vx0;
    float *Vy0;
    float *Vz0;
} FluidCube;

// Function prototypes
FluidCube* FluidCubeCreate(int size, float diffusion, float viscosity, float dt);
void FluidCubeFree(FluidCube *cube);
void FluidCubeStep(FluidCube *cube);
void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, float amount);
void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ);
static void set_bnd(int b, float *x, int N);
static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N);
static void diffuse(int b, float *x, float *x0, float diff, float dt, int iter, int N);
static void advect(int b, float *d, float *d0, float *velocX, float *velocY, float *velocZ, float dt, int N);
static void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter, int N);

// FluidCube functions
FluidCube* FluidCubeCreate(int size, float diffusion, float viscosity, float dt)
{
    FluidCube *cube = malloc(sizeof(*cube));
    
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = calloc(size * size * size, sizeof(float));
    cube->density = calloc(size * size * size, sizeof(float));
    
    cube->Vx = calloc(size * size * size, sizeof(float));
    cube->Vy = calloc(size * size * size, sizeof(float));
    cube->Vz = calloc(size * size * size, sizeof(float));
    
    cube->Vx0 = calloc(size * size * size, sizeof(float));
    cube->Vy0 = calloc(size * size * size, sizeof(float));
    cube->Vz0 = calloc(size * size * size, sizeof(float));
    
    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    free(cube->s);
    free(cube->density);
    free(cube->Vx);
    free(cube->Vy);
    free(cube->Vz);
    free(cube->Vx0);
    free(cube->Vy0);
    free(cube->Vz0);
    free(cube);
}

static void set_bnd(int b, float *x, int N) {
    // ... (same as before, no changes)
     for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0, N)] = b == 3 ? -x[IX(i, j, 1, N)] : x[IX(i, j, 1, N)];
            x[IX(i, j, N-1, N)] = b == 3 ? -x[IX(i, j, N-2, N)] : x[IX(i, j, N-2, N)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0, k, N)] = b == 2 ? -x[IX(i, 1, k, N)] : x[IX(i, 1, k, N)];
            x[IX(i, N-1, k, N)] = b == 2 ? -x[IX(i, N-2, k, N)] : x[IX(i, N-2, k, N)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0, j, k, N)] = b == 1 ? -x[IX(1, j, k, N)] : x[IX(1, j, k, N)];
            x[IX(N-1, j, k, N)] = b == 1 ? -x[IX(N-2, j, k, N)] : x[IX(N-2, j, k, N)];
        }
    }

    // Handle corners
    x[IX(0, 0, 0, N)] = 0.33f * (x[IX(1, 0, 0, N)] + x[IX(0, 1, 0, N)] + x[IX(0, 0, 1, N)]);
    x[IX(0, N-1, 0, N)] = 0.33f * (x[IX(1, N-1, 0, N)] + x[IX(0, N-2, 0, N)] + x[IX(0, N-1, 1, N)]);
    x[IX(0, 0, N-1, N)] = 0.33f * (x[IX(1, 0, N-1, N)] + x[IX(0, 1, N-1, N)] + x[IX(0, 0, N-2, N)]);
    x[IX(0, N-1, N-1, N)] = 0.33f * (x[IX(1, N-1, N-1, N)] + x[IX(0, N-2, N-1, N)] + x[IX(0, N-1, N-2, N)]);
    x[IX(N-1, 0, 0, N)] = 0.33f * (x[IX(N-2, 0, 0, N)] + x[IX(N-1, 1, 0, N)] + x[IX(N-1, 0, 1, N)]);
    x[IX(N-1, N-1, 0, N)] = 0.33f * (x[IX(N-2, N-1, 0, N)] + x[IX(N-1, N-2, 0, N)] + x[IX(N-1, N-1, 1, N)]);
    x[IX(N-1, 0, N-1, N)] = 0.33f * (x[IX(N-2, 0, N-1, N)] + x[IX(N-1, 1, N-1, N)] + x[IX(N-1, 0, N-2, N)]);
    x[IX(N-1, N-1, N-1, N)] = 0.33f * (x[IX(N-2, N-1, N-1, N)] + x[IX(N-1, N-2, N-1, N)] + x[IX(N-1, N-1, N-2, N)]);

}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N) {
    // ... (same as before, no changes)
    float cRecip = 1.0f / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m, N)] =
                        (x0[IX(i, j, m, N)]
                            + a*(    x[IX(i+1, j, m, N)] + x[IX(i-1, j, m, N)]
                                    + x[IX(i, j+1, m, N)] + x[IX(i, j-1, m, N)]
                                    + x[IX(i, j, m+1, N)] + x[IX(i, j, m-1, N)])) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse(int b, float *x, float *x0, float diff, float dt, int iter, int N) {
    // ... (same as before, no changes)
     float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, float *d, float *d0, float *velocX, float *velocY, float *velocZ, float dt, int N) {
    // ... (same as before, no changes)
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);
    
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                float i0, i1, j0, j1, k0, k1;
                float tmp;

                i0 = (i - dtx * velocX[IX(i, j, k, N)]) - 0.5f;
                i1 = i0 + 1.0f;
                j0 = (j - dty * velocY[IX(i, j, k, N)]) - 0.5f;
                j1 = j0 + 1.0f;
                k0 = (k - dtz * velocZ[IX(i, j, k, N)]) - 0.5f;
                k1 = k0 + 1.0f;

                tmp = d0[IX((int)i0, (int)j0, (int)k0, N)];
                tmp += (i1 - i0) * (d0[IX((int)i1, (int)j0, (int)k0, N)] - tmp);
                tmp += (j1 - j0) * (d0[IX((int)i0, (int)j1, (int)k0, N)] - tmp);
                tmp += (k1 - k0) * (d0[IX((int)i0, (int)j0, (int)k1, N)] - tmp);

                d[IX(i, j, k, N)] = tmp;
            }
        }
    }
    set_bnd(b, d, N);
}

static void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter, int N) {
    // ... (same as before, no changes)
 for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k, N)] = -0.5f * (
                    (velocX[IX(i+1, j, k, N)] - velocX[IX(i-1, j, k, N)]) +
                    (velocY[IX(i, j+1, k, N)] - velocY[IX(i, j-1, k, N)]) +
                    (velocZ[IX(i, j, k+1, N)] - velocZ[IX(i, j, k-1, N)])
                ) / N;
            }
        }
    }
    set_bnd(0, div, N);
    lin_solve(0, p, div, 1, 6, iter, N);

    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k, N)] -= 0.5f * (p[IX(i+1, j, k, N)] - p[IX(i-1, j, k, N)]) * N;
                velocY[IX(i, j, k, N)] -= 0.5f * (p[IX(i, j+1, k, N)] - p[IX(i, j-1, k, N)]) * N;
                velocZ[IX(i, j, k, N)] -= 0.5f * (p[IX(i, j, k+1, N)] - p[IX(i, j, k-1, N)]) * N;
            }
        }
    }
}

void FluidCubeStep(FluidCube *cube)
{
    int N = cube->size;
    float visc = cube->visc;
    float diff = cube->diff;
    float dt = cube->dt;
    
    diffuse(1, cube->Vx0, cube->Vx, visc, dt, MAX_SOLVER_ITERATIONS, N);
    diffuse(2, cube->Vy0, cube->Vy, visc, dt, MAX_SOLVER_ITERATIONS, N);
    diffuse(3, cube->Vz0, cube->Vz, visc, dt, MAX_SOLVER_ITERATIONS, N);
    
    project(cube->Vx0, cube->Vy0, cube->Vz0, cube->Vx, cube->Vy, MAX_SOLVER_ITERATIONS, N);
    
    advect(1, cube->Vx, cube->Vx0, cube->Vx0, cube->Vy0, cube->Vz0, dt, N);
    advect(2, cube->Vy, cube->Vy0, cube->Vx0, cube->Vy0, cube->Vz0, dt, N);
    advect(3, cube->Vz, cube->Vz0, cube->Vx0, cube->Vy0, cube->Vz0, dt, N);
    
    SWAP(cube->Vx0, cube->Vx);
    SWAP(cube->Vy0, cube->Vy);
    SWAP(cube->Vz0, cube->Vz);
    
    diffuse(0, cube->s, cube->density, diff, dt, MAX_SOLVER_ITERATIONS, N);
    advect(0, cube->density, cube->s, cube->Vx, cube->Vy, cube->Vz, dt, N);
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, float amount)
{
    cube->density[IX(x, y, z, cube->size)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ)
{
     cube->Vx[IX(x, y, z, cube->size)] += amountX;
    cube->Vy[IX(x, y, z, cube->size)] += amountY;
    cube->Vz[IX(x, y, z, cube->size)] += amountZ;
}

int main() {
    // Fluid parameters
    int size = 50; // Size of the cube
    float diffusion = 0.0f;
    float viscosity = 0.0f;
    float dt = 0.01f;

    // Create the fluid cube
    FluidCube *cube = FluidCubeCreate(size, diffusion, viscosity, dt);

    // Create the SFML window
    sfRenderWindow *window = sfRenderWindow_create((sfVideoMode){size * 10, size * 10, 32}, "Fluid Simulation", sfResize | sfClose, NULL);
    sfRenderWindow_setFramerateLimit(window, 60); // Limit to 60 FPS

    while (sfRenderWindow_isOpen(window)) {
        // Handle events
        sfEvent event;
        while (sfRenderWindow_pollEvent(window, &event)) {
            if (event.type == sfEvtClosed)
                sfRenderWindow_close(window);
        }

        // Update the fluid simulation
        FluidCubeStep(cube);

        // Add density and velocity in the middle
        FluidCubeAddDensity(cube, size / 2, size / 2, size / 2, 10.0f);
        FluidCubeAddVelocity(cube, size / 2, size / 2, size / 2, 1.0f, 0.0f, 0.0f); // Example velocity

        // Clear window
        sfRenderWindow_clear(window, sfBlack);

        // Draw the fluid
        for (int x = 1; x < size - 1; x++) {
            for (int y = 1; y < size - 1; y++) {
                float density = cube->density[IX(x, y, size / 2, cube->size)]; // Use a slice of the cube for 2D visualization
                if (density > 0) {
                    sfColor color = sfColor_fromRGB(density * 5, density * 5, density * 5); // Scale for visualization
sfRectangleShape *rect = sfRectangleShape_create();
                    sfRectangleShape_setFillColor(rect, color);
                    sfRectangleShape_setSize(rect, (sfVector2f){10, 10}); // Adjust size
                    sfRectangleShape_setPosition(rect, (sfVector2f){x * 10, y * 10}); // Position in the window
                    sfRenderWindow_drawRectangleShape(window, rect, NULL);
                    sfRectangleShape_destroy(rect);
                }
            }
        }

        // Display the contents of the window
        sfRenderWindow_display(window);
    }

    // Clean up
    FluidCubeFree(cube);
    sfRenderWindow_destroy(window);
    return 0;
}

