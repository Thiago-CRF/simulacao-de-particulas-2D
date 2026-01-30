#include "raylib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// gcc -Wall -Wextra -g3 -std=c2x -lm particle.c -o output/particle_sim -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 && ./output/particle_sim

#define WIDTH 900
#define HEIGHT 600
#define FPS 100

// definições da simulação
#define PARTICLE_NUM 15
#define INIT_RADIUS 20
#define INIT_VELOCITY 5

typedef struct
{
   float x, y, rad;
   float velocity_x, velocity_y;
} Particle;

void initialize_particles(Particle *particles)
{
    srand(time(NULL));
    /* range é pra calcular um numero aleatório entre um minimo e máximo 
    sendo a velocidade entre (-INIT_VELOCITY) e INIT_VELOCITY
    e o raio entre 5 e INIT_RADIUS*/
    int vel_range = (INIT_VELOCITY - (-INIT_VELOCITY) + 1);
    int rad_range = (INIT_RADIUS - (5) + 1);
    /* inicia a particula com posições aleatórias e com velocidade aleatória 
        entre (-velocidade inicial) e velocidade inicial */
    for(int i=0; i<PARTICLE_NUM; i++)
    {
        particles[i] = (Particle){(rand()%WIDTH), (rand()%HEIGHT), 
            (rand()%rad_range + (5)), 
            (rand()%vel_range + (-INIT_VELOCITY)), 
            (rand()%vel_range + (-INIT_VELOCITY))};
    }
}

void update_particle(Particle *particle)
{
    particle->x += particle->velocity_x;
    particle->y += particle->velocity_y;

    float x = particle->x;
    float y = particle->y;
    float rad = particle->rad;

    /* colisão com as paredes, corrige a posição 
    se for pra fora de parede e reflete a particula */
    if((x-rad) < 0)
    {
        particle->x = rad;
        particle->velocity_x = -particle->velocity_x;
    }
    if((x+rad) > WIDTH)
    {
        particle->x = WIDTH - rad;
        particle->velocity_x = -particle->velocity_x;
    }
    if((y-rad) < 0)
    {
        particle->y = rad;
        particle->velocity_y = -particle->velocity_y;
    }
    if((y+rad) > HEIGHT)
    {
        particle->y = HEIGHT - rad;
        particle->velocity_y = -particle->velocity_y;
    }
}

void draw_particles(Particle *particles)
{
    for(int i=0; i<PARTICLE_NUM; i++)
    {
        DrawCircle(particles[i].x, particles[i].y, particles[i].rad, WHITE);
    }
}

void update_environment(Particle *particles)
{
    for(int i=0; i<PARTICLE_NUM; i++)
    {
        update_particle(&particles[i]);
    }
}

int main()
{
    InitWindow(WIDTH, HEIGHT, "Simulação de particulas");
    SetTargetFPS(FPS);

    Particle *particles = (Particle*)malloc(PARTICLE_NUM * sizeof(Particle));
    initialize_particles(particles);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);
        DrawFPS(10,10);

        update_environment(particles);

        draw_particles(particles);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}