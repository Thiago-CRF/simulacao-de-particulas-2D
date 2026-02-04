#include "raylib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// gcc -Wall -Wextra -g3 -std=c2x -lm particle.c -o output/particle_sim -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 && ./output/particle_sim

#define WIDTH 900
#define HEIGHT 600
#define FPS 200

// definições da simulação
#define PARTICLE_NUM 8
#define INIT_RADIUS 20
#define INIT_VELOCITY 120
#define DENSITY 1

typedef struct
{
   float x, y, rad;
   double velocity_x, velocity_y;
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

void wall_particle_collison(Particle *particle)
{
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

void particle_particle_collision(Particle *current_particle, Particle *particle_array, int pos)
{
    Vector2 curVect = (Vector2){current_particle->x, current_particle->y};
    Vector2 otherVect;

    // verifica a colisão da particula atual com cada uma das outras particulas(particle_array[i])
    for(int i=pos+1; i<PARTICLE_NUM; i++)
    {
        // pula se estiver comparando a particula com ela mesma
        if(current_particle == &particle_array[i])
            continue;
        
        otherVect = (Vector2){particle_array[i].x, particle_array[i].y};
        
        bool check = CheckCollisionCircles(curVect, current_particle->rad, 
                        otherVect, particle_array[i].rad);
        if(check)
        {
            // calculo da direção das particulas rebatendo

            /* 
            Matemática da colisão:
            em 1D: velocidade_final(pós colisão)_A e velocidade_final_B é calculado assim:
            
            mA e mB = massa de A e de B. vamos usar o raio pra calcular a massa.
            massa = densidade*volume
            (fazer cálculo da massa do circulo com base no raio e densidade)
            vA2: velocidade final de A. vA1: velocidade inicial de A;
            vB2: velocidade final de B. vB1: velocidade inicial de B;

            vA2 = (((mA - mB) / (mA + mB))*vA1) + (((2 * mB) / (mA + mB))*vB1)
            vB2 = (((2 * mA) / (mA + mB))*vA1) + (((mB - mA) / (mA + mB))*vB1)

            Em 2D temos que decompor a velocidade em x e y pra fazer o calculo do vetor
            e calcular o momento angular da colisão

            a velocidade geral de cada particula tem que ser dividida em em duas velocidades
            perpendiculares, uma tangente ao ponto de colisão e outra perpendicular a linha 
            tangente ao ponto de colisão.

            preciso calcular o angulo de deslocamento inicial das particula em colisão
            usando o x e y da velocidade.
            usando a função atan2 eu tenho o angulo da direção de x e y, entre -PI e PI

            angulo_inicial = atan2(velocity_x, velocity_y) obs: atan2(vel_y, vel_x), poias atan2 recebe os valores assim
            
            e também preciso calcular o angulo de contato entre as particulas, pra tagente.
            pra isso

            encontrei um jeito de calcular o angulo de contato assim:
            ang_ctt = atan2((y2-y1), (x2-x1))
            sendo x1, x2 e y1, y2 as posições da particula 1 e 2

            no fim, pra calcular as velocidades x e y das particulas 1 e 2 vai ser:

            vFx_1' = ((vAx * cos(ang_ctt) + vAy * sen(ang_ctt)) * (mA-mB)) + (2*mB * (vBx * cos(ang_ctt) + vBy * sen(ang_ctt)))
            vFx_1 = ((vFx_1' / (mA+mB)) * cos(ang_ctt)) + ((vAy * cos(ang_ctt)) - (vAx * sen(ang_ctt)) * cos(ang_ctt + PI/2))

            vFy_1' = ((vAx * cos(ang_ctt) + vAy * sen(ang_ctt)) * (mA-mB)) + (2*mB * (vBx * cos(ang_ctt) + vBy * sen(ang_ctt)))
            vFy_1 = ((vFy_1' / (mA+mB)) * sen(ang_ctt)) + ((vAy * cos(ang_ctt)) - vAx * sen(ang_ctt)) * sen(ang_ctt + PI/2)
        
            v1 e v2 são as velocidades escalares das particulas, hypot(vel_x, vel_y)
            */

            // primeiro tem que fazer com que as particulas fiquem em uma posição pré colisão, pra que não tenha bug de calcular
            // mais de uma vez por conta de uma estar minimamente dentro da outra. Mas antes pego a posição dela pra calcular o angulo de colisão depois
            float x_p1 = current_particle->x;
            float y_p1 = current_particle->y;
            float x_p2 = particle_array[i].x;
            float y_p2 = particle_array[i].y;
            
            // calculo de sobreposição das particulas
            float dx = x_p2-x_p1;
            float dy = y_p2-y_p1;
            float dist = sqrt((dx*dx) + (dy*dy));

            float overlap = (current_particle->rad + particle_array[i].rad) - dist;
            if(overlap > 0)
            {
                float move_x = (dx/dist) * (overlap / 2.0f);
                float move_y = (dy/dist) * (overlap / 2.0f);

                current_particle->x -= move_x;
                current_particle->y -= move_y;
                particle_array[i].x += move_x;
                particle_array[i].y += move_y;

                x_p1 = current_particle->x;
                y_p1 = current_particle->y;
                x_p2 = particle_array[i].x;
                y_p2 = particle_array[i].y;
            }

           // primeiro, calculo a massa
            double m1 = DENSITY * PI * pow(current_particle->rad, 2)/1000;
            double m2 = DENSITY * PI * pow(particle_array[i].rad, 2)/1000;

            // depois, calculo de angulos
            // angulos iniciais de p1 e p2
            double ang_p1 = atan2(current_particle->velocity_y, current_particle->velocity_x);
            double ang_p2 = atan2(particle_array[i].velocity_y, particle_array[i].velocity_x);

            // angulo da linha que liga os centros das particulas, em relação ao eixo x
            double ang_ctt = atan2((y_p2-y_p1),(x_p2-x_p1));

            // calculo da velocidade escalar
            double v1 = hypot(current_particle->velocity_x, current_particle->velocity_y);
            double v2 = hypot(particle_array[i].velocity_x, particle_array[i].velocity_y);

            // por fim, calculo das velocidades após colisão
            double velF_x1 = (((v1 * cos(ang_p1-ang_ctt) * (m1-m2)) + (2*m2*v2 * cos(ang_p2-ang_ctt))) / (m1+m2)) * (cos(ang_ctt) + v1*sin(ang_p1-ang_ctt) * cos(ang_ctt + PI/2));
            double velF_y1 = (((v1 * cos(ang_p1-ang_ctt) * (m1-m2)) + (2*m2*v2 * cos(ang_p2-ang_ctt))) / (m1+m2)) * (sin(ang_ctt) + v1*sin(ang_p1-ang_ctt) * sin(ang_ctt + PI/2));
            double velF_x2 = (((v2 * cos(ang_p2-ang_ctt) * (m2-m1)) + (2*m1*v1 * cos(ang_p1-ang_ctt))) / (m2+m1)) * (cos(ang_ctt) + v2*sin(ang_p2-ang_ctt) * cos(ang_ctt + PI/2));
            double velF_y2 = (((v2 * cos(ang_p2-ang_ctt) * (m2-m1)) + (2*m1*v1 * cos(ang_p1-ang_ctt))) / (m2+m1)) * (sin(ang_ctt) + v2*sin(ang_p2-ang_ctt) * sin(ang_ctt + PI/2));

            current_particle->velocity_x = velF_x1;
            current_particle->velocity_y = velF_y1;
            particle_array[i].velocity_x = velF_x2;
            particle_array[i].velocity_y = velF_y2; 
            printf("colisao");
        }
    }
}

void update_movement(Particle *particle)
{
    double delta_t = GetFrameTime();

    particle->x += particle->velocity_x * delta_t;
    particle->y += particle->velocity_y * delta_t;
}

void draw_particles(Particle *particles)
{
    for(int i=0; i<PARTICLE_NUM; i++)
    {
        DrawCircle(particles[i].x, particles[i].y, particles[i].rad, WHITE);
    }
}

void update_particles(Particle *particles)
{
    for(int i=0; i<PARTICLE_NUM; i++)
    {
        wall_particle_collison(&particles[i]);
        particle_particle_collision(&particles[i], particles, i);
        update_movement(&particles[i]);
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

        update_particles(particles);
        draw_particles(particles);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}