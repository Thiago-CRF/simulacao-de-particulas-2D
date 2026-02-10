Simulação simples de um "particulas" circulares em 2D feita em C usando a biblioteca externa 'Raylib'

Variáveis de definição da simulação estão no começo do código, nos #define.

Tamanho, posição e velocidades iniciais das particulas são aleatórios com base nos valores iniciais no #define
Velocidade é entre: -INIT_VELOCITY e +INIT_VELOCITY
Raio da particula é entre: 5 e INIT_RADIUS
Posição é qualquer posição da tela

Precisa-se instalar a biblioteca no computador 'https://github.com/raysan5/raylib/wiki'  
Para compilar, precisa usar as flags da biblioteca raylib. 

Tem um makefile com as flags da biblioteca no repositório (-lraylib -lGL -lm -lpthread -ldl -lrt -lX11), é mais prático para compilar o programa (Só usar "make run" no diretório da pasta).

Na pasta /output tem o programa executavel para linux.