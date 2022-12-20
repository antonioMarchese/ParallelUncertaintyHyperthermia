#include <random>
#include <chrono>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <thread>

const auto num_of_threads = std::thread::hardware_concurrency();

#define real        double
#define TAM         100
#define N           100
#define M           100
#define K           2
#define STEPS       6000
#define INTERVAL    600
#define EXPERIMENTS 100

using namespace std;

class velocidade {
  public:
    real angulo;
    real magnitude;
};

class initial_velocity {
  public:
    real vx;
    real vy;
};

// I'll have to create a MPI_Type, so I took off the methods and implemented functions in order to let the 'velocidade' class only with real params
real velocidade_x(real magnitude, real angulo, real incerteza_angulo, real incerteza_magnitude){
      return magnitude*incerteza_magnitude*cos(angulo*incerteza_angulo);
}

real velocidade_y(real magnitude, real angulo, real incerteza_angulo, real incerteza_magnitude){
  return magnitude*incerteza_magnitude*sin(angulo*incerteza_angulo);
}

void saveVTK(real *save, int savedStep, int comm_size)
{
    FILE *t;
    char filename[30];
    sprintf(filename, "./resultados/results%d.vtk", savedStep);
    t = fopen(filename, "w");
    fprintf(t, "# vtk DataFile Version 3.0\n");
    fprintf(t, "results.vtk\n");
    fprintf(t, "ASCII\n");
    fprintf(t, "DATASET RECTILINEAR_GRID\n");
    fprintf(t, "DIMENSIONS %d %d 1\n", TAM, TAM);

    fprintf(t, "X_COORDINATES %d double\n", TAM);
    for(int i=0; i<TAM; i++)
    {
        fprintf(t, "%f ", i*(0.1/TAM));
    }
    fprintf(t, "\n");

    fprintf(t, "Y_COORDINATES %d double\n", TAM);
    for(int j=0; j<TAM; j++)
    {
        fprintf(t, "%f ", j*(0.1/TAM));
    }
    fprintf(t, "\n");

    fprintf(t, "Z_COORDINATES 1 double\n");
    fprintf(t, "0");
    fprintf(t, "\n");

    fprintf(t, "POINT_DATA %d \n", TAM*TAM*1);
    fprintf(t, "FIELD FieldData 1 \n");
    fprintf(t, "Temperatura 1 %d double \n", TAM*TAM*1);
      for(int j=0; j<TAM; j++)
      {
        for(int i=0; i<TAM; i++)
        {
            fprintf(t, "%f \n", save[TAM*TAM*savedStep + TAM*j + i]/comm_size);
        }
      }
    fprintf(t, "\n");
    fclose(t);
}

void save_velocity_VTK(initial_velocity **save)
{
    FILE *t;
    char filename[30];
    sprintf(filename, "./resultados/velocidade_incial.vtk");
    t = fopen(filename, "w");
    fprintf(t, "# vtk DataFile Version 3.0\n");
    fprintf(t, "results.vtk\n");
    fprintf(t, "ASCII\n");
    fprintf(t, "DATASET RECTILINEAR_GRID\n");
    fprintf(t, "DIMENSIONS %d %d 1\n", TAM, TAM);

    fprintf(t, "X_COORDINATES %d double\n", TAM);
    for(int i=0; i<TAM; i++)
    {
        fprintf(t, "%f ", i*(0.1/TAM));
    }
    fprintf(t, "\n");

    fprintf(t, "Y_COORDINATES %d double\n", TAM);
    for(int j=0; j<TAM; j++)
    {
        fprintf(t, "%f ", j*(0.1/TAM));
    }
    fprintf(t, "\n");

    fprintf(t, "Z_COORDINATES 1 double\n");
    fprintf(t, "0");
    fprintf(t, "\n");

    fprintf(t, "POINT_DATA %d \n", TAM*TAM*1);
    fprintf(t, "VECTORS velocity double \n");
      for(int j=0; j<TAM; j++)
      {
        for(int i=0; i<TAM; i++)
        {
            fprintf(t, "%f %f 0\n", save[j][i].vx, save[j][i].vy);
        }
      }
    fprintf(t, "\n");
    fclose(t);
}

int main()
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    float pi = 3.14159265;
    std::default_random_engine generator;
    // obtain a seed from the timer
    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 = d.count();
    generator.seed (seed2);
    std::uniform_real_distribution<real> theta_distribution(0.7, 1.3);
    std::uniform_real_distribution<real> magnitude_distribution(0.5, 1.5);
    std::uniform_real_distribution<real> velocity_distribution(1.5*pow(10, -3), 4.5*pow(10, -3)); // v -> tumor
    std::uniform_real_distribution<real> theta(0, pi*2); // 0 -> vel < 0; 1 -> vel > 0
    std::uniform_real_distribution<real> is_capilar(0, 1);

    char str[40];
    velocidade *v = new velocidade[N*M];
    initial_velocity **initial_vel = new initial_velocity*[TAM];
    real p, p_b, c, c_b, *e, *qm, Ta = 37, capilar, incerteza_magnitude, incerteza_angulo, vx, vy;
    real hx = 0.001, ht = 0.5, *kappa, kip, kim, kjp, kjm, upwind_i, upwind_j, **Calor_nano, A = 0.5*pow(10, 6), r0 = 3.1*pow(10, -3);
    int i, j, x, y, n, lessX, lessY, plusX, plusY, pele, k, _swap = 0, **tecido, savedSteps = STEPS/INTERVAL, thread_count = num_of_threads;
    // MPI params
    int my_rank, comm_size, local_n;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //  Create the datatype
    MPI_Datatype velocity_type;

    // MPI Datatype params
    int block_lenght[2] = {1, 1};
    MPI_Aint displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(v, &base_address);
    MPI_Get_address(&v->angulo, &displacements[0]);
    MPI_Get_address(&v->magnitude, &displacements[1]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);

    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(2, block_lenght, displacements, types, &velocity_type);
    MPI_Type_commit(&velocity_type);

    //  Each process must be responsible for calculating (EXPERIMENTS/comm_size) experiments
    local_n = EXPERIMENTS/comm_size;

    real *u = new real[K*M*N], *u_new = new real[local_n*savedSteps*M*N], *mean = new real[savedSteps*M*N], aux, *sdv = new real[savedSteps*M*N];
    // 'u' -> um vetor de tamanho K*M*N, representando a matriz U, com
    //  N -> numero de linhas
    //  M -> numero de colunas
    //  K -> profundidade
    // 'u_new' = um vetor de tamanho local_n*savedSteps*M*N para salvar 'savedSteps' passos de tempo para cada experimento (local_n)
    //  'mean' = um vetor de tamanho savedSteps*M*N para salvar a media de 'savedSteps' passos de tempo
    //  'sdv' = um vetor de tamanho savedSteps*M*N para salvar o desvio padrão de 'savedSteps' passos de tempo

    real *final_mean = new real[savedSteps*M*N], *final_sdv = new real[savedSteps*M*N];
    // Assuming that these values dont change for healthy and diseased tissue
    p = 1000;
    p_b = 1000;
    c = 4000;
    c_b = 4000;
    e = new real[2];
    e[0] = 0.02;
    e[1] = 0.01;
    qm = new real[2];
    qm[0] = 420;
    qm[1] = 4200;
    kappa = new real[2];
    kappa[0] = 0.50;
    kappa[1] = 0.55;

    // Allocating memory

    tecido = new int*[TAM]; // 0 = healthy; 1 = diseased
    Calor_nano = new real*[TAM];

    for (i = 0; i < TAM; i++)
    {
        initial_vel[i] = new initial_velocity[TAM];
        Calor_nano[i] = new real[TAM];
        tecido[i] = new int[TAM];
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            Calor_nano[i][j] = A*exp(-(pow((i*hx - 0.05), 2) + pow((j*hx - 0.05), 2))/pow(r0, 2));
            tecido[i][j] = ((i*hx >= 0.04 + hx/2 && i*hx <= 0.06 + hx/2) && (j*hx >= 0.04 + hx/2 && j*hx <= 0.06 + hx/2)) ? 1 : 0;
        }
    }

    //  Only the master process will determine the magnitude and angle values of the velocity vector
    if (my_rank == 0)
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < M; j++)
            {
                is_capilar(generator) >= 0.5 ? capilar = 0 : capilar = 1;
                v[N*i + j].angulo = theta(generator);

                tecido[i][j] == 0 ? v[N*i + j].magnitude = (capilar*2*pow(10, -3)) : v[N*i + j].magnitude = (capilar*velocity_distribution(generator));

                initial_vel[i][j].vx = velocidade_x(v[N*i + j].magnitude, v[N*i + j].angulo, 1, 1);
                initial_vel[i][j].vy = velocidade_y(v[N*i + j].magnitude, v[N*i + j].angulo, 1, 1);
            }
        }
    }

    // MPI_Bcast(void* send_buf, int send_count, MPI_Datatype send_type, int source_proc, MPI_Comm);
    MPI_Bcast(v, N*M, velocity_type, 0, MPI_COMM_WORLD);

    // Saves the initial velocity
    if (my_rank == 0)
        save_velocity_VTK(initial_vel);
#   pragma omp parallel num_threads(thread_count) \
        default(none) private(k,_swap, j, i, n, x, y, vx, vy, upwind_i, upwind_j, lessX, plusX, kim, kip, lessY, plusY, kjm, kjp, pele) \
            shared(p, p_b, c, c_b, e, qm, u, local_n, u_new, v, tecido, hx, ht, Ta, Calor_nano, kappa, generator, incerteza_magnitude, incerteza_angulo, magnitude_distribution, theta_distribution, savedSteps)
    for (n = 0; n < local_n; n++)
    {
        _swap = 0;
        // Restarting the initial matrix
        if (omp_get_thread_num() == 0)
        {
            for (k = 0; k < K; k++)
                for (i = 0; i < N; i++)
                    for (j = 0; j < M; j++)
                        u[N*M*k + N*i + j] = Ta;
        }
        incerteza_angulo = theta_distribution(generator);
        incerteza_magnitude = magnitude_distribution(generator);
        // EDP
        for (k = 0; k < STEPS; k++)
        {
        #   pragma omp for
            for (i = 0; i < TAM; i++)
            {
                for (j = 0; j < TAM; j++)
                {
                    pele = tecido[i][j];
                    x = i;
                    y = j;
                    vx = velocidade_x(v[N*i + j].magnitude, v[N*i + j].angulo, incerteza_angulo, incerteza_magnitude);
                    vy = velocidade_y(v[N*i + j].magnitude, v[N*i + j].angulo, incerteza_angulo, incerteza_magnitude);

                    x == 0         ? lessX = x : lessX = x - 1;
                    x == (TAM - 1) ? plusX = x : plusX = x + 1;

                    x == 0         ? kim = kappa[pele] : kim = 2*kappa[tecido[lessX][y]]*kappa[pele]/(kappa[tecido[lessX][y]] + kappa[pele]);
                    x == (TAM - 1) ? kip = kappa[pele] : kip = 2*kappa[tecido[plusX][y]]*kappa[pele]/(kappa[tecido[plusX][y]] + kappa[pele]);

                    y == 0         ? lessY = y : lessY = y - 1;
                    y == (TAM - 1) ? plusY = y : plusY = y + 1;

                    y == 0         ? kjm = kappa[pele] : kjm = 2*kappa[tecido[x][lessY]]*kappa[pele]/(kappa[tecido[x][lessY]] + kappa[pele]);
                    y == (TAM - 1) ? kjp = kappa[pele] : kjp = 2*kappa[tecido[x][plusY]]*kappa[pele]/(kappa[pele] + kappa[tecido[x][plusY]]);

                    vx > 0         ? upwind_i = e[pele]*p_b*c_b*vx*(u[N*M*_swap + N*x + y] - u[N*M*_swap + N*lessX + y])/hx : upwind_i = e[pele]*p_b*c_b*vx*(u[N*M*_swap + N*plusX + y] - u[N*M*_swap + N*x + y])/hx;
                    vy > 0         ? upwind_j = e[pele]*p_b*c_b*vy*(u[N*M*_swap + N*x + y] - u[N*M*_swap + N*x + lessY])/hx : upwind_j = e[pele]*p_b*c_b*vy*(u[N*M*_swap  + N*x + plusY] - u[N*M*_swap + N*x + y])/hx;

                      u[N*M*(!_swap) + N*x + y] = (ht/(p*c*(1 - e[pele]) + p_b*c_b*e[pele]))*(
                      (kip*u[N*M*_swap + N*plusX + y] + kjp*u[N*M*_swap + N*x + plusY] - (kip+kim+kjp+kjm)*u[N*M*_swap + N*x + y] + kim*u[N*M*_swap + N*lessX + y] + kjm*u[N*M*_swap + N*x + lessY])/(pow(hx, 2))
                      + qm[pele]*(1 - e[pele])
                      + Calor_nano[x][y]
                      - upwind_i
                      - upwind_j)
                      + u[N*M*_swap + N*x + y];
                }
            }
            _swap = !_swap;

            if (k % INTERVAL == 0 && omp_get_thread_num() == 0)
                for (i = 0; i < TAM; i++)
                    for (j = 0; j < TAM; j++)
                        u_new[n*N*M*savedSteps + N*M*(k/INTERVAL) + N*i + j] = u[N*M*_swap + N*i + j];
        }
    }
    //  END EXPERIMENTS LOOP

    //  Calculating SDV and MEAN
    for (k = 0; k < savedSteps; k++)
    {
        for (i = 0; i < TAM; i++)
        {
            for (j = 0; j < TAM; j++)
            {
                aux = 0;
                for (n = 0; n < local_n; n++)
                    aux += u_new[n*N*M*savedSteps + N*M*k + N*i + j];

                aux = aux/local_n;
                mean[M*N*k + N*i + j] = aux;
            }

        }
    }

    for (k = 0; k < savedSteps; k++)
    {
        for (i = 0; i < TAM; i++)
        {
            for (j = 0; j < TAM; j++)
            {
                aux = 0;
                for (n = 0; n < local_n; n++)
                    aux += pow((u_new[n*N*M*savedSteps + N*M*k + N*i + j] - mean[M*N*k + N*i + j]), 2);

                aux = sqrt(aux/local_n);
                sdv[M*N*k + N*i + j] = aux;
            }
        }
    }

    // Now we need to send the local vectors of 'sdv' and 'mean' to the master process -> Reduce
    // MPI_Reduce (void* input_data_p, void* output_data_p, int count, MPI_Datatype datatype, MPI_Op operator, int dest_process, MPI_Comm comm)
    // if count > 1 => MPI_Reduce can operate with arrays instead of scalars
    MPI_Reduce (mean, final_mean, savedSteps*M*N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce (sdv, final_sdv, savedSteps*M*N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        for (k = 0; k < savedSteps; k++)
        {
            sprintf(str, "./resultados/hipertermia_heterogeneo_np_%d.txt", k+1);
            FILE *saida = fopen(str, "w");
            if (saida == NULL)
            {
                printf("Ocorreu um erro na criação do arquivo %s.\n", str);
                return 1;
            }
            for (i = 0; i < TAM; i++)
            {
                for (j = 0; j < TAM; j++)
                    fprintf(saida, "%lf ", final_mean[M*N*k + N*i + j]/comm_size);

                fprintf(saida, "\n");
            }
            fclose(saida);
            saveVTK(final_mean, k, comm_size);
        }
        for (k = 0; k < savedSteps; k++)
        {
            sprintf(str, "./resultados/hipertermia_sdv_%d.txt", k+1);
            FILE *saida = fopen(str, "w");
            for (i = 0; i < TAM; i++)
            {
                for (j = 0; j < TAM; j++)
                    fprintf(saida, "%lf ", final_sdv[M*N*k + N*i + j]/comm_size);

                fprintf(saida, "\n");
            }
            fclose(saida);
        }
    }

    delete[] u;
    delete[] u_new;
    delete[] sdv;
    delete[] mean;
    MPI_Finalize();
    //  Freeing memory
    if (my_rank == 0)
    {
        for (i = 0; i < TAM; i++)
        {
            delete[] tecido[i];
            delete[] Calor_nano[i];
            delete[] initial_vel[i];
        }
        delete[] final_mean;
        delete[] final_sdv;
        delete[] v;
        delete[] initial_vel;
        delete[] tecido;
        delete[] qm;
        delete[] e;
        delete[]Calor_nano;
        delete[] kappa;
        printf("Programa executado com sucesso.\n");
    }


    return 0;
}
