#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TERMINATOR -1

void fseq(int n) {
    int c = 0;

    int* prime = (int*)malloc(n * n * sizeof(int));

    for(int i=2; i<=n*n; i++){
        prime[i] = 1;
    }
    // double sqrt_n = sqrt(n);
    for (int i = 2; c<n && i < n*n; i++){
    
        if (prime[i]==1)
        {
            for (int j = i+i; j <= n*n; j= j+i)
            {
                prime[j] = 0;
            }

            c++;
        }
    }
}


int main(int argc, char** argv) {
    int rank, size;
    float t1, t2, t3, t4;
    int n, current_prime, nb;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Enter the nb of primes to find: \n");
        scanf("%d", &n);

        t3 = MPI_Wtime();
    	fseq(n);
    	t4 = MPI_Wtime();
        
        t1 = MPI_Wtime();
        int start = 3;
        current_prime = 2;
        printf("Process %d found the prime: %d\n", rank, current_prime);

        for (int nb = start; nb <= n * n; nb += 2) {
            if (nb % current_prime != 0) {
                MPI_Send(&nb, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            }
        }

        nb = TERMINATOR;
        MPI_Send(&nb, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);

    } else {
        int first_nb = 1;

        while (1) {
            MPI_Recv(&nb, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            if (nb == TERMINATOR) {
                //not last p
                if (rank < size - 1) {
                    MPI_Send(&nb, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                }
                break;
            }

            if (first_nb) {
                current_prime = nb;
                printf("Process %d: Prime found: %d\n", rank, current_prime);
                first_nb = 0;
            } else {
                if (nb % current_prime != 0) {
                    if (rank < size - 1) {
                        MPI_Send(&nb, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                    }
                }
            }
        }
    }

    t2 = MPI_Wtime();

    if (rank == 0) {
        float tp = t2-t1;
        float tseq = t4-t3;
        float Sp = tseq/tp;
        float E = Sp/size;
        printf("Sequential time: %f seconds\n", tseq);
        printf("Parallel Execution time: %f seconds\n", tp);
        printf("Speedup: %f\n", Sp);
        printf("Efficiency: %f\n", E);
    }

    MPI_Finalize();
    return 0;
}
