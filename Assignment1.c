#include <mpi.h>
#include <stdio.h>
#include <math.h>

int steps = 0; //comp steps

// Function to evaluate the curve (y = f(x))
float f(float x) {
    steps++;
    return x * x ; // Example: y = x^2
}

// Function to compute the area of a trapezoid
float trapezoid_area(float a, float b, float d) { 
    float area = 0;
    for (float x = a; x < b; x+=d) {
        area += f(x) + f(x+d);
    }
    
    return area * d / 2.0f;
}

//specialized for seq so the comp steps don't add up
float f2(float x) {
    return x * x ; // Example: y = x^2
}

float trapSeq(float a, float b, int n) {
    float h = (b - a) / n;
    float area = 0.0;

    // Calculate the area
    for (int i = 0; i < n; i++) {
        area += (f2(a + i * h) + f2(a + (i + 1) * h)) * h / 2.0;
    }

    return area;
}

int main(int argc, char** argv) {
    int rank, size;
    float a = 0.0f, b = 1.0f;  // Limits of integration
    int n;
    float start, end, local_area, total_area;
    float t1, t2, t3, t4; //for mpiTime
    
    MPI_Init(&argc, &argv); // Initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get number of processes

    if (rank == 0) {
        // Get the number of intervals from the user
        printf("Enter the number of intervals: ");
        for (int i = 1; i <= 10; i++) {}
        scanf("%d", &n);
    	t3 = MPI_Wtime();
    	float seqtrap = trapSeq(a,b,n);
    	t4 = MPI_Wtime();
    }
    
    // Broadcast the number of intervals to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Calculate the interval size for each process
    float d = (b - a) / n; // delta
    float region = (b - a)/ size;
    
    // Calculate local bounds for each process
    start = a + rank * region;
    end = start + region;
    
    //start
    t1 = MPI_Wtime();
    
    // Each process calculates the area of its subinterval
    local_area = trapezoid_area(start, end, d);
    
    // Reduce all local areas to the total area on the root process
    MPI_Reduce(&local_area, &total_area, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //end
    t2 = MPI_Wtime();

    int Tsteps = 0;//total
    MPI_Reduce(&steps, &Tsteps, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("The total area under the curve is: %f\n", total_area);
        float texec = t2-t1;
        float tseq = t4-t3;
        float Sp = tseq/texec;
        float E = Sp/size;
        printf("Sequential time: %f seconds\n", tseq);
        printf("Parallel Execution time: %f seconds\n", texec);
        printf("Speedup: %f\n", Sp);
        printf("Efficiency: %f\n", E);
        printf("Total computational steps: %d\n", Tsteps);
        
    }
    
    MPI_Finalize(); // Finalize MPI
    return 0;
}
