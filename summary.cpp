#include <mpi.h>
#include <iostream>

using namespace std;

void summary(int N, int my_rank, int commsize) {
    double count = 0;
    if (N > commsize) {
        double end = (my_rank + 1) * N / commsize - 1;
        double k = (my_rank)*N / commsize;
        
        if (my_rank == 0)
            k = 1;
        if (my_rank == commsize - 1)
            end = N;
        for (; k <= end; k++)
            count += 1 / k;
        
    }
    else {
        if (my_rank < N)
            count = 1 / ((double)my_rank + 1);
    }
    MPI_Send(&count, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    
    int my_rank, commsize;
    MPI_Init(&argc, &argv);
    MPI_Status* status = new MPI_Status;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int N = atoi(argv[1]);
    
    if (my_rank > 0)
        summary(N, my_rank, commsize);
    if (my_rank == 0) {
        double val = 0;
        if (N > commsize) {
            double end = (my_rank + 1) * N / commsize - 1;
            double k = 1;
            if (my_rank == commsize - 1)
                end = N;
            for (; k <= end; k++)
                val += 1 / k;
        }
        else
            val = 1;

        double sum = val;
        val = 0;
        for (int i = 1; i < commsize; i++) {
            MPI_Recv(&val, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, status);
            sum += val;
        }
        cout << sum << endl;
    }
    MPI_Finalize();
    if (my_rank == 0) {
        double sum = 0;
        for (int i = 1; i <= N; i++)
            sum += 1 / (double)i;
        cout << sum;
    }
}