#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int my_rank, commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int val = 2;
    if (my_rank != 0) {
        MPI_Recv(&val, 1, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "I'm proccess " << my_rank << " got number " << val << " from proccess " << my_rank - 1 << endl;
        val *= val;
    }

    MPI_Send(&val, 1, MPI_INT, (my_rank + 1) % commsize, 0, MPI_COMM_WORLD);
    if (my_rank == 0){
        MPI_Recv(&val, 1, MPI_INT, commsize - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "I'm proccess " << my_rank << " got number " << val << " from proccess " << commsize - 1;
    }
    MPI_Finalize();
}