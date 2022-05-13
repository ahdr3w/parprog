//!!!! output: 3.0378e-07 sec



#include <mpi.h>
#include <iostream>
#include <time.h>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int my_rank, commsize;
    double t1, t2;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int buf1 = 1;
    int buf2 = 0;
    t1 = MPI_Wtime();
    for (int i = 0; i < 5000; i++) {
        if (my_rank == 0)
            MPI_Send(&buf1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        if (my_rank == 1)
            MPI_Recv(&buf2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (my_rank == 0) {
        t2 = MPI_Wtime();
        cout << "time : " << (t2 - t1) / 5000 << " sec";
    }
    MPI_Finalize();
}