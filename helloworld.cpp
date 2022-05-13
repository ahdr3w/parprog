#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    int commsize, my_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    cout << "Hello world rank: " << my_rank << ", commsize: " << commsize << endl;
    MPI_Finalize();
}
