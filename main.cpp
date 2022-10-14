#include <iostream>
#include <mpi.h>

const double EPS = 1E-9;

int main(int argc, char **argv) {
    int errCode;

    if ((errCode = MPI_Init(&argc, &argv)) != MPI_SUCCESS) {
        std::cerr << "MPI_Init error: code " << errCode << std::endl;
        MPI_Abort(MPI_COMM_WORLD, errCode);
    }
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Hello, MPI world! I am " << rank << " of " << size << std::endl;
    MPI_Finalize();
    return 0;
}
