#include "HamiltonSpace.h"
#include <mpi.h>
#include <memory>

using namespace Hamilton_Space;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    auto hs = new HamiltonSpace();
    delete hs;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
