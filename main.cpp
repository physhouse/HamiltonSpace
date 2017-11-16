#include "HamiltonSpace.h"
#include <mpi.h>
#include <memory>

using namespace Hamilton_Space;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    auto hs = std::make_shared<HamiltonSpace>();

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Wrapping Up..\n");
    MPI_Finalize();
    return 0;
}
