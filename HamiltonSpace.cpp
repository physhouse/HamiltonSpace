#include <HamiltonSpace.h>
#include <Parser.h>

using namespace HamiltonSpace {

HamiltonSpace::HamiltonSpace()
{
    input = std::make_share();
}

void HamiltonSpace::parseInput()
{
    VAR_BEGIN
        GET_INT(input->nsteps)
        GET_INT(input->atomPerProc)
        GET_REAL(input->boxLength)
        GET_REAL(input->beta)
        GET_REAL(input->dt)
        GET_REAL(input->mass)
    VAR_END
}


}
