#include "HamiltonSpace.h"
#include "Parser.h"
#include <cstdio>

using namespace Hamilton_Space;

HamiltonSpace::HamiltonSpace()
{
    input = std::make_shared<InputManager>();
    parseInput();

    atom = std::make_shared<Atom>(input);

    messenger = std::make_shared<Communicator>(input);
    messenger->setup(input->cutNeighbor, atom);

    atom->genInitialConfig(input);
    exec();
}

HamiltonSpace::~HamiltonSpace()
{
    printf("[HamiltonSpace] HamiltonSpace object get destructed\n");
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
        GET_REAL(input->cutNeighbor)
    VAR_END
}

void HamiltonSpace::exec()
{
    messenger->generateGhosts(atom);
    for (int i=0; i<100; i++)
    {
       atom->propagate();
       messenger->communicateGhosts(atom);
       messenger->reverseCommunicateGhosts(atom);
    }
    messenger->exchangeAtoms(atom);
    atom->enforcePBC();
    messenger->generateGhosts(atom);
    atom->printFrame();
}
