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
    atom->printFrame();
    
    balancer = std::make_shared<RCBTree>(atom);
    balancer->buildDistributedRCBTree();
    printf("finished setup\n");

    messengerRCB = std::make_shared<CommRCB>(input);
    messengerRCB->setup(input->cutNeighbor, atom, balancer->getMyRCBTreeNode());
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
        GET_INT(input->numAtoms)
        GET_REAL(input->boxLength)
        GET_REAL(input->beta)
        GET_REAL(input->dt)
        GET_REAL(input->mass)
        GET_REAL(input->cutNeighbor)
    VAR_END
}

void HamiltonSpace::exec()
{
    // Testing the Grid Communicator
    //messenger->generateGhosts(atom);
    //for (int i=0; i<100; i++)
    //{
    //   atom->propagate();
    //   messenger->communicateGhosts(atom);
    //   messenger->reverseCommunicateGhosts(atom);
    //}
    //messenger->exchangeAtoms(atom);
    //atom->enforcePBC();
    //messenger->generateGhosts(atom);
    //atom->printFrame();

    messengerRCB->generateGhosts(atom);
    atom->printFrame();
}
