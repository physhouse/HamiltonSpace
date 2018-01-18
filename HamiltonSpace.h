#ifndef _HAMILTONSPACE_H_
#define _HAMILTONSPACE_H_

#include <memory>
#include "Atom.h"
#include "Communicator.h"
#include "InputManager.h"
#include "CommRCB.h"
#include "RCBTree.h"

namespace Hamilton_Space {

class HamiltonSpace
{
public:
    HamiltonSpace();
    ~HamiltonSpace();
    void parseInput();
    void exec();

private:
    std::shared_ptr<InputManager> input;
    std::shared_ptr<Atom> atom;
    //std::shared_ptr<class NeighborList> neighbor;
    std::shared_ptr<Communicator> messenger;
    std::shared_ptr<CommRCB> messengerRCB;
    std::shared_ptr<RCBTree> balancer;
    //std::shared_ptr<class Integrator> manager;
};

}

#endif
