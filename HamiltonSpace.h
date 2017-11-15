#ifndef _HAMILTONSPACE_H_
#define _HAMILTONSPACE_H_

#include <memory>
#include <Atom.h>

namespace HamiltonSpace {

class HamiltonSpace
{
public:
    HamiltonSpace(int argc, char** argv);
    ~HamiltonSpace();
    parseInput();
    exec();

private:
    std::shared_ptr<class InputManager> input;
    std::shared_ptr<class Atom> atom;
    //std::shared_ptr<class NeighborList> neighbor;
    std::shared_ptr<class Communication> messenger;
    std::shared_ptr<class Integrator> manager;
};

}

#endif
