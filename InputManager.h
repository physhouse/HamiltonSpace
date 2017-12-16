#ifndef _INPUTMANAGER_H_
#define _INPUTMANAGER_H_

#include "Type.h"

namespace Hamilton_Space {

struct InputManager
{
    HS_float boxLength;
    HS_float beta;
    int      nsteps;
    HS_float dt;
    int      numAtoms;
    HS_float mass;
    HS_float cutNeighbor;
};

}

#endif
