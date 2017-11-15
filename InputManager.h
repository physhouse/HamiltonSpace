#ifndef _INPUTMANAGER_H_
#define _INPUTMANAGER_H_

namespace HamiltonSpace {

#include <Type.h>

struct InputManager
{
    HS_float boxLength;
    HS_float beta;
    int      nsteps;
    HS_float dt;
    int      atomPerProc;
    HS_float mass;
}

}

#endif
