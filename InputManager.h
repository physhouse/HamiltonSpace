#ifndef _INPUTMANAGER_H_
#define _INPUTMANAGER_H_

namespace HamiltonSpace {

#include <Type.h>

class InputManager
{
public:
    HS_float boxLength;
    HS_float beta;
    int      nsteps;
    HS_float dt;
}

}

#endif
