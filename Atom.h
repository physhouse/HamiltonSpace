#ifndef _HAMILTONSPACE_H_
#define _HAMILTONSPACE_H_

#include <memory>
#include <vector>
#include <Atom.h>
#include <Type.h>

namespace HamiltonSpace {

class Atom
{
public:
    Atom();
    ~Atom();
    void packSendAtoms(int first, int last, HS_float lo, HS_float hi, int* count, std::vector<HS_float> &buffer);

public:
    int natoms;		// total number of atoms
    int nlocal, nghost; // nlocal->number of host atoms, nghost->number of ghost atoms
    HS_float mass;	// For Lennard-Jones Particle with uniform mass

    std::vector<std::vector<HS_float> > x;
    std::vector<std::vector<HS_float> > v;
    std::vector<std::vector<HS_float> > f;
};

}

#endif
