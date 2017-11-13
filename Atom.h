#ifndef _HAMILTONSPACE_H_
#define _HAMILTONSPACE_H_

#include <memory>
#include <vector>
#include <Atom.h>
#include <Type.h>

namespace HamiltonSpace {

struct Box
{
    std::vector<std::vector<HS_float> > range;
    HS_float lengthx;
    HS_float lengthy;
    HS_float lengthz;
};

class Atom
{
public:
    Atom();
    ~Atom();
    void packSendAtoms(int first, int last, int dim, HS_float lo, HS_float hi, int* count, HS_float* buffer);
    void unpackRecvAtoms(int count, HS_float* buffer);
    void packExchange(HS_float* buffer, int count, int dimDirectionIndex);
    void unpackExchange(int count, HS_float* buffer);
    void swap(int i, int j);

public:
    int nall;		// total number of atoms
    int nlocal, nghost; // nlocal->number of host atoms, nghost->number of ghost atoms
    HS_float mass;	// For Lennard-Jones Particle with uniform mass

    Box box;

    std::vector<std::vector<HS_float> > x;
    std::vector<std::vector<HS_float> > v;
    std::vector<std::vector<HS_float> > f;
};

}

#endif
