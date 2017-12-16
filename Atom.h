#ifndef _ATOM_H_
#define _ATOM_H_

#include <memory>
#include <vector>
#include "Type.h"
#include "InputManager.h"

namespace Hamilton_Space {

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
    Atom(std::shared_ptr<InputManager> input);
    ~Atom();
    void genInitialConfig(std::shared_ptr<InputManager> input);
    void clearGhost();
   
    void packSendAtoms(int first, int last, int dim, HS_float lo, HS_float hi, int pbcFlag, int* count, HS_float* buffer, int* sendList);
    void unpackRecvAtoms(int count, HS_float* buffer);
    void packExchange(HS_float* buffer, int* count, int dimDirectionIndex);
    void unpackExchange(int count, HS_float* buffer);
    void packCommunicateGhosts(int count, int dim, int pbcFlag, HS_float* buffer, int* sendlist);
    void unpackCommunicateGhosts(int count, int startIndex, HS_float* buffer);
    void packReverseCommunication(int count, int startIndex, HS_float* buffer);
    void unpackReverseCommunication(int count, int* sendlist, HS_float* buffer);

    void swap(int i, int j);
 
    void enforcePBC();
    void printFrame();  // For Test
    void propagate();   // For Test

public:
    int nall;		// total number of atoms
    int nlocal, nghost; // nlocal->number of host atoms, nghost->number of ghost atoms
    int natoms;		// total number of atoms in the system
    HS_float mass;	// For Lennard-Jones Particle with uniform mass

    Box box;

    std::vector<std::vector<HS_float> > x;
    std::vector<std::vector<HS_float> > v;
    std::vector<std::vector<HS_float> > f;
};

}

#endif
