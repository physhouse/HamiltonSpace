#ifndef _RCBTREE_H_
#define _RCBTREE_H_

#include "Atom.h"
#include "Type.h"
#include <mpi.h>

namespace Hamilton_Space {

class RCBTree
{
public:
    RCBTree(std::shared_ptr<class Atom>);
    ~RCBTree();
    void buildDistributedRCBTree();
    int  findCutDimension(HS_float*, HS_float*);
    void swap(int, int);

    void printFrame();
    
    RCBTreeNode getMyRCBTreeNode();

    // Datastructure for the cut information
    struct MiddleCut {
        int totalLower, totalUpper;    // Total number of particles in the lower/upper division
        HS_float maxLower, minUpper;   // The maximum coordinates of particles in lower division and minimum in upper
        int countLower, countUpper;    // Number of particles at that position
        int procLower, procUpper;      // which processor owns the particle
    };

private:
    int nprocs;		// MPI Communication 
    int rank;		// MPI Communication
    
    HS_float** particles;  // Store the particle coordinates information
    HS_float** velocities; // Store the particle velocities information
    int* mark;
    int* lowerList;
    int* upperList;

    MPI_Datatype cut_type;
    MPI_Op cut_op;

    int numParticles;
    std::shared_ptr<class Atom> atom;
    RCBTreeNode rcbinfo;

    HS_float* bufferSend;
    HS_float* bufferRecv;
};

}

#endif
