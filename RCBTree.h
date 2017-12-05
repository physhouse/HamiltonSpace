#ifndef _RCBTREE_H_
#define _RCBTREE_H_

#include "Atom.h"
#include "Type.h"

namespace HamiltonSpace {

struct RCBTreeNode
{
    int dimension;
    HS_float cut;
}

class RCBTree
{
public:
    RCBTree(std::shared_ptr<Atom>);
    ~RCBTree();
    void buildDistributedRCBTree();
    int  findCutDimension(HS_float*, HS_float*);

    // Datastructure for the cut information
    struct MiddleCut {
        int totalLower, totalUpper;    // Total number of particles in the lower/upper division
        HS_float maxLower, minUpper;   // The maximum coordinates of particles in lower division and minimum in upper
        int countLower, countUpper;    // Number of particles at that position
        int procLower, procUpper;      // which processor owns the particle
    }

private:
    int nprocs;		// MPI Communication 
    int rank;		// MPI Communication
    
    std::shared_ptr<Atom> atom;
    RCBTreeNode* rcbinfo;
    std::vector<RCBTreeNode> tree;
}

}

#endif
