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

private:
    int nprocs;		// MPI Communication 
    int rank;		// MPI Communication
    
    std::shared_ptr<Atom> atom;
    RCBTreeNode* rcbinfo;
    std::vector<RCBTreeNode> tree;
}

}

#endif
