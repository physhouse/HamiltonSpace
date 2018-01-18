#ifndef _COMMRCB_H_
#define _COMMRCB_H_

// This class implements the RCB tree type Domain Decomposition Parallelization
// General Molecular Dynamics and Monte Carlo Simulations
// 

#include <vector>
#include <memory>
#include "Atom.h"
#include "Type.h"
#include "InputManager.h"
#include "RCBTree.h"

// TODO: Design the Inheritance Hierarchy for Communicators
// 01/04/18
//
namespace Hamilton_Space {

// This data structure keeps the information of swap for each dimension for a specific direction
struct RCBSwap
{
    int dim; 			// Dimension of this swap'
    int numOverlaps;		// Number of overlapping domains I need to send particles to
    int numRecvOverlaps;		// Number of overlapping domains I need to recv particles from

    std::vector<int> sendToProc;	// For each communication, the IDs of proc I will send my atoms to
    std::vector<int> recvFromProc;	// For each communication, the IDs of proc I will receive atoms from
    std::vector<int> sendNum;
    std::vector<int> recvNum;
    std::vector<int> pbcFlags;		// Whether any PBC is enforced on this swap
                                        // It is a bit flag
    std::vector<HS_float*> slablo;
    std::vector<HS_float*> slabhi;

    std::vector<HS_float*> bufferSend;
    std::vector<HS_float*> bufferRecv;

    std::vector<int*> sendList;
    std::vector<int>  recvStartIndex;
};

struct RCBExchange
{
    int numExchange;
    std::vector<int> sendToProc;	// For each communication, the IDs of proc I will send my atoms to
    std::vector<int> recvFromProc;	// For each communication, the IDs of proc I will receive atoms from
    std::vector<int> sendNum;
    std::vector<int> recvNum;
    std::vector<HS_float*> bufferSend;
    std::vector<HS_float*> bufferRecv;
};

class CommRCB
{
public:
    CommRCB(std::shared_ptr<InputManager> input);
    ~CommRCB();
    void clear();
    void setup(HS_float, std::shared_ptr<class Atom>, RCBTreeNode);
    void generateGhosts(std::shared_ptr<class Atom>);			// Generate the ghost lists after building neighborlist
    void exchangeAtoms(std::shared_ptr<class Atom>);			// Exchange Atoms to neighbors after building neighborlist
    void communicateGhosts(std::shared_ptr<class Atom>);		// communicate ghost list information each step
    void reverseCommunicateGhosts(std::shared_ptr<class Atom>);	// communicate ghost atom forces each step if using Newton 3rd

    // Walk the tree
    void setupTree(RCBTreeNode);
    void findBoxOverlapProcs(HS_float*, HS_float*, std::vector<int>& overlap, int start, int end);
    void boxIntersection(HS_float* ghostLower, HS_float* ghostUpper, int targetBox, HS_float* slablo, HS_float* slabhi);
    bool boxTouch(int, int, int);

private:
    int rank;
    int nprocs;
    int numSwaps;				// Total number of ghost atoms communication swap
    std::vector<RCBSwap> swap;			// Swap data structure
    std::vector<RCBExchange> exchange;		// Exchange data structure
    RCBTreeNode* tree;				// The RCB binary tree for recording the cut and dimension information
    std::vector<HS_float> boxRange;
};

}
#endif
