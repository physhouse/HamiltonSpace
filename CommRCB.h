#ifndef _COMMRCB_H_
#define _COMMRCB_H_

// This class implements the RCB tree type Domain Decomposition Parallelization
// General Molecular Dynamics and Monte Carlo Simulations
// 
// 2017-12-18

#include <vector>
#include <memory>
#include "Atom.h"
#include "Type.h"
#include "InputManager.h"
#include "RCBTree.h"

namespace Hamilton_Space {

struct Swap
{
    int sendToProc;		// For each communication, the IDs of proc I will send my atoms to
    int recvFromProc;		 // For each communication, the IDs of proc I will receive atoms from
    int sendNum;
    int recvNum;
    int pbcFlags;		// Whether any PBC is enforced on this swap
                                                // It is a bit flag
    int dim; 			// Dimension of this swap
    HS_float slablo;
    HS_float slabhi;

    HS_float* bufferSend;
    HS_float* bufferRecv;

    int* sendList;
    int recvStartIndex;
};

class CommRCB
{
public:
    CommRCB(std::shared_ptr<InputManager> input);
    ~Communicator();
    void setup(HS_float, std::shared_ptr<class Atom>);
    void generateGhosts(std::shared_ptr<class Atom>);			// Generate the ghost lists after building neighborlist
    void exchangeAtoms(std::shared_ptr<class Atom>);			// Exchange Atoms to neighbors after building neighborlist
    void communicateGhosts(std::shared_ptr<class Atom>);		// communicate ghost list information each step
    void reverseCommunicateGhosts(std::shared_ptr<class Atom>);	// communicate ghost atom forces each step if using Newton 3rd

private:
    std::vector<int> processorGrid;
    int rank;
    int numSwaps;				// Total number of ghost atoms communication swap
    std::vector<Swap> swap;  			// Swap data structure

    // The Map to store neighbor Proc IDs, this is used for exchange only
    std::vector<std::vector<int> > neighbor;	

    // The buffer for atom exchange
    HS_float* bufferExchangeSend;
    HS_float* bufferExchangeRecv;
};

}
#endif

