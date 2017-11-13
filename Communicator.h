#ifndef _COMMUNICATOR_H_
#define _COMMUNICATOR_H_

// This class implements the Grid-type Domain Decomposition Parallelization of General Molecular Dynamics and Monte Carlo Simulations
// 
// 2017-11-06

#include <map>
#include <vector>
#include <Atom.h>
#include <string>

namespace HamiltonSpace {

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
};

class Communicator
{
public:
    Communicator();
    ~Communicator();
    void setup(HS_Float, std::shared_ptr<Atom>);
    void generateGhosts(std::shared_ptr<Atom>);			// Generate the ghost lists after building neighborlist
    void exchangeAtoms(std::shared_ptr<Atom>);			// Exchange Atoms to neighbors after building neighborlist
    void communicateGhosts(std::shared_ptr<Atom>;		// communicate ghost list information each step
    void reverseCommunicateGhosts(std::shared_ptr<Atom>);	// communicate ghost atom forces each step if using Newton 3rd

private:
    int rank;
    int numSwaps;				// Total number of ghost atoms communication swap
    std::vector<Swap> swap;  			// Swap data structure

    // The Map to store neighbor Proc IDs, this is used for exchange only
    std::vector<std::vector<int> > neighbors(3, std::vector<int>(2));	

    // The buffer for atom exchange
    HS_float* bufferExchangeSend;
    HS_float* bufferExchangeRecv;
};

}
#endif

