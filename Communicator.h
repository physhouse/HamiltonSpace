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

class Communicator
{
public:
    Communicator();
    ~Communicator();
    void setup(HS_Float, Atom &);
    void exchangeAtoms(Atom &);
    void communicateGhosts(Atom &);
    void reverseCommunicateGhosts(Atom &);

private:
    int rank;
    int numSwaps;				// Total number of data communication swap
    std::vector<std::vector<int> > neighbors(3, std::vector<int>(2));	// The Map to store neighbor Proc IDs, this is used for exchange only

    // These two vectors are used for generating ghost atom list only
    std::vector<int> sendToProc;		// For each communication, the IDs of proc I will send my atoms to
    std::vector<int> recvFromProc;		// For each communication, the IDs of proc I will receive atoms from
    std::vector<int> pbcFlags;			// Whether any PBC is enforced on this swap
						// It is a bit flag
    

    // Buffer arrays for the exchange
    HS_Float* bufferSend;
    HS_Float* bufferRecv;
    
};

}
#endif

