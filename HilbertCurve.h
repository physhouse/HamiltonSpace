/* This subsection includes the code for generating the hillbert 
 * space-filling curve load balancing data structures
 * 
 * Brief Introduction of how the code works
 * 
 * Note that: The Grid-type spatial decomposition is always performed first
 * 	      This is also true for RCBTree load balancer
 *
 * Step0: The simulation domain is divided into N (preferably power of 2) cells
 * 	  Each cell is mapped into a 1D coordinate defined by the HSF curve
 * Step1: The communicator should accumulate the number of particles in each cell
 * 	  And then the master will determine how to divide the HSF curve so that
 * 	  each partition has approximately equal amount of particles
 * Step2: A key function Map(), which maps the particles of the current processors
 * 	  to have markers either denotes inside the local domain or should be migrated
 * Step3: Perform the particle migration
 * Step4: Find out the neighboring cells that could provide ghost particles to the current
 *        sub-domain, by implementing the FindInterestingCells() and FindInterestingRanks()
 *        Then generate the ghost particles lists to perform communication with at each 
 *        time step
 *
 *  @Author: Yining Han (ynhan@uchicago.edu)
 *  https://physhouse.github.com
 *
 */

#ifndef _HILBERT_CURVE_H_
#define _HILBERT_CURVE_H_

namespace Hamilton_Space {

// The data structure for the individual hilbert cell
struct HillbertCell
{
    GridIndex		gid;
    HilbertCellIndex    id;
    HilbertWork		workload;
};

// The section of HSFC that the current processor owns
struct HilbertSection
{
    int owner;
    HilbertCellIndex begin;
    HilbertCellIndex end;   
};

// The server for the main program to get HSFC information
class HilbertCurve
{
public:
    HilbertCurve();
    ~HilbertCurve();
    void SetupHilbertCurve();
    void GenerateHilbertSection(); // Load Balancer to distribute workloads

    inline void CellIndex2HilbertIndex(int3 gridPoint, HilbertIndex &index);     

    // The following routines are used for trading particles to their correct position
    void MapParticles();    // Map the particles to mark their markers, either as local or to be migrated

    // The following routines are used for ghost list generation
    // For Hook Up
    void FindInterestingCells();  // Find the IDs of potentially interested cells
    void FindInterestingRanks();  // Find the procs that I need to talk to when building ghost lists

private:
    struct CellInfo
    {
        int  owner;
	int  countHookups;
        int* Hookups;
        HilbertCell cell;
    };
    std::vector<HilbertCells> cells;
    HilbertSection section;

};

}

#endif
