#include "RCBTree.h"
#include <mpi.h>

using namespace HamiltonSpace;

RCBTree::RCBTree(std::shared_ptr<Atom> p) : atom(p)
{
    // Determining the total size of MPI communication spaces
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
    tree.resize(nprocs);  // Tree records the RCB information for each processors
    rcbinfo = new RCBTreeNode; // rcbinfo records the dimension of the current cut, as well as the cut position
}

void RCBTree:buildDistributedRCBTree()
{
    // Recursively construct the RCB Tree
    // for each iteraction, split the processors into two groups and generate a new communication subgroup
    // The lower subgroup maintains information below the cut, wheras the upper one maintains information above the cut

    // The root starts from the whole simulation box
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    int procLow = 0;
    int procHigh = nprocs - 1;

    HS_float lowerBound[3];
    HS_float upperBound[3];

    // Root: The whole simulation box
    lowerBound[0] = 0.0;
    lowerBound[1] = 0.0;
    lowerBound[2] = 0.0;
    upperBound[0] = atom->box.lengthx;
    upperBound[1] = atom->box.lengthy;
    upperBound[2] = atom->box.lengthz;

    //  cutLow             cutMiddle          cutHigh
    //      _____________________________________
    //     |                  :                  |
    //     |                  :                  |
    //     |__________________:__________________|
    //         proclo<->mid     mid<->prochigh
    //

    HS_float cutLow, cutMiddle, cutHigh;

    for (int i=0; i<atom->nlocal; i++)
    {
        particles[i][0] = x[i][0];
        particles[i][1] = x[i][1];
        particles[i][2] = x[i][2];
    }

    // Recursively build the binary tree
    // In this simple version, the load balance is ensured such that each node has nearly equal number of atoms
    // The attempted cut is made by linear interpolation
    while (procLow != procHigh)
    {
        // Procedure for building binary trees
        // 1. determine how many processors are there in each communication subgroup
        // 2. Determine which dimension of the 3 to cut, based on the heuristics that the longest axis should be cut
        // 3. attempt to cut the space according to boxlenght * procMiddle / nprocs
        
        int procMiddle = procLow + (procHigh - procLow) / 2 + 1; 
        // For odd number of processors give the lower half one more processor

        // Find the longest axis, cut according to that dimension        
        int cutDim = findCutDimension(lowerBound, upperBound);
        cutLow = lowerBound[cutDim];
        cutHigh = upperBound[cutDim];
        // The cut Middle point is assumed to be determined 
        // based on the assumtion that procs owns same number of atoms
        cutMiddle = (procMiddle - procLow) * (upperBound[cutDim] - lowerBound[cutDim]) / (procHigh - procLow + 1);

        // The MiddleCut information
        MiddleCut midme, mid;
    
        while (1)
        {
            // PHASE I: Generate the MiddleCut Information
            midme.maxLower = -VERY_BIG_FLOAT;
            midme.minUpper = VERY_BIG_FLOAT;
            midme.countLower = 0;
            midme.countUpper = 0;
         
            int indexLower, indexUpper;
            for (int i=0; i<numParticles; i++)
            {
                HS_float x = particle[i][cutDim];
                if (x < cutMiddle)
                {
                    mark[i] = 0; // Marker that marks particle in the lower half
                    midme.totalLower++;
                    if (closeEnough(x, midme.maxLower))
                    {
                        midme.countLower++;
                    }
                    else if (x > midme.maxLower)
                    {
                        midme.maxLower = x;
                        midme.countLower = 1;
                        indexLower = i;
                    }
                }
                else
                {
                    mark[i] = 1; // Marker that marks particle in the lower half
                    midme.totalUpper++;
                    if (closeEnough(x, midme.minUpper))
                    {
                        midme.countUpper++;
                    }
                    else if (x < midme.minUpper)
                    {
                        midme.minUpper = x;
                        midme.countUpper = 1;
                        indexUpper = i; 
                    }
                } 
            }
            // Gather the MiddleCut information from all the processors
            MPI_Allreduce(&midme, &mid, 1, TYPE_CUT, OP_CUT, comm); // TODO: Implement OP_CUT

            // PHASE II: Load Balancing
            // Basic Idea: for the unbalanced partition, sacrifice the boundary layer(or a single particle) of the fatter partition
            // mark them as leaving particles and adjust the cut plane to be the position of those particles
            // iterate until equipartition
            if (mid.totalLower < targetLower)  // Indicating the lower part need to be expanded
            {
                cutMiddle = mid.minUpper;
                if (mid.countLower == 1)  // Only move one particle
                {
                    if (mid.totalLower + mid.countLower < targetLower)
                    {
                        if (rank == midme.procUpper)
                        {
                            mark[indexUpper] = 0; // In the upper Ranks, mark the attribution of moving particles to be 0
                        }
                    }
                    else break;
                }
                else  // Moving multiple particles
                {
                    if (mid.totalLower + mid.countLower < targetLower)
                    {
                        // TODO: 12/05 Yining
                    }
                }
            }
            else if (mid.totalUpper < targetUpper) // Indicating the upper part need to be expanded
            {
            }
            else break; // Lucky: Even Partition Achieved
        }

        // Mark the particles to 0 (lower than cut) and 1 (higher than cut)
        // This partitioning is used for communication
        int numUpperPartLocal = 0;
        for (int i=0; i<numParticles; i++)
        {
            //TODO  12/01/Yining
	    if (particles[i][cutDim] < cutMiddle) mark[i] = 0;
            else
	    {
		mark[i] = 1;
		numUpperPartLocal++;
            }
        }
        // Get the number of Atoms in the upper part
        MPI_Allreduce(&numUpperPart, &numUpper, 1, MPI_INT, MPI_SUM, comm);
    }
}

// helper functions
int RCBTree::findCutDimension(HS_float* lowerBound, HS_float* upperBound)
{
    int cutDim = 0;
    HS_float cutLength = upperBound[0] - lowerBound[0];
    if ((upperBound[1] - lowerBound[1]) > cutLength) 
    {
        cutDim = 1;
        cutLength = upperBound[1] - lowerBound[1];
    }
    if ((upperBound[2] - lowerBound[2]) > cutLength) 
    {
        cutDim = 2;
        cutLength = upperBound[2] - lowerBound[2];
    }
    return cutDim;
}

bool RCBTree:closeEnough(double x, double y)
{
    return ; //TODO: Float Point Number comparison
}
