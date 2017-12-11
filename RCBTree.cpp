#include "RCBTree.h"
#include "Atom.h"
#include "Type.h"
#include "Memory.h"
#include <algorithm>
#include <mpi.h>

using namespace Hamilton_Space;

void cutMerge(void *, void *, int *, MPI_Datatype *);

RCBTree::RCBTree(std::shared_ptr<class Atom> p)
{
    // Determining the total size of MPI communication spaces
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create MPI data type
    MPI_Type_contiguous(sizeof(MiddleCut), MPI_CHAR, &cut_type);
    MPI_Type_commit(&cut_type);
    MPI_Op_create(cutMerge, 1, &cut_op);
     
    atom = p;

    tree.resize(nprocs);  // Tree records the RCB information for each processors
    rcbinfo = new RCBTreeNode; // rcbinfo records the dimension of the current cut, as well as the cut position

    allocateMatrix2D(particles, BUFFMAX, 3);
    allocateArray1D(mark, BUFFMAX);
    allocateArray1D(lowerList, BUFFMAX);
    allocateArray1D(upperList, BUFFMAX);
}

RCBTree::~RCBTree()
{
    destroy(particles);
    destroy(mark);
    destroy(lowerList);
    destroy(upperList);
}

void RCBTree::buildDistributedRCBTree()
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

    numParticles = atom->nlocal;
    for (int i=0; i<numParticles; i++)
    {
        particles[i][0] = atom->x[i][0];
        particles[i][1] = atom->x[i][1];
        particles[i][2] = atom->x[i][2];
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

        // total number of particles
        int totalNumParticles = 0;
        MPI_Allreduce(&numParticles, &totalNumParticles, 1, MPI_INT, MPI_SUM, comm);
        int targetLower = static_cast<int> (totalNumParticles * (procMiddle - procLow) / (procHigh - procLow + 1));
        int targetUpper = totalNumParticles - targetLower;

        // Find the longest axis, cut according to that dimension        
        int cutDim = findCutDimension(lowerBound, upperBound);
        cutLow = lowerBound[cutDim];
        cutHigh = upperBound[cutDim];
        // The cut Middle point is assumed to be determined 
        // based on the assumtion that procs owns same number of atoms
        cutMiddle = (procMiddle - procLow) * (upperBound[cutDim] - lowerBound[cutDim]) / (procHigh - procLow + 1);

        // The MiddleCut information
        MiddleCut midme, mid;
        bool breakFlag = false;
    
        // Approach 1: Gradual Change of Middle CutPlain\
        // TODO: Recursive determination of middle CutPlain
        while (1)
        {
            // PHASE I: Generate the MiddleCut Information
            midme.maxLower = -HS_INFINITY;
            midme.minUpper = HS_INFINITY;
            midme.countLower = 0;
            midme.countUpper = 0;
         
            int indexLower, indexUpper;
            for (int i=0; i<numParticles; i++)
            {
                HS_float x = particles[i][cutDim];
                if (x < cutMiddle)
                {
                    mark[i] = 0; // Marker that marks particle in the lower half
                    midme.totalLower++;
                    if (closeEnough(x, midme.maxLower))
                    {
                        lowerList[midme.countLower] = i;
                        midme.countLower++;
                    }
                    else if (x > midme.maxLower)
                    {
                        midme.maxLower = x;
                        midme.countLower = 1;
                        lowerList[0] = i;
                    }
                }
                else
                {
                    mark[i] = 1; // Marker that marks particle in the lower half
                    midme.totalUpper++;
                    if (closeEnough(x, midme.minUpper))
                    {
                        lowerList[midme.countUpper] = i;
                        midme.countUpper++;
                    }
                    else if (x < midme.minUpper)
                    {
                        midme.minUpper = x;
                        midme.countUpper = 1;
                        upperList[0] = i;
                    }
                } 
            }
            // Gather the MiddleCut information from all the processors
            MPI_Allreduce(&midme, &mid, 1, cut_type, cut_op, comm); // TODO: Implement OP_CUT

            // PHASE II: Load Balancing
            // Basic Idea: for the unbalanced partition, sacrifice the boundary layer(or a single particle) of the fatter partition
            // mark them as leaving particles and adjust the cut plane to be the position of those particles
            // iterate until equipartition
            if (mid.totalLower < targetLower)  // Indicating the lower part need to be expanded
            {
                cutMiddle = mid.minUpper;
                if (mid.countUpper == 1)  // Only move one particle
                {
                    if (mid.totalLower + mid.countUpper < targetLower)
                    {
                        if (rank == midme.procUpper)
                        {
                            mark[upperList[0]] = 0; // In the upper Ranks, mark the attribution of moving particles to be 0
                        }
                    }
                    else breakFlag = true;
                }
                else  // Moving multiple particles
                {
                    if (mid.totalLower + mid.countUpper < targetLower) // All boundary particles need to be sacrificed
                    {
                       for (int i=0; i<midme.countUpper; i++)
                       {
                           mark[upperList[i]] = 0; // Move all to the lower partition
                       } 
                        
                    }
                    else  // Enough for balance, continue
                    {
                        // Do a scan over all the ranks, determine how many to move
                        int localCount = 0;
                        int countBefore;
                        if (closeEnough(midme.minUpper, mid.minUpper)) localCount = midme.countUpper;
                        MPI_Scan(&localCount, &countBefore, 1, MPI_INT, MPI_SUM, comm);
                        if (mid.totalLower + countBefore < targetLower)
                        {
                            int numMove = std::min(targetLower - countBefore, midme.countUpper);
                            for (int i=0; i<numMove; i++)
                            {
                                mark[upperList[i]] = 0;
                            }
                        }
                        breakFlag = true;
                    }
                }
            }
            else if (mid.totalLower < targetLower) // Indicating the upper part need to be expanded
            {
                cutMiddle = mid.maxLower;
                if (mid.countLower == 1)  // Only move one particle
                {
                    if (mid.totalUpper + mid.countLower < targetUpper)
                    {
                        if (rank == midme.procLower)
                        {
                            // In the lower Ranks, mark the attribution of moving particles to be 1
                            // So that it can be migrated to upper partition
                            mark[lowerList[0]] = 1; 
                        }
                    }
                    else breakFlag = true;
                }
                else  // Moving multiple particles
                {
                    if (mid.totalUpper + mid.countLower < targetLower) // All boundary particles need to be sacrificed
                    {
                       for (int i=0; i<midme.countLower; i++)
                       {
                           mark[lowerList[i]] = 0; // Move all to the lower partition
                       } 
                        
                    }
                    else  // Enough for balance, continue
                    {
                        // Do a scan over all the ranks, determine how many to move
                        int localCount = 0;
                        int countBefore;
                        if (closeEnough(midme.maxLower, mid.maxLower)) localCount = midme.countLower;
                        MPI_Scan(&localCount, &countBefore, 1, MPI_INT, MPI_SUM, comm);
                        if (mid.totalUpper + countBefore < targetUpper)
                        {
                            int numMove = std::min(targetUpper - countBefore, midme.countLower);
                            for (int i=0; i<numMove; i++)
                            {
                                mark[lowerList[i]] = 0;
                            }
                        }
                        breakFlag = true;
                    }
                }
            }
            else breakFlag = true; // Lucky: Even Partition Achieved
  
            if (breakFlag) break;
        }

        // Mark the particles to 0 (lower than cut) and 1 (higher than cut)
        // This partitioning is used for communication
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

// Merge cut structure of two inputs
void cutMerge(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
    RCBTree::MiddleCut *cut1 = (RCBTree::MiddleCut *) in;
    RCBTree::MiddleCut *cut2 = (RCBTree::MiddleCut *) inout;

    cut2->totalLower += cut1->totalLower;
    if (cut2->maxLower <cut1->maxLower)
    {
        cut2->maxLower = cut1->maxLower;
        cut2->countLower = cut1->countLower;
        cut2->procLower = cut1->procLower;
    }
    else if (closeEnough(cut2->maxLower, cut1->maxLower))
    {
        cut2->countLower += cut1->countLower;
        if (cut1->procLower < cut2->procLower) cut2->procLower = cut1->procLower;
    }

    cut2->totalUpper += cut1->totalUpper;
    if (cut2->minUpper > cut1->minUpper)
    {
        cut2->minUpper = cut1->minUpper;
        cut2->countUpper = cut1->countUpper;
        cut2->procUpper = cut1->procUpper;
    }
    else if (closeEnough(cut2->minUpper, cut1->minUpper))
    {
        cut2->countUpper += cut1->countUpper;
        if (cut1->procUpper < cut2->procUpper) cut2->procUpper = cut1->procUpper;
    }
}
