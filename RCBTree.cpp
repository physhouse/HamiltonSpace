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
    printf("[RCBTree] RCBTree Load Balancer\n");
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create MPI data type
    MPI_Type_contiguous(sizeof(MiddleCut), MPI_CHAR, &cut_type);
    MPI_Type_commit(&cut_type);
    MPI_Op_create(cutMerge, 1, &cut_op);
     
    atom = p;

    //tree.resize(nprocs);  // Tree records the RCB information for each processors
    rcbinfo.cut = HS_INFINITY; // rcbinfo records the dimension of the current cut, as well as the cut position
    rcbinfo.dimension = -1; // Default: The rcbinfo have dimension -1 and infinity for cut, meaning not a valid infomation

    allocateMatrix2D(particles, BUFFMAX, 3);
    allocateMatrix2D(velocities, BUFFMAX, 3);
    allocateArray1D(mark, BUFFMAX);
    allocateArray1D(lowerList, BUFFMAX);
    allocateArray1D(upperList, BUFFMAX);
    allocateArray1D(bufferSend, BUFFMAX);
    allocateArray1D(bufferRecv, BUFFMAX);
}

RCBTree::~RCBTree()
{
    destroy(particles);
    destroy(velocities);
    destroy(mark);
    destroy(lowerList);
    destroy(upperList);
    destroy(bufferSend);
    destroy(bufferRecv);
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
    upperBound[0] = atom->box.length[0];
    upperBound[1] = atom->box.length[1];
    upperBound[2] = atom->box.length[2];

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
        //printf("[Rank %d] local %d global %d\n", rank, numParticles, totalNumParticles);
        int targetLower = static_cast<int> (totalNumParticles * (procMiddle - procLow) / (procHigh - procLow + 1));
        int targetUpper = totalNumParticles - targetLower;

        // Find the longest axis, cut according to that dimension        
        int cutDim = findCutDimension(lowerBound, upperBound);

        // Note: cutLow denotes the lowerbound of the slab that I would like to examine
        //       cutHigh denotes the upperbound of the slab that I would like to examine
        //       By assumption, cutMiddle = 1/2 * (cutLow + cutHigh)
        cutLow = lowerBound[cutDim];  
        cutHigh = upperBound[cutDim];
        // The cut Middle point is assumed to be determined 
        // based on the assumtion that procs owns same number of atoms
        //cutMiddle = (procMiddle - procLow) * (upperBound[cutDim] - lowerBound[cutDim]) / (procHigh - procLow + 1);

        // The MiddleCut information
        MiddleCut midme, mid;
        bool breakFlag = false;
    
        // Approach 1: Determining the Middle CutPlain through recursive search
        while (1)
        {
            cutMiddle = 0.5 * (cutLow + cutHigh);
            // PHASE I: Generate the MiddleCut Information
            midme.maxLower = -HS_INFINITY;
            midme.minUpper = HS_INFINITY;
            midme.countLower = 0;
            midme.countUpper = 0;
            midme.totalLower = 0;
            midme.totalUpper = 0;
         
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
                cutLow = mid.minUpper;
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
                        if ((mid.totalLower + countBefore < targetLower))
                        {
                            int numMove = std::min(targetLower - countBefore, localCount);
                            for (int i=0; i<numMove; i++)
                            {
                                mark[upperList[i]] = 0;
                            }
                        }
                        breakFlag = true;
                    }
                }
            }
            else if (mid.totalUpper < targetUpper) // Indicating the upper part need to be expanded
            {
                cutHigh = mid.maxLower;
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
                    if (mid.totalUpper + mid.countLower < targetUpper) // All boundary particles need to be sacrificed
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
                            int numMove = std::min(targetUpper - countBefore, localCount);
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
            printf("[Rank %d] cutMiddle = %lf totalUpper %d totalLower %d\n", rank, cutMiddle, mid.totalUpper, mid.totalLower);
            if (breakFlag) break;
        } // End of while(1) loop

        // Store the cut infomation only if I am the procMiddle
        if (rank == procMiddle)
        {
           rcbinfo.cut = cutMiddle; 
           rcbinfo.dimension = cutDim;
        }

        // Mark the particles to 0 (lower than cut) and 1 (higher than cut)
        // This partitioning is used for communication
        //
        // Find out the processor I would like to exchange particles with
        int myBuddy, activeSide;
        if (rank < procMiddle)
        {
             myBuddy = rank + procMiddle - procLow;
             activeSide = 0; 	// Lower than the cutoff -> proc 0->procMiddle-1
             upperBound[cutDim] = cutMiddle;
        }
        else
        {
             myBuddy = rank - procMiddle + procLow;
             activeSide = 1;	// Higher than the cutoff -> proc procMiddle->procHigh-1
             lowerBound[cutDim] = cutMiddle;
        }

        // Partition the list based on the particle marks
        // particles[0->nkeep-1] = particles remains here, particles[nkeep->numParticles] = particles leaving this processor
        int nkeep = numParticles;
        int i = 0;

        MPI_Barrier(MPI_COMM_WORLD);
        printf("[%d] buddy %d side %d nkeep %d lo %d hi %d\n", rank, myBuddy, activeSide, nkeep, procLow, procHigh);
        MPI_Barrier(MPI_COMM_WORLD);

        while (i < nkeep)
        {
            if (mark[i] != activeSide)
            {
                swap(i, nkeep-1);
                nkeep--;
            }
            else i++;
        }


        int sendCount = 0;
        for (int i=nkeep; i<numParticles; i++)
        {
            bufferSend[6*sendCount] = particles[i][0];
            bufferSend[6*sendCount + 1] = particles[i][1];
            bufferSend[6*sendCount + 2] = particles[i][2];
            bufferSend[6*sendCount + 3] = velocities[i][0];
            bufferSend[6*sendCount + 4] = velocities[i][1];
            bufferSend[6*sendCount + 5] = velocities[i][2];
            sendCount++;
        }
 
        // Start the communication process
        MPI_Status status;
        MPI_Request request;
        int recvCount;
        MPI_Irecv(&recvCount, 1, MPI_INT, myBuddy, 0, MPI_COMM_WORLD, &request);
        MPI_Send(&sendCount, 1, MPI_INT, myBuddy, 0, MPI_COMM_WORLD);
        MPI_Wait(&request, &status);
        printf("rank %d send %d recv %d\n", rank, sendCount, recvCount);

        MPI_Irecv(bufferRecv, 6*recvCount, MPI_DOUBLE, myBuddy, 0, MPI_COMM_WORLD, &request);
        MPI_Send(bufferSend, 6*sendCount, MPI_DOUBLE, myBuddy, 0, MPI_COMM_WORLD);
        MPI_Wait(&request, &status);

        numParticles = nkeep + recvCount;
        for (int i=0; i<recvCount; i++)
        {
            particles[i + nkeep][0] = bufferRecv[6*i];
            particles[i + nkeep][1] = bufferRecv[6*i + 1];
            particles[i + nkeep][2] = bufferRecv[6*i + 2];
            velocities[i + nkeep][0] = bufferRecv[6*i + 3];
            velocities[i + nkeep][1] = bufferRecv[6*i + 4];
            velocities[i + nkeep][2] = bufferRecv[6*i + 5];
        }

        // Split the communicator to 1/2 of its original size
        int split;
        if (rank < procMiddle)
        {
            procHigh = procMiddle - 1;
            split = 0;
        }       
        else
        {
            procLow = procMiddle;
            split = 1;
        }

        MPI_Comm comm_half;
        MPI_Comm_split(comm, split, rank, &comm_half);
        MPI_Comm_free(&comm);
        comm = comm_half;
        
    } // End of while (proc!=procLow && proc!=procHigh) loop
    printFrame();
    MPI_Barrier(comm);
    printf("[Rank %d] lower %lf %lf %lf upper %lf %lf %lf\n", rank, lowerBound[0],lowerBound[1],lowerBound[2],upperBound[0],upperBound[1],upperBound[2]);

    // Update the domain infomation to Atom class
    rcbinfo.split[0] = atom->box.range[0][0] = lowerBound[0];
    rcbinfo.split[1] = atom->box.range[1][0] = lowerBound[1];
    rcbinfo.split[2] = atom->box.range[2][0] = lowerBound[2];
    rcbinfo.split[3] = atom->box.range[0][1] = upperBound[0];
    rcbinfo.split[4] = atom->box.range[1][1] = upperBound[1];
    rcbinfo.split[5] = atom->box.range[2][1] = upperBound[2];

    // Update the particle coordinate information
    atom->nlocal = numParticles;
    atom->nghost = 0;
    atom->nall = atom->nlocal + atom->nghost;
    for (int i=0; i<numParticles; i++)
    {
        atom->x[i][0] = particles[i][0];
        atom->x[i][1] = particles[i][1];
        atom->x[i][2] = particles[i][2];
        atom->v[i][0] = velocities[i][0];
        atom->v[i][1] = velocities[i][1];
        atom->v[i][2] = velocities[i][2];
    }
}

// Return the RCB tree node that the current processor is keeping
RCBTreeNode RCBTree::getMyRCBTreeNode()
{
    return rcbinfo;
}

// helper functions
void RCBTree::swap(int i, int j)
{
    double temp;
    int tempMark;
    for (int dim = 0; dim < 3; dim++) 
    {
        temp = particles[i][dim];
        particles[i][dim] = particles[j][dim];
        particles[j][dim] = temp;
        temp = velocities[i][dim];
        velocities[i][dim] = velocities[j][dim];
        velocities[j][dim] = temp;
    }
    tempMark = mark[i];
    mark[i] = mark[j];
    mark[j] = temp;
}

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
// Stragety:
//     maxLower coutLower procLower -> should be dealed with together
//     minUpper etc. similar
//
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

// For Debugging purpose
void RCBTree::printFrame()
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    char filename[100];
    sprintf(filename, "rcb_p%d.dat", me);
    FILE* fp = fopen(filename, "w");
    fprintf(fp, "%s\nnparticles %d\n", filename,  numParticles);
    for (int i=0; i<numParticles; i++)
    {
        fprintf(fp, "%lf %lf %lf %d\n", particles[i][0], particles[i][1], particles[i][2], mark[i]);
    }
    fclose(fp);
}
