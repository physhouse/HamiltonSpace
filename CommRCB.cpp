// This class implements the Recursive Binary Tree Domain Decomposition Parallelization of General Particle Based 
// Molecular Dynamics and Monte Carlo Simulations
// It can also be extended in Mesh based simulation techniques

#include "CommRCB.h"
#include <mpi.h>
#include <algorithm>
#include <cstring>
#include <cassert>

using namespace Hamilton_Space;

CommRCB::CommRCB(std::shared_ptr<InputManager> input)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    printf("rank = %d nprocs = %d\n", rank, nprocs);
    tree = new RCBTreeNode [nprocs]; 

    boxRange.resize(3);
    swap.resize(6);
    exchange.resize(3);
    printf("[CommRCB] Communicator of RCB Decomposition Constructed\n");
}

CommRCB::~CommRCB()
{
    clear();
    delete[] tree;
    printf("[CommRCB] CommRCB object get desctructed\n");
}

void CommRCB::clear()
{
    for (const auto& item : swap)
    {
        for (int i=0; i<item.numOverlaps; i++)
        {
            delete[] item.bufferSend[i];
            delete[] item.sendList[i];
        }
        for (int i=0; i<item.numRecvOverlaps; i++)
        {
            delete[] item.bufferRecv[i];
        }
    }
}

void CommRCB::setup(HS_float cutoff, std::shared_ptr<Atom> atom, RCBTreeNode rcbone)
{
    for (int dim=0; dim<3; dim++) boxRange[dim] = atom->box.length[dim];

    setupTree(rcbone);
    printf("[CommRCB] RCB Tree Constructed\n");
    printf("[%d] cut = %lf dim = %d\n", rank, tree[rank].cut, tree[rank].dimension);
    // Determining which boxes we need to communicate with to generate the ghost particle list
    // For each dimension, look at both up and down positions, construct the ghost box 
    // Find the processors that overlap with the ghost box, add them to the overlap list
    // If the current processor itself lies in the list, take care of this special case
    // one do not need to worry about the "corner" blocks, the reason is similar to the grid like tessellation

    bool oneBox, twoBox;
    HS_float ghostLower[3], ghostUpper[3];
    HS_float ghostSecondLower[3], ghostSecondUpper[3];

    for (int idim = 0; idim < 3; idim++)
    {
        for (int dir = 0; dir < 2; dir++)
        {
           ghostLower[0] = atom->box.range[0][0];
           ghostLower[1] = atom->box.range[1][0];
           ghostLower[2] = atom->box.range[2][0];
           ghostUpper[0] = atom->box.range[0][1];
           ghostUpper[1] = atom->box.range[1][1];
           ghostUpper[2] = atom->box.range[2][1];
    
           // dir == 0 refers to send to lower part (recv from upper part)
           if (dir == 0) 
           {
               ghostLower[idim] -= cutoff;
               ghostUpper[idim] = atom->box.range[idim][0];
           }
           else
           {
               ghostUpper[idim] += cutoff; 
               ghostLower[idim] = atom->box.range[idim][1];
           }
  
           oneBox = true;
           twoBox = false;
           // Tell if wrap over the peoriodic boundaries
           if (dir == 0 && ghostLower[idim] < 0.0)
           {
               twoBox = true;
               std::memcpy(ghostSecondLower, ghostLower, sizeof(ghostSecondLower));
               std::memcpy(ghostSecondUpper, ghostUpper, sizeof(ghostSecondUpper));
               ghostSecondLower[idim] += atom->box.length[idim];
               ghostSecondUpper[idim] = atom->box.length[idim];
               if (closeEnough(ghostUpper[idim], 0.0))
               {
                   oneBox = false;
               }
           }
           else if (dir == 1 && ghostUpper[idim] > atom->box.length[idim])
           {
               twoBox = true;
               std::memcpy(ghostSecondLower, ghostLower, sizeof(ghostSecondLower));
               std::memcpy(ghostSecondUpper, ghostUpper, sizeof(ghostSecondUpper));
               ghostSecondLower[idim] = 0.0;
               ghostSecondUpper[idim] -= atom->box.length[idim];
               if (closeEnough(ghostLower[idim], atom->box.length[idim]))
               {
                   oneBox = false;
               }
           }

           if (oneBox)
           {
               ghostLower[idim] = std::max(ghostLower[idim], 0.0);
               ghostUpper[idim] = std::min(ghostUpper[idim], atom->box.length[idim]);
           }
           printf("[%d] one %s two %s ghostLower[%d] = %lf %lf %lf  ghostUpper[%d] = %lf %lf %lfghost2Lower[%d] = %lf ghost2Upper[%d] = %lf\n", rank, oneBox?"true":"false", twoBox?"true":"false", idim, ghostLower[0], ghostLower[1], ghostLower[2], idim, ghostUpper[0], ghostUpper[1], ghostUpper[2], idim, ghostSecondLower[idim], idim, ghostSecondUpper[idim]);

           std::vector<int> overlapBoxOne;
           std::vector<int> overlapBoxTwo;
           printf("[%d] start %d end %d\n", rank, 0, nprocs-1);
           // Traverse the tree overlapping the first box 
           if (oneBox)  findBoxOverlapProcs(ghostLower, ghostUpper, overlapBoxOne, 0, nprocs - 1);
           // Deal with the second box, deal with pbc issues right here
           if (twoBox)  findBoxOverlapProcs(ghostSecondLower, ghostSecondUpper, overlapBoxTwo, 0, nprocs - 1);

           int iSwap = 2*idim + dir;
           swap[iSwap].dim = idim;
           swap[iSwap].numOverlaps = 0;
           if (oneBox) swap[iSwap].numOverlaps += overlapBoxOne.size();
           if (twoBox) swap[iSwap].numOverlaps += overlapBoxTwo.size();

           printf("[%d] Dimension %d Direction %d numOverlaps = %d\n", rank, idim, dir, swap[iSwap].numOverlaps);

           int totalOverlap = swap[iSwap].numOverlaps;
           swap[iSwap].sendToProc.resize(totalOverlap);

           if (dir == 0) 
           {
	       swap[iSwap+1].numRecvOverlaps = totalOverlap;
	       swap[iSwap+1].recvFromProc.resize(totalOverlap);
               swap[iSwap+1].recvNum.resize(totalOverlap);
               swap[iSwap+1].bufferRecv.resize(totalOverlap);
               swap[iSwap+1].recvStartIndex.resize(totalOverlap);
           }
           else
	   {
	       swap[iSwap-1].numRecvOverlaps = totalOverlap;
	       swap[iSwap-1].recvFromProc.resize(totalOverlap);
               swap[iSwap-1].recvNum.resize(totalOverlap);
               swap[iSwap-1].bufferRecv.resize(totalOverlap);
               swap[iSwap-1].recvStartIndex.resize(totalOverlap);
           }

           swap[iSwap].sendNum.resize(totalOverlap);
           swap[iSwap].pbcFlags.resize(totalOverlap);
           swap[iSwap].slablo.resize(totalOverlap);
           swap[iSwap].slabhi.resize(totalOverlap);

           swap[iSwap].bufferSend.resize(totalOverlap);
           swap[iSwap].sendList.resize(totalOverlap);
          
           int nOverlap = 0;
           MPI_Barrier(MPI_COMM_WORLD);
           if (oneBox)
           {
               for (auto &item : overlapBoxOne)
               {
                   printf("[%d] boxone overlap with %d\n", rank, item);
                   // Notice that in the first case overlap cannot be rank by definition
                   // Generate the swap data structure and push it back to the swap vector
                   // increment nSwap to record the swap vector size explicitly
                   swap[iSwap].sendToProc[nOverlap] = item;
                   swap[iSwap].bufferSend[nOverlap] = new HS_float [BUFFMAX];

                   if (dir == 0) swap[iSwap+1].bufferRecv[nOverlap] = new HS_float [BUFFMAX];
                   else swap[iSwap-1].bufferRecv[nOverlap] = new HS_float [BUFFMAX];

                   swap[iSwap].sendList[nOverlap] = new int [BUFFMAX];

                   swap[iSwap].slablo[nOverlap] = new HS_float [3];
                   swap[iSwap].slabhi[nOverlap] = new HS_float [3];
                   printf("[%d] iSwap = %d here\n", rank, iSwap);

                   // Reasoning about this two line of code:
                   // Note that dir == 0 refers to the situation that you send particles to lower part
                   // In the next swap which represent the same dimension but with dir == 1 we should receive particels
                   // from the processor we send to at this round to balance the communication
                   // Similar analysis hold for the case which dir == 1
                   if (dir == 0) swap[iSwap + 1].recvFromProc[nOverlap] = item;
                   else swap[iSwap - 1].recvFromProc[nOverlap] = item;
                   // No pbc wrap is happening this round
                   swap[iSwap].pbcFlags[nOverlap] = 0;
 
                   // Determine slablo slabhi for this round of send
                   // This is quite different from the normal uniform tessellation
                   // The first step relies on find the intersection of the proposed ghost box and
                   // the actual RCB Box you will sent particles to
                   
                   boxIntersection(ghostLower, ghostUpper, item, swap[iSwap].slablo[nOverlap], swap[iSwap].slabhi[nOverlap]);
                   // Determining the sending box
                   if (dir == 0)
                   {
                       swap[iSwap].slablo[nOverlap][idim] = atom->box.range[idim][0];
                       swap[iSwap].slabhi[nOverlap][idim] = std::min(atom->box.range[idim][1], tree[item].split[3+idim] + cutoff);
                   }
                   else
                   {
                       swap[iSwap].slablo[nOverlap][idim] = std::max(atom->box.range[idim][0], tree[item].split[idim] - cutoff);
                       swap[iSwap].slabhi[nOverlap][idim] = atom->box.range[idim][1];
                   }

                   // To take the corners into account Very important information
                   // If send in the x axis, just send particles in the sendbox
                   // If send in the y axis, you need to send the ghost particles received in the previous communication, which was accepted from x axis trade
                   // to your neighbors in the y axis, thus we need to enlarge the sending box to take these corner areas into account explicitly
                   // similarly, sending in z axis, we need to include the ghost particles in box y and x axis communications
                   if (idim >= 1)
                   {
                       swap[iSwap].slablo[nOverlap][0] -= cutoff;
                       swap[iSwap].slabhi[nOverlap][0] += cutoff;
                       printf("[debug %d] modified lo %lf hi %lf\n", rank, swap[iSwap].slablo[nOverlap][0], swap[iSwap].slabhi[nOverlap][0]);
                   }
                   if (idim == 2)
                   {
                       swap[iSwap].slablo[nOverlap][1] -= cutoff;
                       swap[iSwap].slabhi[nOverlap][1] += cutoff;
                       printf("[debug %d] modified lo %lf hi %lf\n", rank, swap[iSwap].slablo[nOverlap][1], swap[iSwap].slabhi[nOverlap][1]);
                   }

                   printf("[%d] dim %d  slablo %lf %lf %lf slabhi %lf %lf %lf\n", rank, idim, swap[iSwap].slablo[nOverlap][0],swap[iSwap].slablo[nOverlap][1],swap[iSwap].slablo[nOverlap][2],swap[iSwap].slabhi[nOverlap][0],swap[iSwap].slabhi[nOverlap][1],swap[iSwap].slablo[nOverlap][2]);
	           nOverlap++;
               }
           }

           int nOverlap2 = 0;
           if (twoBox)
           {
               for (auto &item : overlapBoxTwo)
               {
                   printf("[%d] boxtwo overlap with %d\n", rank, item);
                   // Note that in the second ghost box, the proc myself could probably sit in the overlapping list
                   // We need to move it to the last of the overlap list
                   if (item == rank) std::swap(overlapBoxTwo[nOverlap2], overlapBoxTwo[overlapBoxTwo.size()-1]);
                   swap[iSwap].sendToProc[nOverlap + nOverlap2] = item;
                   swap[iSwap].bufferSend[nOverlap + nOverlap2] = new HS_float [BUFFMAX];
 
                   if (dir == 0) swap[iSwap+1].bufferRecv[nOverlap + nOverlap2] = new HS_float [BUFFMAX];
                   else swap[iSwap-1].bufferRecv[nOverlap + nOverlap2] = new HS_float [BUFFMAX];

                   swap[iSwap].sendList[nOverlap + nOverlap2] = new int [BUFFMAX];

                   swap[iSwap].slablo[nOverlap + nOverlap2] = new HS_float [3];
                   swap[iSwap].slabhi[nOverlap + nOverlap2] = new HS_float [3];
                   printf("[%d] iSwap = %d here\n", rank, iSwap);

                   if (dir == 0) swap[iSwap + 1].recvFromProc[nOverlap + nOverlap2] = item;
                   else swap[iSwap - 1].recvFromProc[nOverlap + nOverlap2] = item;

                   // pbc wrap is happening this round
                   swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_ANY_FLAG;
                   if (dir == 0)
                   {
                       if (idim == 0) swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_POS_X;
                       if (idim == 1) swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_POS_Y;
                       if (idim == 2) swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_POS_Z;
                   }
                   else
                   {
                       if (idim == 0) swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_NEG_X;
                       if (idim == 1) swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_NEG_Y;
                       if (idim == 2) swap[iSwap].pbcFlags[nOverlap + nOverlap2] |= PBC_NEG_Z;
                   }
                   boxIntersection(ghostSecondLower, ghostSecondUpper, item, swap[iSwap].slablo[nOverlap + nOverlap2], swap[iSwap].slabhi[nOverlap + nOverlap2]);
                   // Determining the sending box
                   if (dir == 0)
                   {
                       swap[iSwap].slablo[nOverlap + nOverlap2][idim] = atom->box.range[idim][0];
                       swap[iSwap].slabhi[nOverlap + nOverlap2][idim] = std::min(atom->box.range[idim][1], tree[item].split[3+idim] - atom->box.length[idim] + cutoff);
                   }
                   else
                   {
                       swap[iSwap].slablo[nOverlap + nOverlap2][idim] = std::max(atom->box.range[idim][0], tree[item].split[idim] + atom->box.length[idim] - cutoff);
                       swap[iSwap].slabhi[nOverlap + nOverlap2][idim] = atom->box.range[idim][1];
                   }

                   if (idim >= 1)
                   {
                       swap[iSwap].slablo[nOverlap + nOverlap2][0] -= cutoff;
                       swap[iSwap].slabhi[nOverlap + nOverlap2][0] += cutoff;
                   }
                   if (idim == 2)
                   {
                       swap[iSwap].slablo[nOverlap + nOverlap2][1] -= cutoff;
                       swap[iSwap].slabhi[nOverlap + nOverlap2][1] += cutoff;
                   }
                   printf("[%d] dim %d  slablo %lf %lf %lf slabhi %lf %lf %lf\n", rank, idim, swap[iSwap].slablo[nOverlap + nOverlap2][0],swap[iSwap].slablo[nOverlap + nOverlap2][1],swap[iSwap].slablo[nOverlap + nOverlap2][2],swap[iSwap].slabhi[nOverlap + nOverlap2][0],swap[iSwap].slabhi[nOverlap + nOverlap2][1],swap[iSwap].slablo[nOverlap + nOverlap2][2]);
	           nOverlap2++;
               }
           }
        }
    }
    
    printf("[CommRCB] Ghost Communication Buffer Allocation finished\n");

    // For the exchange phase
    for (int dim = 0; dim < 3; dim++)
    {
        int numExchange = 0;
        std::vector<int> touchProcs;
        for (auto item : swap[2*dim].sendToProc)
        {
            if (boxTouch(item, dim, 0))
    	    {
                touchProcs.push_back(item);
                numExchange++;
            }
        }
        for (auto item : swap[2*dim+1].sendToProc)
        {
            if (boxTouch(item, dim, 1))
    	    {
                touchProcs.push_back(item);
                numExchange++;
            }
        }

        exchange[dim].numExchange = numExchange;
        exchange[dim].sendToProc.resize(numExchange);
        exchange[dim].recvFromProc.resize(numExchange);
        exchange[dim].bufferSend.resize(numExchange);
        exchange[dim].bufferRecv.resize(numExchange);

        for (int i=0; i<numExchange; i++)
        {
            exchange[dim].sendToProc[i] = touchProcs[i];
            exchange[dim].recvFromProc[i] = touchProcs[i];
            exchange[dim].bufferSend[i] = new HS_float [EXCHANGEMAX];
            exchange[dim].bufferRecv[i] = new HS_float [EXCHANGEMAX];
        }
    }

    // Total number of Swaps need to be performed
    printf("[CommRCB] CommRCB setup finished\n");
         
}

/* Routine to generate the ghost atom lists */
void CommRCB::generateGhosts(std::shared_ptr<Atom> atom)
{
    atom->clearGhost(); // Clear up the ghost atom lists
    printf("[%d] nlocal = %d nghost = %d\n", rank, atom->nlocal, atom->nghost);
    
    int first = 0;		// The first atom need to be examined
    int last = atom->nlocal;	// The last atom need to be examined
    int count;
    int previousEnd = atom->nlocal;

    for (int i=0; i<6; i++)
    {
        MPI_Request request[swap[i].numRecvOverlaps];
        int nsend[swap[i].numOverlaps];
        int nrecv[swap[i].numRecvOverlaps];

        for (int j=0; j<swap[i].numRecvOverlaps; j++)
        {
            MPI_Irecv(&nrecv[j], 1, MPI_INT, swap[i].recvFromProc[j], 0, MPI_COMM_WORLD, &request[j]);
        }
       
        for (int j=0; j<swap[i].numOverlaps; j++)
        {
            atom->packSendAtomsRCB(first, last, swap[i].dim, swap[i].slablo[j], swap[i].slabhi[j], swap[i].pbcFlags[j], (int *)&(swap[i].sendNum[j]), swap[i].bufferSend[j], swap[i].sendList[j]);
            nsend[j] = swap[i].sendNum[j];
            MPI_Send(&nsend[j], 1, MPI_INT, swap[i].sendToProc[j], 0, MPI_COMM_WORLD);
        }
        MPI_Waitall(swap[i].numRecvOverlaps, request, MPI_STATUSES_IGNORE);

        for (int j=0; j<swap[i].numRecvOverlaps; j++)
        {
            swap[i].recvNum[j] = nrecv[j];
            MPI_Irecv(swap[i].bufferRecv[j], swap[i].recvNum[j], MPI_DOUBLE, swap[i].recvFromProc[j], 0, MPI_COMM_WORLD, &request[j]);
        }
        for (int j=0; j<swap[i].numOverlaps; j++)
        {
            MPI_Send(swap[i].bufferSend[j], swap[i].sendNum[j], MPI_DOUBLE, swap[i].sendToProc[j], 0, MPI_COMM_WORLD);
        }
        MPI_Waitall(swap[i].numRecvOverlaps, request, MPI_STATUSES_IGNORE);

        for (int j=0; j<swap[i].numRecvOverlaps; j++)
        {
            atom->unpackRecvAtoms(swap[i].recvNum[j], swap[i].bufferRecv[j]);
            last = last + swap[i].recvNum[j] / 3;
        }
        //for (int j=0; j<totalOverlaps; j++)
        //{
        //    atom->packSendAtomsRCB(first, last, swap[i].dim, swap[i].slablo[j], swap[i].slabhi[j], swap[i].pbcFlags[j], (int *)&(swap[i].sendNum[j]), swap[i].bufferSend[j], swap[i].sendList[j]);
        //    nsend = swap[i].sendNum[j];
        //    printf("DEBUGGING: NLOCAL %d NGHOST %d SENDTO %d RECVFROM %d\n", atom->nlocal, atom->nghost, swap[i].sendToProc[j], swap[i].recvFromProc[j]);

        //    MPI_Irecv(&nrecv, 1, MPI_INT, swap[i].recvFromProc[j], 0, MPI_COMM_WORLD, &request);
        //    MPI_Send(&nsend, 1, MPI_INT, swap[i].sendToProc[j], 0, MPI_COMM_WORLD);
        //    MPI_Wait(&request, &status); 
        //    swap[i].recvNum[j] = nrecv;
        //    printf("[%d] nsend = %d nrecv = %d\n", rank, nsend, nrecv);

        //    printf("[%d]:send %d recv %d lo %lf hi %lf sendTo %d recvFrom %d nlocal %d nghost %d\n", rank, swap[i].sendNum[j], swap[i].recvNum[j], swap[i].slablo[j][swap[i].dim], swap[i].slabhi[j][swap[i].dim], swap[i].sendToProc[j], swap[i].recvFromProc[j], atom->nlocal, atom->nghost);

        //    MPI_Irecv(swap[i].bufferRecv[j], swap[i].recvNum[j], MPI_DOUBLE, swap[i].recvFromProc[j], 0, MPI_COMM_WORLD, &request);
        //    MPI_Send(swap[i].bufferSend[j], swap[i].sendNum[j], MPI_DOUBLE, swap[i].sendToProc[j], 0, MPI_COMM_WORLD);
      
        //    MPI_Wait(&request, &status); 

        //    atom->unpackRecvAtoms(swap[i].recvNum[j], swap[i].bufferRecv[j]);

        //    first = 0;
        //    last = last + swap[i].recvNum[j]/3;
        //}
      
        // Generate the startRecvIndex for forward Communication, Exclusive Scan
        swap[i].recvStartIndex[0] = previousEnd;
        for (int j=1; j<swap[i].numRecvOverlaps; j++)
        {
            swap[i].recvStartIndex[j] = swap[i].recvStartIndex[j-1] + swap[i].recvNum[j-1] / 3;
        }
        previousEnd = swap[i].recvStartIndex[swap[i].numRecvOverlaps - 1] + swap[i].recvNum[swap[i].numRecvOverlaps - 1] / 3;

        //if (i == 2) 
        //{
        //    atom->printFrame();
        //    MPI_Barrier(MPI_COMM_WORLD);
        //    exit(-1);
        //}
    }
} 

/* Routine to migrate atoms to neighboring proc */
/* This function is only called when rebuilding the neighbor lists together with generateGhosts() */
void CommRCB::exchangeAtoms(std::shared_ptr<Atom> atom)
{
    atom->clearGhost();

    MPI_Status status;
    int nsend, nrecv;
   
    for (int idim = 0; idim < 3; idim++)
    {
        printf("Dealing with Dim %d\n", idim+1);
        atom->packExchangeRCB(exchange[idim].bufferSend, exchange[idim].sendNum, idim, tree);

        MPI_Request request[exchange[idim].numExchange];
        for (int i=0; i<exchange[idim].numExchange; i++)
        {
            MPI_Irecv(&exchange[idim].recvNum[i], 1, MPI_INT, exchange[idim].recvFromProc[i], 0, MPI_COMM_WORLD, &request[i]);
        }
        for (int i=0; i<exchange[idim].numExchange; i++)
        {
            MPI_Send(&exchange[idim].sendNum[i], 1, MPI_INT, exchange[idim].sendToProc[i], 0, MPI_COMM_WORLD);
        }
        MPI_Waitall(exchange[idim].numExchange, request, MPI_STATUSES_IGNORE);

        for (int i=0; i<exchange[idim].numExchange; i++)
        {
            MPI_Irecv(exchange[idim].bufferRecv[i], exchange[idim].recvNum[i], MPI_DOUBLE, exchange[idim].recvFromProc[i], 0, MPI_COMM_WORLD, &request[i]);
        }
        for (int i=0; i<exchange[idim].numExchange; i++)
        {
            MPI_Send(&exchange[idim].bufferSend[i], exchange[idim].sendNum[i], MPI_DOUBLE, exchange[idim].sendToProc[i], 0, MPI_COMM_WORLD);
        }
        MPI_Waitall(exchange[idim].numExchange, request, MPI_STATUSES_IGNORE);

        atom->unpackExchangeRCB(exchange[idim].bufferRecv, exchange[idim].recvNum, exchange[idim].numExchange);
    }

}

/* Update the ghost particle coordinates for each time step
 * utilizing all the information when initializing the ghost atom lists
 * */
void CommRCB::communicateGhosts(std::shared_ptr<Atom> atom)
{
    MPI_Status status;
  
    for (int i=0; i<6; i++)
    {
        int totalOverlaps = swap[i].numOverlaps;
        MPI_Request request[totalOverlaps];

        /* Note that we need to set up all the recv requests first before actually going into set up the send */
        /* The purpose of this is to prevent dead lock */
        for (int j=0; j<totalOverlaps; j++)
        {
            MPI_Irecv(swap[i].bufferRecv[j], swap[i].recvNum[j], MPI_DOUBLE, swap[i].recvFromProc[j], 0, MPI_COMM_WORLD, &request[i]);
        }    
        for (int j=0; j<totalOverlaps; j++)
        {
	    atom->packCommunicateGhosts(swap[i].sendNum[j]/3, swap[i].dim, swap[i].pbcFlags[j], swap[i].bufferSend[j], swap[i].sendList[j]);
            MPI_Send(swap[i].bufferSend[j], swap[i].sendNum[j], MPI_DOUBLE, swap[i].sendToProc[j], 0, MPI_COMM_WORLD);
        }
        MPI_Waitall(totalOverlaps, request, MPI_STATUS_IGNORE);

        for (int j=0; j<totalOverlaps; j++)
        {
            atom->unpackCommunicateGhosts(swap[i].recvNum[j]/3, swap[i].recvStartIndex[j], swap[i].bufferRecv[j]);
        }
    }
}

/* Backward communication of ghost particle forces, Implementation requires the reversed data structures
 * This is only used when the Newton's third law is applied when evaluating the pairwise forces
 * Otherwise just turned it off */
void CommRCB::reverseCommunicateGhosts(std::shared_ptr<Atom> atom)
{
    MPI_Status status;

    for (int i=0; i<6; i++)
    {
        int totalOverlaps = swap[i].numOverlaps;
        MPI_Request request[totalOverlaps];

        /* Note that we need to set up all the recv requests first before actually going into set up the send */
        /* The purpose of this is to prevent dead lock */
        for (int j=0; j<totalOverlaps; j++)
        {
            MPI_Irecv(swap[i].bufferRecv[j], swap[i].sendNum[j], MPI_DOUBLE, swap[i].recvFromProc[j], 0, MPI_COMM_WORLD, &request[i]);
        }    
        for (int j=0; j<totalOverlaps; j++)
        {
	    atom->packReverseCommunication(swap[i].recvNum[j]/3, swap[i].recvStartIndex[j], swap[i].bufferSend[j]);
            MPI_Send(swap[i].bufferSend[j], swap[i].recvNum[j], MPI_DOUBLE, swap[i].sendToProc[j], 0, MPI_COMM_WORLD);
        }
        MPI_Waitall(totalOverlaps, request, MPI_STATUS_IGNORE);

        for (int j=0; j<totalOverlaps; j++)
        {
            atom->unpackReverseCommunication(swap[i].sendNum[j]/3, swap[i].sendList[j], swap[i].bufferRecv[j]);
        }
    }
}


/* Setup the RCB Tree data structure for determining communication patterns */
void CommRCB::setupTree(RCBTreeNode node)
{
    RCBTreeNode rcbone;
    rcbone.cut = node.cut;
    rcbone.dimension = node.dimension;
    std::memcpy(rcbone.split, node.split, sizeof(node.split));
    MPI_Allgather(&rcbone, sizeof(RCBTreeNode), MPI_CHAR,
 	 	  tree, sizeof(RCBTreeNode), MPI_CHAR, MPI_COMM_WORLD); 
}

/* Find the domains that overlaps with the regtangular box defined by boxlo and boxhi */
/* Recursive Determination 
 * start => the start index of procs that one would like to explore
 * end => the end index of procs that one would like to explore in
 * overlap => vector stroring the overlap procs */
void CommRCB::findBoxOverlapProcs(HS_float* boxlo, HS_float* boxhi, std::vector<int>& overlap, int start, int end)
{
    // printf("[%d] calling start %d end %d\n", rank, start, end);
    // If at the node, add to overlap vector
    if (start == end)
    {
        //printf("[%d] overlaps %d\n", rank, start);
        overlap.push_back(start);
    } 
    else
    {
        int procMiddle = start + (end - start) / 2 + 1;
        HS_float cut = tree[procMiddle].cut;
        int dim = tree[procMiddle].dimension;
        //printf("[%d] dim = %d boxlo[dim] = %lf boxhi[dim] = %lf cut = %lf start = %d end = %d middle = %d\n", rank, dim, boxlo[dim], boxhi[dim], cut, start, end, procMiddle);
       
        if (boxlo[dim] < cut) findBoxOverlapProcs(boxlo, boxhi, overlap, start, procMiddle - 1);
        if (boxhi[dim] > cut) findBoxOverlapProcs(boxlo, boxhi, overlap, procMiddle, end);
    }
}

/* Find the intersection of Boxes defined by ghostBox and the domain box of the "targetBox"th RCB domain
 * store the results in slablo and slabhi vector */
void CommRCB::boxIntersection(HS_float* ghostLower, 
                              HS_float* ghostUpper, 
                              int targetBox, 
                              HS_float* slablo, 
                              HS_float* slabhi)
{
    /* By contraction, the split[0-2] represents lower bound while split[3-5] denotes upper bound */
    slablo[0] = std::max(tree[targetBox].split[0], ghostLower[0]);
    slablo[1] = std::max(tree[targetBox].split[1], ghostLower[1]);
    slablo[2] = std::max(tree[targetBox].split[2], ghostLower[2]);
    slabhi[0] = std::min(tree[targetBox].split[3], ghostUpper[0]);
    slabhi[1] = std::min(tree[targetBox].split[4], ghostUpper[1]);
    slabhi[2] = std::min(tree[targetBox].split[5], ghostUpper[2]);
    printf("[slab %d]%lf %lf %lf %lf %lf %lf\n", rank, slablo[0], slabhi[0], slablo[1], slabhi[1], slablo[2], slabhi[2]);
    printf("[tree %d]%lf %lf %lf %lf %lf %lf\n", rank, tree[targetBox].split[0], tree[targetBox].split[3], tree[targetBox].split[1], tree[targetBox].split[4], tree[targetBox].split[2], tree[targetBox].split[5]);
    printf("[ghos %d]%lf %lf %lf %lf %lf %lf\n", rank, ghostLower[0], ghostUpper[0], ghostLower[1], ghostUpper[1], ghostLower[2], ghostUpper[2]);
}

/* Tell if the two box touches with me */
bool CommRCB::boxTouch(int p1, int dim, int dir)
{
    if (dir == 0)
    {
        if (closeEnough(tree[rank].split[dim], tree[p1].split[3+dim])) return true;
        else if (closeEnough(tree[rank].split[dim], 0.0) && closeEnough(tree[p1].split[3+dim], boxRange[dim])) return true;
        else return false;
    }
    else
    {
        if (closeEnough(tree[rank].split[3+dim], tree[p1].split[dim])) return true;
        else if (closeEnough(tree[p1].split[dim], 0.0) && closeEnough(tree[rank].split[3+dim], boxRange[dim])) return true;
        else return false;
    }
}
