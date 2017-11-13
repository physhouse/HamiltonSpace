// This class implements the Grid-type Domain Decomposition Parallelization of General Molecular Dynamics and Monte Carlo Simulations
// 
// 2017-11-06

#include <Communicator.h>
#include <mpi.h>

using namespace HamiltonSpace;

Communicator::Communicator()
{
    HS_float cutoff = InputManager->cutNeighbor;
    bufferExchangeSend = new HS_float(EXCHANGEMAX);
    bufferExchangeRecv = new HS_float(EXCHANGEMAX);
}

Communicator::~Communicator()
{
    for (auto &item : swaps)
    {
        delete[] item.bufferSend;
        delete[] item.bufferRecv;
    }

    delete[] bufferExchangeSend;
    delete[] bufferExchangeRecv;
}

void Communicator::setup(HS_float cutoff, std::shared_ptr<Atom> atom)
{
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Step 1: Determine the box dimension of the current spatial Domain
    // Algorithm aims at finding the spatial decomposition such that the total
    // surface area of these domains reaches the minimum, under such situations,
    // the particles near the boundaries are at the lowest ratio

    std::vector<int> processorGrid(3);
    std::vector<HS_float> boxLength[3];
    boxLength[0] = atom->box.lengthx;
    boxLength[1] = atom->box.lengthy;
    boxLength[2] = atom->box.lengthz;
 
    HS_float bestSurface = INFINITY;
    for (int nx=1; nx<=nproc; nx++)
    {
        for (int ny=1; ny<=nproc/nx; ny++)
        {
            HS_float surf = boxLength[0] * boxLength[1] / nx / ny + boxLength[0] * boxLength[2] / nx / nz 
                          + boxLength[1] * boxLength[2] / ny / nz;

            if (surf < bestSurface)
            {
                processorGrid[0] = nx;
                processorGrid[1] = ny;
                processorGrid[2] = nz;

		bestSurface = surf;
            }
        }
    }
    
    assert(processorGrid[0]*processorGrid[1]*processorGird[2] == nprocs);
    printf("[Communicator] Grid Partitioning Initialized (nx, ny, nz) = (%d, %d, %d)\n", 
	   processorGrid[0], processorGrid[1], processorGrid[2]);

    // Step 2: Generate the Grid Configuration for the Processors
    // Determining All the adjacent neighbors processor IDs and store them in the neighbors map
    // This neighbors map is used for determinig the particle migration step of communication
    
    // int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[],
    //                     const int periods[], int reorder, MPI_Comm *comm_cart)
    MPI_Comm GridWorld;
    bool periods[3] = {true, true, true};
    bool reorder = false;
    MPI_Cart_create(MPI_COMM_WORLD, 3, processorGrid.data(), periods, reorder, &GridWorld);

    int my3DLocation[3];
    MPI_Cart_get(GridWorld, 3, processorGrid.data(), periods, my3DLocation);
    
    int up, down, left, right, front, back;
    MPI_Cart_shift(GridWorld, 0, 1, &left, &right);
    MPI_Cart_shift(GridWorld, 1, 1, &back, &front);
    MPI_Cart_shift(GridWorld, 2, 1, &down, &up);
   
    neighbor[0][0] = left;
    neighbor[0][1] = right;
    neighbor[1][0] = back;
    neighbor[1][1] = front;
    neighbor[2][0] = down;
    neighbor[2][1] = up;


    /* Scheme of the neighbors
              UP(z axis)
              |    FRONT (y axis)
              |   /
              |  /
              | /
    LEFT------ME------RIGHT (x axis)
             /|
            / |
           /  |
         BACK |
             DOWN
    *
    */


    // Determining the geometric range of current processor
    for (int idim=0; idim<3; idim++)
    {
        atom->box.xlo = my3DLocation[idim] * boxLength[idim] / processorGrid[idim];
        atom->box.xhi = (my3DLocation[idim]+1) * boxLength[idim] / processorGrid[idim];
    }

    // Step 3: Generate the send and recv processor IDS
    // Determining all the processor IDs to get atoms from when building ghost atoms
    
    /* Determine how many cells I need at each direction*/
    std::vector<int> need(3,0);
    need[0] = static_cast<int>(cutoff * processorGrid[0] / boxLength[0] + 1);
    need[1] = static_cast<int>(cutoff * processorGrid[1] / boxLength[1] + 1);
    need[2] = static_cast<int>(cutoff * processorGrid[2] / boxLength[2] + 1);

    /* The total number of swaps need to be performed */
    int maxSwap = 2 * (need[0] + need[1] + need[2]);
    swap.resize(maxSwap);
    for (auto &item : swap)
    {
        item.bufferSend = new HS_float(BUFFMAX);
        item.bufferRecv = new HS_float(BUFFMAX);
    }
    
    /* Determining whether a PBC wrapping is needed for each swap 
       Also determine the slab region */

    /* Note: It's important to understand the algorithm for generating ghost lists first
     *
     *  PBC Concerns: When all the boundary processors owns the necessary peoriodic images in the ghost particle limits
     *                The particle coordinates are already wrapped and no need to worry about minimum image at all, this is advantages of domain decomp
     *                Only do the PBC correction when SENDing atoms across the boundary
     *                Once it is treated while sending, no need to do it again while receiving
     *
     *  For each swap, just send one block of ghost atoms to the target neighbor processor
     *  The first swap send all the local particles of this proc to the target
     *  The second swap send the ghost particles that was received from the first swap to the target
     *  The third swap send the ghost particles that was received from the second swap to the target... etc
     *
     *  slablo[swap] and slabhi[swap] marks the coordinate range of particles that should be sent on the current swap
     *  they are the 
     *
     *  For the first half, send to the left and receive from right, like the situation
     *
     *  |   |  Me  |     |     |     |
     *
     *      | Cutoff  | Only need particles within these range
     *        
     *  therefore, only the particles with x < currentbox.lo + cutoff need to be sended to the left box
     *  Similar patterns for the second half while sending to the right and receiving from left
     *
     *  PBC Concerns: Since only sending to the left, only need to think about PBC while "me" is on the left end of this dimension
     *
     */

    int iSwap = 0;
    int sourceSlab;	
    int lo, hi;
    for (int idim = 0; idim < 3; idim++)
    {
        for (int j = 0; j < 2 * need[idim]; j++)
        {
            if (j % 2 == 0)  // First-half: Receive from right and send to left
            {
                swap[iSwap].sendToProc = neighbor[idim][0];
                swap[iSwap].recvFromProc = neighbor[idim][1];
                sourceSlab = my3DLocation[idim] + j / 2;

                lo = sourceSlab * boxLength[idim] / processorGrid[idim];
                hi = std::min(atom.box.lo[idim] + cutoff, (sourceSlab + 1) * boxLength[idim] / processorGrid[idim]);
                
                if (my3DLocation[idim] == 0)
                {
                    swap[iSwap].pbcFlags |= PBC_ANY_FLAG;
                    if (idim == 0) swap[iSwap].pbcFlags |= PBC_POS_X;
                    if (idim == 1) swap[iSwap].pbcFlags |= PBC_POS_Y;
                    if (idim == 2) swap[iSwap].pbcFlags |= PBC_POS_Z;
                }
            }
            else // Second-half: Receive from left and send to right
            {
                swap[iSwap].sendToProc = neighbor[idim][1];
                swap[iSwap].recvFromProc = neighbor[idim][1];
                sourceSlab = my3DLocation[idim] - j / 2;
                lo = std::max(atom.box.hi[idim] - cutoff, sourceSlab * boxLength[idim] / processorGrid[idim]);
                hi = (sourceSlab + 1) * boxLength[idim] / processorGrid[idim];

                if (my3DLocation[idim] == processorGrid[idim] - 1)
                {
                    swap[iSwap].pbcFlags |= PBC_ANY_FLAG;
                    if (idim == 0) swap[iSwap].pbcFlags |= PBC_NEG_X;
                    if (idim == 1) swap[iSwap].pbcFlags |= PBC_NEG_Y;
                    if (idim == 2) swap[iSwap].pbcFlags |= PBC_NEG_Z;
                }
            }
 
            swap[iSwap].dim = idim;
            swap[iSwap].slablo = lo;
            swap[iSwap].slabhi = hi;
            iSwap++;
        }
    }

    // Total number of Swaps need to be performed
    nSwaps = iSwap;
         
}

/* Routine to generate the ghost atom lists */
void Communicator::generateGhosts(std::shared_ptr<Atom> atom)
{
    atom->clearGhost(); // Clear up the ghost atom lists
    
    int first = 0;
    int last = atom->nlocal;
    int count;
    MPI_Status status;
    MPI_Request request;

    for (const auto& item : swap)
    {
        atom->packSendAtoms(first, last, item.dim, item.slablo, item.slabhi, item.pbcFlag, &item.sendNum, item.bufferSend);
       
        MPI_Send(&item.sendNum, 1, MPI_INT, item.sendToProc, 0, MPI_COMM_WORLD);
        MPI_Recv(&item.recvNum, 1, MPI_INT, item.recvFromProc, 0, MPI_COMM_WORLD, &status);
       
        MPI_Irecv(item.bufferRecv, item.recvNum, MPI_DOUBLE, item.recvFromProc, 0, MPI_COMM_WORLD, &request);
        MPI_Send(item.bufferSend, sendNum[iSwap], MPI_DOUBLE, item.sendToProc, 0, MPI_COMM_WORLD);
      
        MPI_Wait(&request, &status); 

        atom->unpackRecvAtoms(item.recvNum, item.bufferRecv);

        first = last;
        last = last + item.recvNum;
    }
} 

/* Routine to migrate atoms to neighboring proc */
/* This function is only called when rebuilding the neighbor lists together with generateGhosts() */
void Communicator::exchangeAtoms(std::shared_ptr<Atom> atom)
{
    MPI_Status status;
    MPI_Request request;
    int nsend, nrecv;
    for (int idim = 0; idim < 3; idim++)
    {
        if (processorGrid[idim] > 1)
        {
            /* Packing strategy:
             * (1) When atom is about to leave this box, exchange it with atom at nlocal, nlocal--
             * (2) Buffer the information of the leaving atom, nExchange++
             * (3) After collecting the sending buffer, MPI_Send to the corresponding neighbors
             */
            atom->packExchange(bufferExchangeSend, &nsend, 2*idim);

            MPI_Send(&nsend, 1, MPI_INT, neighbor[idim][0], 0, MPI_COMM_WORLD);
            MPI_Recv(&nrecv, 1, MPI_INT, neighbor[idim][1], 0, MPI_COMM_WORLD, &status);

            MPI_Irecv(bufferExchangeRecv, nrecv, MPI_DOUBLE, neighbor[idim][1], 0, MPI_COMM_WORLD, &request);
            MPI_Send(bufferExchangeSend, nsend, MPI_DOUBLE, neighbor[idim][0], 0, MPI_COMM_WORLD);

            /* Unpacking strategy:
             * Just check the incoming atoms and append it at the end, nlocal += nrecv
             */
            atom->unpackExchange(nrecv, bufferExchangeRecv);
          
            // Round 2: Send to right
            atom->packExchange(bufferExchangeSend, &nsend, 2*idim + 1);

            MPI_Send(&nsend, 1, MPI_INT, neighbor[idim][1], 0, MPI_COMM_WORLD);
            MPI_Recv(&nrecv, 1, MPI_INT, neighbor[idim][0], 0, MPI_COMM_WORLD, &status);

            MPI_Irecv(bufferExchangeRecv, nrecv, MPI_DOUBLE, neighbor[idim][0], 0, MPI_COMM_WORLD, &request);
            MPI_Send(bufferExchangeSend, nsend, MPI_DOUBLE, neighbor[idim][1], 0, MPI_COMM_WORLD);

            atom->unpackExchange(nrecv, bufferExchangeRecv);
        }
    }
}

void Communicator::communicate()
{
}


void Communicator::communicateGhosts(std::shared_ptr<Atom> atom)
{
}

void Communicator::reverseCommunicateGhosts(std::shared_ptr<Atom> atom)
{
}

}

