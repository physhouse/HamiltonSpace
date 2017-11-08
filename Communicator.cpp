#ifndef _COMMUNICATOR_H_
#define _COMMUNICATOR_H_

// This class implements the Grid-type Domain Decomposition Parallelization of General Molecular Dynamics and Monte Carlo Simulations
// 
// 2017-11-06

#include <mpi.h>

using namespace HamiltonSpace;

Communicator::Communicator()
{
}

Communicator::~Communicator()
{
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
    HS_float lx, ly, lz;
    boxLength[0] = atom->xhi - atom->xlo;
    boxLength[1] = atom->yhi - atom->ylo;
    boxLength[2] = atom->zhi - atom->zlo;
 
    HS_float bestSurface = INFINITY;
    for (int nx=1; nx<=nproc; nx++)
    {
        for (int ny=1; ny<=nproc/nx; ny++)
        {
            HS_float surf = lx * ly / nx / ny + lx * lz / nx / nz + ly * lz / ny / nz;
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
    MPI_Cart_shift(GridWorld, 1, 1, &front, &back);
    MPI_Cart_shift(GridWorld, 2, 1, &up, &down);
   
    neighbor[0][0] = left;
    neighbor[0][1] = right;
    neighbor[1][0] = front;
    neighbor[1][1] = back;
    neighbor[2][0] = down;
    neighbor[2][1] = up;


    /* Scheme of the neighbors
              UP
              |    FRONT
              |   /
              |  /
              | /
    LEFT------ME------RIGHT
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
    need[0] = static_cast<int>(cutoff * processorGrid[0] / lx + 1);
    need[1] = static_cast<int>(cutoff * processorGrid[1] / lx + 1);
    need[2] = static_cast<int>(cutoff * processorGrid[2] / lx + 1);

    /* The total number of swaps need to be performed */
    int maxSwap = 2 * (need[0] + need[1] + need[2]);
    sendToProc.resize(maxSwap);
    recvFromProc.resize(maxSwap);
    pbcFlags.resize(maxSwap);
    
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
                sendToProc[iSwap] = neighbor[idim][0];
                recvFromProc[iSwap] = neighbor[idim][1];
                sourceSlab = my3DLocation[idim] + j / 2;
                lo = sourceSlab * boxLength[idim] / processorGrid[idim];
                hi = std::min(atom.box.lo[idim] + cutoff, (sourceSlab + 1) * boxLength[idim] / processorGrid[idim]);
                
                if (my3DLocation[idim] == 0)
                {
                    pbcFlags[iSwap] |= PBC_ANY_FLAG;
                    if (idim == 0) pbcFlags[iSwap] |= PBC_POS_X;
                    if (idim == 1) pbcFlags[iSwap] |= PBC_POS_Y;
                    if (idim == 2) pbcFlags[iSwap] |= PBC_POS_Z;
                }
            }
            else // Second-half: Receive from left and send to right
            {
                sendToProc[iSwap] = neighbor[idim][1];
                recvFromProc[iSwap] = neighbor[idim][0];
                sourceSlab = my3DLocation[idim] - j / 2;
                lo = std::max(atom.box.hi[idim] - cutoff, sourceSlab * boxLength[idim] / processorGrid[idim]);
                hi = (sourceSlab + 1) * boxLength[idim] / processorGrid[idim];

                if (my3DLocation[idim] == processorGrid[idim] - 1)
                {
                    pbcFlags[iSwap] |= PBC_ANY_FLAG;
                    if (idim == 0) pbcFlags[iSwap] |= PBC_NEG_X;
                    if (idim == 1) pbcFlags[iSwap] |= PBC_NEG_Y;
                    if (idim == 2) pbcFlags[iSwap] |= PBC_NEG_Z;
                }
            }
 
            slablo[iSwap] = lo;
            slabhi[iSwap] = hi;
            iSwap++;
        }
    }

    // Total number of Swaps need to be performed
    nSwaps = iSwap;
         
}

/* Routine to generate the ghost atom lists */
void Communicator::generateGhosts(std::shared_ptr<Atom> atom)
{
    atom->nghost = 0; // Clear up the ghost atom lists
    
    for 
} 

/* Routine to migrate atoms to neighboring proc */
void Communicator::exchangeAtoms(std::shared_ptr<Atom> atom)
{
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
#endif

