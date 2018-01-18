#include "Atom.h"
#include "InputManager.h"
#include "Type.h"
#include "Random.h"
#include "mpi.h"
#include "RCBTree.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace Hamilton_Space;

Atom::Atom(std::shared_ptr<InputManager> input)
{
    mass = input->mass;
    natoms = input->numAtoms;
    nlocal = 0;
    nghost = 0;
    nall = nlocal + nghost;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    box.range.resize(3);
    box.length.resize(3);
    box.length[0] = box.length[1] = box.length[2] = input->boxLength;
    for (auto& item: box.range) item.resize(2);

    x.resize(MAX_ARRAY);
    v.resize(MAX_ARRAY);
    f.resize(MAX_ARRAY);
    for (int i=0; i < MAX_ARRAY; i++)
    {
        x[i].resize(3);
        v[i].resize(3);
        f[i].resize(3);
    }
}

Atom::~Atom()
{
    printf("[Atom] Atom object get destructed\n");
}

void Atom::genInitialConfig(std::shared_ptr<InputManager> input)
{
    FILE* fp = fopen("initial.dat", "r");
    HS_float tx, ty, tz;
    HS_float factor = 1.0 / std::sqrt(mass * input->beta);
    for (int i=0; i<natoms; i++)
    {
        fscanf(fp, "%lf %lf %lf", &tx, &ty, &tz);
        bool bx = (tx > box.range[0][1]) || (tx < box.range[0][0]);
        bool by = (ty > box.range[1][1]) || (ty < box.range[1][0]);
        bool bz = (tz > box.range[2][1]) || (tz < box.range[2][0]);
        if (!(bx || by || bz))
        {
            x[nlocal][0] = tx;
            x[nlocal][1] = ty;
            x[nlocal][2] = tz;
            v[nlocal][0] = factor * randomGauss();
            v[nlocal][1] = factor * randomGauss();
            v[nlocal][2] = factor * randomGauss();
            nlocal++;
        }
        //printf("here %d local %d\n", i, nlocal);
    }
    fclose(fp);
    nall = nlocal + nghost;
    /* Default: Generate Simple Cublc Lattice */

    //int nperdim = std::cbrt(nlocal);

    //std::vector<HS_float> gridSize(3);
    //for (int idim = 0; idim < 3; idim++) {
    //    gridSize[idim] = (box.range[idim][1] - box.range[idim][0]) / nperdim;
    //}

    //int nx, ny, nz;
    //HS_float factor = 1.0 / std::sqrt(mass * input->beta);
    //for (int i=0; i<nlocal; i++)
    //{
    //    nz = i / nperdim / nperdim;
    //    ny = (i - nz * nperdim * nperdim) / nperdim;
    //    nx = i % nperdim;

    //    x[i][0] = (nx + 0.5) * gridSize[0] + box.range[0][0];
    //    x[i][1] = (ny + 0.5) * gridSize[1] + box.range[1][0];
    //    x[i][2] = (nz + 0.5) * gridSize[2] + box.range[2][0];
    //   
    //    /* Initial Velocity from Maxwell Boltzmann Distribution */ 
    //    v[i][0] = factor * randomGauss();
    //    v[i][1] = factor * randomGauss();
    //    v[i][2] = factor * randomGauss();
    //}
}

void Atom::clearGhost()
{
    nghost = 0;
    nall = nlocal;
}

/* Wrap up the atoms to be sent in bufferSend */
void Atom::packSendAtoms(int first,
			 int last,
                         int dim,
			 HS_float lo,
			 HS_float hi,
                         int pbcFlag,
			 int* count,
			 HS_float* bufferSend,
			 int* sendList)
{
    *count = 0;
    bool pbcAny = pbcFlag & PBC_ANY_FLAG;
    int pbcx, pbcy, pbcz;
    pbcx = pbcy = pbcz = 0;
    
    int nsend = 0; 
    //printf("%d %d pbcAny = %d\n", first, last, pbcAny);

    if (pbcAny)
    {
        if (pbcFlag & PBC_POS_X) pbcx = 1;
        else if (pbcFlag & PBC_NEG_X) pbcx = -1;
        if (pbcFlag & PBC_POS_Y) pbcy = 1;
        else if (pbcFlag & PBC_NEG_Y) pbcy = -1;
        if (pbcFlag & PBC_POS_Z) pbcz = 1;
        else if (pbcFlag & PBC_NEG_Z) pbcz = -1;

        printf("x y z = %d %d %d\n", pbcx, pbcy, pbcz);

        for (int i = first; i <= last; i++)
        {
            if (x[i][dim] < hi && x[i][dim] > lo)
            {
                sendList[nsend++] = i;
                bufferSend[(*count)++] = x[i][0] + pbcx * box.length[0];
                bufferSend[(*count)++] = x[i][1] + pbcy * box.length[1];
                bufferSend[(*count)++] = x[i][2] + pbcz * box.length[2];
                //printf("here %lf %lf %lf count = %d\n", x[i][0], x[i][1], x[i][2], *count);
            }
        }
    }
    else
    {
        for (int i = first; i <= last; i++)
        {
            if (x[i][dim] < hi && x[i][dim] > lo)
            {
                sendList[nsend++] = i;
                bufferSend[(*count)++] = x[i][0];
                bufferSend[(*count)++] = x[i][1];
                bufferSend[(*count)++] = x[i][2];
            }
        }
    }
}

/* Wrap up the atoms to be sent in bufferSend (for RCB load balancer)*/
void Atom::packSendAtomsRCB(int first,
			    int last,
                            int dim,
			    HS_float* lo,
			    HS_float* hi,
                            int pbcFlag,
			    int* count,
			    HS_float* bufferSend,
			    int* sendList)
{
    *count = 0;
    bool pbcAny = pbcFlag & PBC_ANY_FLAG;
    int pbcx, pbcy, pbcz;
    pbcx = pbcy = pbcz = 0;
    
    int nsend = 0; 
    //printf("%d %d pbcAny = %d\n", first, last, pbcAny);
    if (pbcAny)
    {
        if (pbcFlag & PBC_POS_X) pbcx = 1;
        else if (pbcFlag & PBC_NEG_X) pbcx = -1;
        if (pbcFlag & PBC_POS_Y) pbcy = 1;
        else if (pbcFlag & PBC_NEG_Y) pbcy = -1;
        if (pbcFlag & PBC_POS_Z) pbcz = 1;
        else if (pbcFlag & PBC_NEG_Z) pbcz = -1;

        printf("x y z = %d %d %d\n", pbcx, pbcy, pbcz);

        for (int i = first; i <= last; i++)
        {
            bool isSend = true;
            for (int dim = 0; dim<3; dim++)
            {
                if (x[i][dim] > hi[dim] || x[i][dim] < lo[dim]) isSend = false;
            }
            if (isSend)
            {
                //printf("sending\n");
                sendList[nsend++] = i;
                bufferSend[(*count)++] = x[i][0] + pbcx * box.length[0];
                bufferSend[(*count)++] = x[i][1] + pbcy * box.length[1];
                bufferSend[(*count)++] = x[i][2] + pbcz * box.length[2];
                //printf("here %lf %lf %lf count = %d\n", x[i][0], x[i][1], x[i][2], *count);
            }
        }
    }
    else
    {
        for (int i = first; i <= last; i++)
        {
            bool isSend = true;
            for (int dim = 0; dim<3; dim++)
            {
                if (x[i][dim] > hi[dim] || x[i][dim] < lo[dim]) isSend = false;
            }
            if (isSend)
            {
                //printf("sending\n");
                sendList[nsend++] = i;
                bufferSend[(*count)++] = x[i][0];
                bufferSend[(*count)++] = x[i][1];
                bufferSend[(*count)++] = x[i][2];
            }
        }
    }
    //printf("count = %d\n", *count);
}
/* Unwrap the receive buffer information in local storage */
void Atom::unpackRecvAtoms(int count,
                           HS_float* buffer)
{
    for (int i=0; i<count/3; i++)
    {
        x[nall + i][0] = buffer[3*i];
        x[nall + i][1] = buffer[3*i + 1];
        x[nall + i][2] = buffer[3*i + 2];
    }
    nghost += count/3;
    nall += count/3;
}

/* Wrap up the atoms to be exchanged in Buffer
 * Carefully deal with the holes left in current x,v,f array */
void Atom::packExchange(HS_float* buffer,
		        int* nsend,
			int dimDirectionIndex)
{
    int i = 0;
    *nsend = 0;
    int dim = dimDirectionIndex / 2;
    int direction = dimDirectionIndex % 2;

    while (i < nlocal)
    {
        bool exchange_down = (direction == 0) && (x[i][dim] < box.range[dim][0]);
        bool exchange_up = (direction == 1) && (x[i][dim] > box.range[dim][1]);
        if (exchange_down || exchange_up)
        {
            
            buffer[(*nsend)++] = x[i][0];
            buffer[(*nsend)++] = x[i][1];
            buffer[(*nsend)++] = x[i][2];
            buffer[(*nsend)++] = v[i][0];
            buffer[(*nsend)++] = v[i][1];
            buffer[(*nsend)++] = v[i][2];
          
            /* Exchange the information of x[i] with x[nlocal-1] then nlocal-- 
             * Discarding the i atom */

            swap(i, nlocal-1);
            nlocal--;
        }
        else
        {
            i++;
        }
    }
}

/* swap the information of atom i and j */
void Atom::swap(int i, int j)
{
    //HS_float temp;  
    // Super weired!!! Compiler Need to use optimize -O2 to make this function behave correctly
    for (int dim = 0; dim < 3; dim++)
    {
        HS_float tempx = x[i][dim];
        x[i][dim] = x[j][dim];
        x[j][dim] = tempx;

        HS_float tempv = v[i][dim];
        v[i][dim] = v[j][dim];
        v[j][dim] = tempv;

        HS_float tempf = f[i][dim];
        f[i][dim] = f[j][dim];
        f[j][dim] = tempf;
    }
}

/* Unpack the incoming atoms */
void Atom::unpackExchange(int count,
			  HS_float* buffer)
{
    for (int i=0; i<count/6; i++)
    {
        x[nlocal + i][0] = buffer[6*i];
        x[nlocal + i][1] = buffer[6*i + 1];
        x[nlocal + i][2] = buffer[6*i + 2];

        v[nlocal + i][0] = buffer[6*i + 3];
        v[nlocal + i][1] = buffer[6*i + 4];
        v[nlocal + i][2] = buffer[6*i + 5];
    }
    nlocal += count/6;
    nall += count/6;
}

void Atom::packCommunicateGhosts(int count, 
				 int dim, 
				 int pbcFlag, 
				 HS_float* bufferSend, 
				 int* sendList)
{
    bool pbcAny = pbcFlag & PBC_ANY_FLAG;
    int pbcx, pbcy, pbcz;
    pbcx = pbcy = pbcz = 0;

    int nsend = 0;
    
    if (pbcAny)
    {
        if (pbcFlag & PBC_POS_X) pbcx = 1;
        else if (pbcFlag & PBC_NEG_X) pbcx = -1;
        if (pbcFlag & PBC_POS_Y) pbcy = 1;
        else if (pbcFlag & PBC_NEG_Y) pbcy = -1;
        if (pbcFlag & PBC_POS_Z) pbcz = 1;
        else if (pbcFlag & PBC_NEG_Z) pbcz = -1;

        for (int j = 0; j <= count; j++)
        {
            int i = sendList[j];
            bufferSend[nsend++] = x[i][0] + pbcx * box.length[0];
            bufferSend[nsend++] = x[i][1] + pbcy * box.length[1];
            bufferSend[nsend++] = x[i][2] + pbcz * box.length[2];
        }
    }
    else
    {
        for (int j = 0; j <= count; j++)
        {
            int i = sendList[j];
            bufferSend[nsend++] = x[i][0];
            bufferSend[nsend++] = x[i][1];
            bufferSend[nsend++] = x[i][2];
        }
    }
}

void Atom::unpackCommunicateGhosts(int count, 
				   int startIndex, 
				   HS_float* buffer)
{
    for (int i=0; i<count; i++)
    {
        x[startIndex + i][0] = buffer[3*i];
        x[startIndex + i][1] = buffer[3*i + 1];
        x[startIndex + i][2] = buffer[3*i + 2];
    }
}

void Atom::packReverseCommunication(int count, 
				    int startIndex, 
                                    HS_float* buffer)
{
    for (int i=0; i<count; i++)
    {
        buffer[3*i]   = f[startIndex + i][0];
        buffer[3*i+1] = f[startIndex + i][1];
        buffer[3*i+2] = f[startIndex + i][2];
    }
}

void Atom::unpackReverseCommunication(int count, 
				      int* sendList, 
				      HS_float* buffer)
{
    for (int i=0; i<count; i++)
    {
        int j = sendList[i];
        f[j][0] += buffer[3*i];
        f[j][1] += buffer[3*i+1];
        f[j][2] += buffer[3*i+2];
    }
}

void Atom::packExchangeRCB(std::vector<HS_float*> &bufferSend,
                           std::vector<int> &sendNum,
                           int dim,
                           RCBTreeNode* tree)
{
    int i = 0;
    while (i < nlocal)
    {
         if (x[i][dim] < 0.0 || x[i][dim] > box.length[dim])
         {
             int proc = findParticleInRCBTree(x[i], tree, 0, nprocs - 1);

             bufferSend[proc][sendNum[proc]++] = x[i][0];
             bufferSend[proc][sendNum[proc]++] = x[i][1];
             bufferSend[proc][sendNum[proc]++] = x[i][2];
             bufferSend[proc][sendNum[proc]++] = v[i][0];
             bufferSend[proc][sendNum[proc]++] = v[i][1];
             bufferSend[proc][sendNum[proc]++] = v[i][2];

             swap(i, nlocal - 1);
             nlocal--;
         }
         else 
         {
             i++;
         }
    }
}

void Atom::unpackExchangeRCB(std::vector<HS_float*> &bufferRecv,
                             std::vector<int> &recvNum,
                             int numOverlaps)
{
    int incr = 0;
    for (int j=0; j<numOverlaps; j++)
    {
        for (int i=0; i<recvNum[j]/6; i++)
        {
             x[nlocal + incr][0] = bufferRecv[j][6*i];
             x[nlocal + incr][1] = bufferRecv[j][6*i + 1];
             x[nlocal + incr][2] = bufferRecv[j][6*i + 2];

             v[nlocal + incr][0] = bufferRecv[j][6*i + 3];
             v[nlocal + incr][1] = bufferRecv[j][6*i + 4];
             v[nlocal + incr][2] = bufferRecv[j][6*i + 5];
             incr++;
        }
    }
    nlocal += incr;
    nall += incr;
}

// To test the parallel Implementation
//
void Atom::enforcePBC()
{
    for (int i=0; i<nlocal; i++)
    {
        if (x[i][0] < 0.0) x[i][0] += box.length[0];
        if (x[i][0] > box.length[0]) x[i][0] -= box.length[0];
        if (x[i][1] < 0.0) x[i][1] += box.length[1];
        if (x[i][1] > box.length[1]) x[i][1] -= box.length[1];
        if (x[i][2] < 0.0) x[i][2] += box.length[2];
        if (x[i][2] > box.length[2]) x[i][2] -= box.length[2];
    }
}

void Atom::printFrame()
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    char filename[100];
    sprintf(filename, "traj_p%d.dat", me);
    FILE* fp = fopen(filename, "w");
    fprintf(fp,"nlocal = %d nghost = %d\n", nlocal, nghost);
    for (int i=0; i<nall; i++)
    {
        fprintf(fp, "%lf %lf %lf %d\n", x[i][0], x[i][1], x[i][2], i<nlocal);
    }
    fclose(fp);
}

void Atom::propagate()
{
    for (int i=0; i<nlocal; i++)
    {
        x[i][0] = x[i][0] + v[i][0] * 0.1;
        x[i][1] = x[i][1] + v[i][1] * 0.1;
        x[i][2] = x[i][2] + v[i][2] * 0.1;
    }
}

// Walk the binary tree to find the proc id of the particle
int Atom::findParticleInRCBTree(std::vector<HS_float> x, RCBTreeNode* tree, int start, int end)
{
    if (start == end)
        return start;
    else
    {
        int procMiddle = start + (end - start) / 2 + 1;
        HS_float cut = tree[procMiddle].cut;
        int dim = tree[procMiddle].dimension;

        if (x[dim] < cut) return findParticleInRCBTree(x, tree, start, procMiddle - 1);
        else  return findParticleInRCBTree(x, tree, procMiddle, end);
    }
}
