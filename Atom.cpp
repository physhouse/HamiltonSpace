#include "Atom.h"
#include "InputManager.h"
#include "Type.h"
#include "Random.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace Hamilton_Space;

Atom::Atom(std::shared_ptr<InputManager> input)
{
    mass = input->mass;
    nlocal = input->atomPerProc;
    nghost = 0;
    nall = nlocal + nghost;

    box.lengthx = box.lengthy = box.lengthz = input->boxLength;
    box.range.resize(3);
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
    /* Default: Generate Simple Cublc Lattice */
    int nperdim = std::cbrt(nlocal);
    HS_float gridSize = box.lengthx / nperdim;
    int nx, ny, nz;
    for (int i=0; i<nlocal; i++)
    {
        nz = i / nperdim / nperdim;
        ny = (i - nz * nperdim * nperdim) / nperdim;
        nx = i % nperdim;

        x[i][0] = (nx + 0.5) * gridSize + box.range[0][0];
        x[i][1] = (ny + 0.5) * gridSize + box.range[1][0];
        x[i][2] = (nz + 0.5) * gridSize + box.range[2][0];
       
        /* Initial Velocity from Maxwell Boltzmann Distribution */ 
        HS_float factor = 1.0 / std::sqrt(mass * input->beta);
        v[i][0] = factor * randomGauss();
        v[i][1] = factor * randomGauss();
        v[i][2] = factor * randomGauss();
    }
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
			 HS_float* bufferSend)
{
    *count = 0;
    bool pbcAny = pbcFlag & PBC_ANY_FLAG;
    int pbcx, pbcy, pbcz;
    pbcx = pbcy = pbcz = 0;
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
                bufferSend[(*count)++] = x[i][0] + pbcx * box.lengthx;
                bufferSend[(*count)++] = x[i][1] + pbcy * box.lengthy;
                bufferSend[(*count)++] = x[i][2] + pbcz * box.lengthz;
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
                bufferSend[(*count)++] = x[i][0];
                bufferSend[(*count)++] = x[i][1];
                bufferSend[(*count)++] = x[i][2];
            }
        }
    }
}

/* Unwrap the receive buffer information in local storage */
void Atom::unpackRecvAtoms(int count,
                           HS_float* buffer)
{
    for (int i=0; i<count/3; i++)
    {
        x[nghost + i][0] = buffer[3*i];
        x[nghost + i][1] = buffer[3*i + 1];
        x[nghost + i][2] = buffer[3*i + 2];
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
        bool exchange_down = (direction == 0) && (x[i][dim] < box.range[dim][direction]);
        bool exchange_up = (direction == 1) && (x[i][dim] > box.range[dim][direction]);
        if (exchange_down || exchange_up)
        {
            
            buffer[*nsend++] = x[i][0];
            buffer[*nsend++] = x[i][1];
            buffer[*nsend++] = x[i][2];
            buffer[*nsend++] = v[i][0];
            buffer[*nsend++] = v[i][1];
            buffer[*nsend++] = v[i][2];
          
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
    int temp;
    for (int dim = 0; dim < 3; dim++)
    {
        temp = x[i][dim];
        x[i][dim] = x[j][dim];
        x[j][dim] = x[i][dim];
        temp = v[i][dim];
        v[i][dim] = v[j][dim];
        v[j][dim] = v[i][dim];
        temp = f[i][dim];
        f[i][dim] = f[j][dim];
        f[j][dim] = f[i][dim];
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
