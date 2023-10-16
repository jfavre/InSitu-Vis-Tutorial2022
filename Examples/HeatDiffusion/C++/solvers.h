#ifndef _SOLVERS_H_INCLUDED_
#define _SOLVERS_H_INCLUDED_

#include <set>
#include <string>
 
typedef struct
{
    int       par_rank;
    int       par_size;
    MPI_Comm  topocomm;
    int       east, west, south, north;
    int       cart_dims[2], rankx, ranky;
    int      bx, by;
    int      resolution; // overall grid size not counting boundary walls.
    int      iter;
    double   gdel, *oldTemp, *Temp;
    unsigned char *Ghost; // ghost-point array
    double    *cx, *cy;
    int      local_extents[6];
    std::string mesh;  // "uniform", "rectilinear", "structured", "unstructured"
    // for unstructured mesh only
    int      *connectivity;
    // for structured and unstructured mesh only
    double    *explicit_cx, *explicit_cy;
    bool verbose;
} simulation_data;

#define BASENAME "."
//#define BASENAME "/scratch/snx3000/jfavre"

#define INCREMENT 10       /* number of steps between convergence check  */
#define MAXSTEPS 1000     /* Maximum number of iterations               */
#define TOL 1e-06          /* Numerical Tolerance */

void SimInitialize(simulation_data *sim);
void MPI_Partition(int PartitioningDimension, simulation_data *sim);
void AllocateGridMemory(simulation_data *sim);
void FreeGridMemory(simulation_data *sim);
void set_initial_bc(simulation_data *sim);
double update_temperature(simulation_data *sim);
void simulate_one_timestep(simulation_data *sim);
void CopyTempValues_2_OldValues(simulation_data *sim);
void WriteFinalGrid(simulation_data *sim);

void exchange_ghost_lines(simulation_data *sim);
void neighbors(simulation_data *sim);
void MPIIOWriteData(const char *filename, simulation_data *sim);
extern MPI_Datatype rowtype, coltype;
#endif
