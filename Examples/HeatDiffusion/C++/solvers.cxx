#include <malloc.h>
#include <string.h>
#include <iostream>
#include <math.h>

#include <mpi.h>
MPI_Datatype rowtype, coltype;

#include "solvers.h"

void SimInitialize(simulation_data *sim)
{
  sim->par_rank = 0;
  sim->par_size = 1;

  sim->bx = sim->by = sim->resolution;
  sim->local_extents[0] = 0;
  sim->local_extents[1] = sim->bx - 1;
  sim->local_extents[2] = 0;
  sim->local_extents[3] = sim->by - 1;
  sim->local_extents[4] = 0;
  sim->local_extents[5] = 0; // 2D problem

  sim->gdel = 1.0;
  sim->iter = 0;
  sim->cx = nullptr;
  sim->cy = nullptr;
  sim->connectivity = nullptr;
  sim->explicit_cx = nullptr;
  sim->explicit_cy = nullptr;
}

void MPI_Partition(int PartitioningDimension, simulation_data *sim)
{
  int coords[2];
  int periods[2]={0,0};
  if(PartitioningDimension == 1)
    sim->cart_dims[1] = 1;

  MPI_Dims_create(sim->par_size, PartitioningDimension, sim->cart_dims);
  fprintf(stdout,"%d: cart_dims[]= %d, %d\n", sim->par_rank, sim->cart_dims[0], sim->cart_dims[1]);

  if(MPI_Cart_create(MPI_COMM_WORLD, 2, sim->cart_dims, periods, 0, &sim->topocomm) != MPI_SUCCESS)
    sim->topocomm = MPI_COMM_WORLD;

  MPI_Comm_rank(sim->topocomm, &sim->par_rank);
  MPI_Comm_size(sim->topocomm, &sim->par_size);
  MPI_Cart_coords(sim->topocomm, sim->par_rank, 2, coords);

  //fprintf(stdout,"%d: Rank_xy[]= %d, %d\n", sim.par_rank, coords[0], coords[1]);
  //fprintf(stdout,"      %2d\n%2d<->[%2d]<->%2d\n      %2d\n", sim.north, sim.west, sim.par_rank, sim.east, sim.south);
  
  sim->rankx = coords[0];
  sim->ranky = coords[1];

// We make no attempt to check that the number of grid points divides evenly
// with the number of MPI tasks.
// rank 0 will display the bottom (southern) boundary wall
// rank (size-1) will display the top (northern) boundary wall
// if run with m=20 and 4 MPI tasks, we will have 10 grid lines per rank
// and VisIt will display a 22x22 grid

  MPI_Bcast(&(sim->resolution), 1, MPI_INT, 0, sim->topocomm);
  sim->bx = sim->resolution / sim->cart_dims[0]; // block size in x
  sim->by = sim->resolution / sim->cart_dims[1]; // block size in y

  sim->local_extents[0] = sim->rankx * sim->bx;
  sim->local_extents[1] = sim->local_extents[0] + sim->bx + 1;
  sim->local_extents[2] = sim->ranky * sim->by;
  sim->local_extents[3] = sim->local_extents[2] + sim->by + 1;
  std::cout << "Extents[" << sim->rankx << ","<<sim->ranky << "] = " << sim->local_extents[0] << "," << sim->local_extents[1] << ", " << sim->local_extents[2] <<","<< sim->local_extents[3] << ", "<<sim->local_extents[4] << "," << sim->local_extents[5] <<std::endl;
  MPI_Type_contiguous(sim->bx+2, MPI_DOUBLE, &rowtype); 
  MPI_Type_commit(&rowtype);

  MPI_Type_vector(sim->by, 1, sim->bx+2, MPI_DOUBLE, &coltype); // count, blocklength, stride,
  MPI_Type_commit(&coltype);
}

void AllocateGridMemory(simulation_data *sim)
{
  int i;
  sim->oldTemp = (double *)calloc(       (sim->bx + 2) * (sim->by + 2), sizeof(double));
  sim->Temp    = (double *)calloc(       (sim->bx + 2) * (sim->by + 2), sizeof(double));
  sim->Ghost   = (unsigned char *)calloc((sim->bx + 1) * (sim->by + 1), sizeof(unsigned char));
  sim->cx      = (double *)malloc(sizeof(double) * (sim->bx + 2));
  sim->cy      = (double *)malloc(sizeof(double) * (sim->by + 2));

  float hsize = 1.0/(sim->resolution+1.0);

  if((sim->mesh == "unstructured") || (sim->mesh == "structured"))
    {
    int i=0;
    sim->explicit_cx = (double *)malloc((sim->bx + 2) * (sim->by + 2) * sizeof(double));
    sim->explicit_cy = (double *)malloc((sim->bx + 2) * (sim->by + 2) * sizeof(double));
    for(int iy=0; iy <= sim->by+1; iy++)
      for(int ix=0; ix <= sim->bx+1; ix++)
        {
        sim->explicit_cx[i] = (ix + sim->rankx * sim->bx) * hsize;
        sim->explicit_cy[i] = (iy + sim->ranky * sim->by) * hsize;
        i++;
        }
/*
    std::cout << "###############################" << std::endl;
    std::cout << "Explicit x, y and z coords of dim " << (sim->bx + 2) * (sim->by + 2) << std::endl;
    for(int d=0; d < (sim->bx + 2) * (sim->by + 2); d++)
    {
      std::cout << sim->explicit_cx[d] << " ";
    }
    std::cout << std::endl;
    for(int d=0; d < (sim->bx + 2) * (sim->by + 2); d++)
    {
      std::cout << sim->explicit_cy[d] << " ";
    }
    std::cout << std::endl;
    std::cout << "###############################" << std::endl;
*/
    }
    
    if(sim->mesh == "unstructured")
    {
    sim->connectivity = (int *)malloc(4 * (sim->bx + 1) * (sim->by + 1) * sizeof(int));
    auto i=0;
    for(int iy=0; iy <= sim->by; iy++)
       for(int ix=0; ix <= sim->bx; ix++)
          {
          sim->connectivity[i+0] = ix + iy*(sim->bx + 2);
          sim->connectivity[i+1] = ix + (iy+1)*(sim->bx + 2);
          sim->connectivity[i+2] = ix + (iy+1)*(sim->bx + 2)+ 1;
          sim->connectivity[i+3] = ix + iy*(sim->bx + 2) + 1;
          i+=4; // quads have 4 indices
          }
    }

//std::cout << "XExtents[" << sim->rankx << ","<<sim->ranky << "] = " << sim->local_extents[0] << "," << sim->local_extents[1] << " vs. " <<(sim->rankx * sim->bx) << "," << ((sim->bx + 1) + sim->rankx * sim->bx)<< std::endl;
  for(i = 0; i < (sim->bx + 2); i++)
    sim->cx[i] = (i + sim->rankx * sim->bx) * hsize;
//std::cout << "YExtents[" << sim->rankx << ","<<sim->ranky << "] = " << sim->local_extents[2] << "," << sim->local_extents[3] << " vs. " <<(sim->ranky * sim->by) << "," << ((sim->by + 1) + sim->ranky * sim->by)<< std::endl;
  for(i = 0; i < (sim->by + 2); i++)
    sim->cy[i] = (i + sim->ranky * sim->by) * hsize;
}

void FreeGridMemory(simulation_data *sim)
{
  free(sim->oldTemp);
  free(sim->Temp);
  free(sim->cx);
  free(sim->cy);
  free(sim->explicit_cx);
  free(sim->explicit_cy);
  free(sim->Ghost);
}

#define DUPLICATECELL 1
#define DUPLICATEPOINT 1
#define HIDDENCELL 32
#define HIDDENPOINT 2
void set_initial_bc(simulation_data *sim)
{
/*********** Boundary Conditions ****************************************
 *  PDE: Laplacian u = 0;      0<=x<=1;  0<=y<=1                        *
 *  B.C.: u(x,0)=sin(pi*x); u(x,1)=sin(pi*x)*exp(-pi); u(0,y)=u(1,y)=0  *
 *  SOLUTION: u(x,y)=sin(pi*x)*exp(-pi*y)                               *
 ************************************************************************/
  int i;
  double x;
  // set to 0 everywhere
  memset(sim->Ghost, 0, sizeof(unsigned char)*(sim->bx+1)*(sim->by+1));
  // set to "1" on bottom line
  memset(sim->Ghost, DUPLICATECELL, sizeof(unsigned char)*(sim->bx+1));
  // set to "1" on top line
  memset(&sim->Ghost[(sim->by+0)*(sim->bx+1)], DUPLICATECELL, sizeof(unsigned char)*(sim->bx+1));
  
  for (i = 1; i < sim->by+0; i++)
    {
    sim->Ghost[i*(sim->bx+1)] = sim->Ghost[i*(sim->bx+1)+sim->by+0] = DUPLICATECELL;
    }
    
  //memset(sim->Temp, 0, sizeof(double)*(sim->bx+2)*(sim->by+2));
  sim->iter = 0; sim->gdel = 1.0;

  if (sim->ranky == 0)
    {
    for (i = 0; i < (sim->bx+2); i++)
      {
      sim->Temp[i] = sin(M_PI*(i+ sim->rankx * sim->bx)/(sim->resolution+1));              /* at y = 0; all x */
      }
    }
  if (sim->ranky == (sim->cart_dims[1]-1)) {
    for (i = 0; i < (sim->bx+2); i++)
      {
      sim->Temp[i+(sim->bx+2)*(sim->by+1)] = sin(M_PI*sim->cx[i])*exp(-M_PI);   // at y = 1; all x
      }
    }
  memset(sim->oldTemp, 0, sizeof(double)*(sim->bx+2)*(sim->by+2));
}

void CopyTempValues_2_OldValues(simulation_data *sim)
{
  //save current solution array to buffer
  memcpy(sim->oldTemp, sim->Temp, sizeof(double)*(sim->bx+2)*(sim->by+2));
}

void simulate_one_timestep(simulation_data *sim)
{
  CopyTempValues_2_OldValues(sim);
/* compute Temp solution according to the Jacobi scheme */
  double del = update_jacobi(sim);
  /* find global max error */

  MPI_Allreduce( &del, &sim->gdel, 1, MPI_DOUBLE, MPI_MAX, sim->topocomm );
  
  exchange_ghost_lines(sim); // update lowest and uppermost grid lines
  sim->iter++;
}

double update_jacobi(simulation_data *sim)
{
  double del = 0.0;

  for(auto j = 1; j < sim->by+1; j++)
    {
    for(auto i = 1; i < sim->bx+1; i++)
      {
      sim->Temp[i+(sim->bx+2)*j] = ( sim->oldTemp[i+(sim->bx+2)*(j+1)] +
                           sim->oldTemp[(i+1)+(sim->bx+2)*j] +
                           sim->oldTemp[(i-1)+(sim->bx+2)*j] +
                           sim->oldTemp[i+(sim->bx+2)*(j-1)] )*0.25;
      del += fabs(sim->Temp[i+(sim->bx+2)*j] - sim->oldTemp[i+(sim->bx+2)*j]); /* find local max error */
      }
    }

  return del;
}

void exchange_ghost_lines(simulation_data *sim)
{
  MPI_Status status;
  /*
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 int dest, int sendtag,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int source, int recvtag, MPI_Comm comm, MPI_Status * status)
                 */
// send my last computed row north and receive from south my south boundary wall
  MPI_Sendrecv(&sim->Temp[0+sim->by*(sim->bx+2)], 1, rowtype, sim->south, 0,
               &sim->Temp[0+0*(sim->bx+2)], 1, rowtype, sim->north, 0,
                sim->topocomm, &status );
// send my first computed row south and receive from north my north boundary wall
  MPI_Sendrecv(&sim->Temp[0 + 1*(sim->bx+2) ], 1, rowtype, sim->north, 1,
               &sim->Temp[0+(sim->by+1)*(sim->bx+2)], 1, rowtype, sim->south, 1,
                sim->topocomm, &status );

// send my last computed column east and receive from west my west boundary wall
  MPI_Sendrecv(&sim->Temp[sim->bx + 1*(sim->bx+2)], 1, coltype, sim->east, 2,
               &sim->Temp[0  + 1*(sim->bx+2)], 1, coltype, sim->west, 2,
                sim->topocomm, &status );
// send my first computed column west and receive from east my east boundary wall
  MPI_Sendrecv(&sim->Temp[1  + 1*(sim->bx+2)], 1, coltype, sim->west, 3,
               &sim->Temp[(sim->bx+1) + 1*(sim->bx+2)], 1, coltype, sim->east, 3,
                sim->topocomm, &status );

}

void neighbors(simulation_data *sim)
{
  MPI_Cart_shift(sim->topocomm, 0, 1,  &sim->west, &sim->east); // in the 0-th dimension
  MPI_Cart_shift(sim->topocomm, 1, 1,  &sim->north, &sim->south);// in the 1-th dimension
}

void MPIIOWriteData(const char *filename, simulation_data *sim)
{
  // global size of array on disk
  char fname[256];
  strcpy(fname, filename);
  strcpy(&fname[strlen(filename)], ".bin");
  int dimuids[2]={sim->resolution+2, sim->resolution+2};
  int ustart[2], ucount[2];
  int disp = 0, offset=0;

  MPI_File      filehandle;
  MPI_Datatype  filetype;

  MPI_File_open(sim->topocomm, fname,
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &filehandle);
  MPI_File_set_size(filehandle, disp);

// write the grid dimensions to allow a restart
/*
  if(sim->par_rank == 0)
    MPI_File_write(filehandle, dimuids, 2, MPI_INT, MPI_STATUS_IGNORE);

  disp = 2 * sizeof(int); // offset because we just wrote 2 integers
*/
  ustart[1] = sim->rankx * (sim->bx);
  ustart[0] = sim->ranky * (sim->by);
  ucount[1] = sim->bx+2;
  ucount[0] = sim->by+2;

 // Create the subarray representing the local block
  MPI_Type_create_subarray(2, dimuids, ucount, ustart,
                           MPI_ORDER_C, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  MPI_File_set_view(filehandle, disp, MPI_DOUBLE,
                    filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(filehandle, sim->Temp, ucount[0]*ucount[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&filehandle);
  MPI_Type_free(&filetype);
}

void WriteFinalGrid(simulation_data *sim)
{
  const char *fname = BASENAME"/Jacobi";
  if(sim->par_rank == 0){
  // first write a header file in BOV format, to enable reading by VisIt
  FILE * fpbov = fopen(BASENAME"/Jacobi.bov", "w");
  fprintf(fpbov,"TIME: %f\n", 0.0); // dummy value 0.0
  fprintf(fpbov,"DATA_FILE: %s.bin\n", fname);
  fprintf(fpbov,"DATA_SIZE: %d %d %d\n", sim->resolution+2, sim->resolution+2, 1); // size of grid in IJK
  fprintf(fpbov,"DATA_FORMAT: DOUBLE\n");
  fprintf(fpbov,"VARIABLE: temperature\n");
  fprintf(fpbov,"DATA_ENDIAN: LITTLE\n");
  fprintf(fpbov,"CENTERING: nodal\n");
  fprintf(fpbov,"BYTE_OFFSET: %d\n", 0); // was 2*(int)sizeof(int));
  fclose(fpbov);

  // first write a header file in XDMF format, to enable reading by ParaView
  FILE * fpxmf = fopen(BASENAME"/Jacobi.xmf", "w");
  fprintf(fpxmf,"<?xml version=\"1.0\" ?>\n");
  fprintf(fpxmf,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fpxmf,"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n");
  fprintf(fpxmf,"  <Domain>\n");
  fprintf(fpxmf,"    <Grid Name=\"Jacobi Mesh\" GridType=\"Uniform\">\n");
  fprintf(fpxmf,"      <Topology TopologyType=\"3DCORECTMESH\" Dimensions=\"1 %d %d\"/>\n", sim->resolution+2, sim->resolution+2);

  fprintf(fpxmf,"      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
  fprintf(fpxmf,"         <DataItem Name=\"Origin\" NumberType=\"Float\" Dimensions=\"3\" Format=\"XML\">0. 0. 0.</DataItem>\n");
  fprintf(fpxmf,"         <DataItem Name=\"Spacing\" NumberType=\"Float\" Dimensions=\"3\" Format=\"XML\">1. 1. 1.</DataItem>\n");
  fprintf(fpxmf,"      </Geometry>\n");
  fprintf(fpxmf,"      <Attribute Name=\"temperature\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fpxmf,"          <DataItem Dimensions=\"1 %d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\">Jacobi.bin</DataItem>\n", sim->resolution+2, sim->resolution+2);
  fprintf(fpxmf,"      </Attribute>\n");
  fprintf(fpxmf,"    </Grid>\n");
  fprintf(fpxmf,"  </Domain>\n");
  fprintf(fpxmf,"</Xdmf>\n");
  fclose(fpxmf);
  }

  // second write the result file in binary
  MPIIOWriteData(fname, sim);

  MPI_Type_free(&rowtype);
  MPI_Type_free(&coltype);

  MPI_Barrier(sim->topocomm);
}
