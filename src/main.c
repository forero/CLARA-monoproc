#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "struct.h"
#include "parse.h"
#include "mtwist.h"
#include "lyman.h"
#include "propagation.h"

int main(int argc, char **argv)
{
  
  /* start MPI */
  ThisProc = 0;
  SizeProc = 1;


  /*initialize Mersene Twister Random Generator*/
  mts_seed32(&(RND_MT_State), All.RandomSeed + ThisProc);


  /* read the simulation setup */
  ReadParameters(argv[1]);


  /*Make basic tests*/
  if(All.TestParallelVel || All.TestParallelVelFast || All.TestFirstScatter ||
      All.TestRND || All.TestPerpVel){
      TestAll();
  }

  /*Make science tests*/
  if(All.NeufeldSlab || All.NeufeldCube || All.ExpandingSphere || All.SimulationCube){
      PropagateAll();
  }


  return 0;    
}
