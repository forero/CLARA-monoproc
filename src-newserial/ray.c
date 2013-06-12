#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "myrand.h"
#include "struct.h"
#include "ray.h"
#include "RND_lyman.h"
#include "propagation.h"

void PhotonFree(lyman_RT_photons *Ph){

  free(Ph->x_out);Ph->x_out=NULL;
  free(Ph->DirX);Ph->DirX=NULL;
  free(Ph->DirY);Ph->DirY=NULL;
  free(Ph->DirZ);Ph->DirZ=NULL;
  free(Ph->PosX);Ph->PosX=NULL;
  free(Ph->PosY);Ph->PosY=NULL;
  free(Ph->PosZ);Ph->PosZ=NULL;
  free(Ph->Intensity);Ph->Intensity=NULL;
  free(Ph->Cell);Ph->Cell=NULL;
  free(Ph->Active);Ph->Active=NULL;
  free(Ph->ScatterHI);Ph->ScatterHI=NULL;
  free(Ph->ScatterDust);Ph->ScatterDust=NULL;
  free(Ph->Wrong);Ph->Wrong=NULL;

  free(Ph);
}

void PhotonListInitialize(lyman_RT_photons *Ph){
    int i;
    double theta, phi;

    for(i=0;i<Ph->N_photons;i++){
      RND_spherical(&(Ph->DirX[i]),
		    &(Ph->DirY[i]),
		    &(Ph->DirZ[i]));

	Ph->PosX[i] = 0.0;
	Ph->PosY[i] = 0.0;
	Ph->PosZ[i] = 0.0;

	if(All.HomogeneousInit){
	  if(All.NeufeldCube){
	    Ph->PosX[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosY[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosZ[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	  }
	  
	  if(All.ExpandingSphere){
	    Ph->PosX[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosY[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosZ[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    do{
	      Ph->PosX[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	      Ph->PosY[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	      Ph->PosZ[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    }while(PropagateIsInside(Ph->PosX[i], Ph->PosY[i], Ph->PosZ[i])==0);	    
	  }	  

	  if(All.NeufeldSlab){
	    Ph->PosX[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosY[i] = 0.0;
	    Ph->PosZ[i] = 0.0;
	  }
	}
    }

    if(All.SimulationCube){
	for(i=0;i<Ph->N_photons;i++){
	    theta = acos(RandFloatUnit());
	    phi = 2.0*PI*RandFloatUnit();
	    Ph->DirX[i] = sin(theta)*cos(phi);
	    Ph->DirY[i] = sin(theta)*sin(phi);
	    Ph->DirZ[i] = cos(theta);
	    Ph->PosX[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosY[i] = 2.0*(RandFloatUnit()-0.5)*All.SlabLength;
	    Ph->PosZ[i] = All.SlabLength/All.Tau;
	}
    }

}

lyman_RT_photons * PhotonListCreate(int N_packages){
  lyman_RT_photons *Ph;
  int i;
  if(N_packages<0){
    fprintf(stderr, "In Ray Create the number of packages is negative\n");    
    exit(0);
  }
  if(!(Ph = malloc(sizeof(lyman_RT_photons)))){
    fprintf(stderr, "Problem in creation of PhotontList\n");
    exit(0);
  }

  /*allocate all the arrays*/
  if(!(Ph->x_out = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (x_out)\n");
    exit(0);
  }
  if(!(Ph->DirX = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Dir)\n");
    exit(0);
  }
  if(!(Ph->DirY = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Dir)\n");
    exit(0);
  }
  if(!(Ph->DirZ = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Dir)\n");
    exit(0);
  }
  if(!(Ph->PosX = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Pos)\n");
    exit(0);
  }
  if(!(Ph->PosY = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Pos)\n");
    exit(0);
  }
  if(!(Ph->PosZ = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Pos)\n");
    exit(0);
  }

  if(!(Ph->Cell = malloc(3*N_packages*sizeof(int)))){
    fprintf(stderr, "Problem in creation of PhotontList (Cell)\n");
    exit(0);
  }


  if(!(Ph->Active = malloc(N_packages*sizeof(int)))){
    fprintf(stderr, "Problem in creation of PhotontList (Active)\n");
    exit(0);
  }

  if(!(Ph->ScatterHI = malloc(N_packages*sizeof(int)))){
    fprintf(stderr, "Problem in creation of PhotontList (ScatterHI)\n");
    exit(0);
  }

  if(!(Ph->ScatterDust = malloc(N_packages*sizeof(int)))){
    fprintf(stderr, "Problem in creation of PhotontList (ScatterDust)\n");
    exit(0);
  }

  if(!(Ph->Wrong = malloc(N_packages*sizeof(int)))){
    fprintf(stderr, "Problem in creation of PhotontList (Wrong)\n");
    exit(0);
  }

  if(!(Ph->Intensity = malloc(N_packages*sizeof(double)))){
    fprintf(stderr, "Problem in creation of PhotontList (Intensity)\n");
    exit(0);
  }

  
  Ph->N_grid_x = 0;
  Ph->N_grid_y = 0;
  Ph->N_grid_z = 0;
  Ph->N_photons = N_packages;

  for(i=0;i<N_packages;i++){
      Ph->x_out[i] = All.InputFrequency;
      Ph->Cell[3*i + 0] = -1;
      Ph->Cell[3*i + 1] = -1;
      Ph->Cell[3*i + 2] = -1;
      Ph->Active[i] = ACTIVE;
      Ph->ScatterHI[i] = 0;
      Ph->ScatterDust[i] = 0;
      Ph->Wrong[i]     = 0;
      Ph->DirX[i]     = 0.0;
      Ph->DirY[i]     = 0.0;
      Ph->DirZ[i]     = 0.0;
      Ph->PosX[i]     = 0.0;
      Ph->PosY[i]     = 0.0;
      Ph->PosZ[i]     = 0.0;
      Ph->Intensity[i] = 1.0;
  }

#ifdef DEBUG
  fprintf(stdout,"Finished the allocation and initialization of %d photon packages\n", N_packages);
#endif

  return Ph;
}
