#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "struct.h"
#include "lyman.h"
#include "propagation.h"
#include "ray.h"
#include "RND_lyman.h"
#include "io.h"
#include "myrand.h"

int PropagateIsInside(double PosX, double PosY, double PosZ)
/*depending on the geometrical consideration of the problem at hand,
  decides if the photon is still inside*/
{
    int is_in;
    double radius;
    is_in = 1;

    if(All.SimulationCube){
	if(fabs(PosX)<All.SlabLength && 
	   fabs(PosY)<All.SlabLength &&
	   PosZ < (All.SlabLength*2.0) && PosZ>0.0){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.NeufeldSlab){
	if(fabs(PosX)<All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.NeufeldCube){
	if(fabs(PosX)<All.SlabLength && 
	   fabs(PosY)<All.SlabLength &&
	   fabs(PosZ)<All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.ExpandingSphere){
	radius = PosX*PosX + PosY*PosY + PosZ*PosZ;
	radius = sqrt(radius);
	if(radius< All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }
    
    return is_in;
}

void PropagateAllSetup(void)
/*makes the general problem setup from the input values*/
{
  double a, nu_doppler, ly_sigma_0, column_HI, v_thermal, ly_sigma_v, x_new;
  
  nu_doppler = CONSTANT_NU_DOPPLER*sqrt(All.Temperature/10000.0); /* in s^-1 */
  a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
  ly_sigma_0 =  LyaCrossSection(0,a);
  
  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
  x_new = 200.0*1.0E5/v_thermal;
  ly_sigma_v = LyaCrossSection(x_new,a);

  column_HI = All.Tau / ly_sigma_0;
  All.SlabLength = All.Tau/(All.NumberDensityHI*ly_sigma_0);
  
#ifdef DEBUG
  fprintf(stdout, "cross section for delta_v=200km/s %e\n", ly_sigma_v);
  fprintf(stdout, "cross section %e\n", ly_sigma_0);
  fprintf(stdout, "ColumnDensinty %e\n", column_HI);
  fprintf(stdout, "Temp %g\n", All.Temperature);
  fprintf(stdout, "a %g\n", a);
  fprintf(stdout, "ly_sigma_0 %g\n", ly_sigma_0);
  fprintf(stdout, "Slab Length [cm] %g\n", All.SlabLength);
  fprintf(stdout, "maximum at (a*tau)^1/3 %g\n", pow(a*All.Tau,0.333));
  fprintf(stdout, "scatterings %g\n", 1.612*All.Tau);
#endif

}

int PropagateStep(double *x_in, double *x_out, double *k_in_x, double *k_in_y, double *k_in_z, double *r_travel, int *status, double a, double n_HI, int n_points) 
/* Calculates: 
   1. new frequency 
   2. new propagation direction, 
   3. distance to travel to next interaction

   Note:
   At this point the frequency must be already the comoving one (not lab frame!)
   If the photon is absorbed by dust, then the status is set to -1;
*/
{
    double lya_sigma;
    double tau;
    double r;
    double N;
    double a_in;
    int i;
    double HIInteractionProb; /*probability of interacting with HI*/
    double dust_sigma;
    double rand_interaction;
    double rand_absorption;
    double nu_doppler;
    double g_recoil;
    double temperature;
    int i_photon;
    
    /*Initalizations*/
    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    temperature = (nu_doppler/CONSTANT_NU_DOPPLER)*(nu_doppler/CONSTANT_NU_DOPPLER)*10000.0;
    g_recoil = PLANCK*nu_doppler/(2.0*BOLTZMANN*temperature);

    /* generates the atom velocity the comoving frame*/
    RND_lyman_atom(k_in_x, k_in_y, k_in_z, x_out, a, g_recoil, n_points);
    
    for(i_photon=0;i_photon<n_points;i_photon++){
      if(status[i_photon]==ACTIVE){

	HIInteractionProb = 1.0;
	rand_absorption = 0.0;
	
	/*get a random number to decide if the thing interact with HI*/
	rand_interaction  = RandFloatUnit();      
	N    = n_HI;

	/* basic status at the interacion point*/
	lya_sigma  = LyaCrossSection(x_in[i_photon], a);
	tau        = LyaTau();
	r          = tau/(lya_sigma*N); /*the distance is in cm*/

	if(All.UseDust){
	  dust_sigma = PI*All.GrainSize*All.GrainSize;
	  All.NumberDensityDust = All.TauDust/(All.SlabLength*dust_sigma);
	  /*probability of interacting with HI*/
	  HIInteractionProb = 
	    All.NumberDensityHI*lya_sigma/(All.NumberDensityHI*lya_sigma + All.NumberDensityDust*dust_sigma);
	  /*	fprintf(stdout, "Interaction prob %e \n", HIInteractionProb);*/
	  /*get a random number to decide if the thing will be absorbed*/
	  rand_absorption = RandFloatUnit();	  
	  /* The traveled distance is modifyed by the presence of dust*/
	  r          = tau/(lya_sigma*N + dust_sigma*All.NumberDensityDust); /*the distance is in cm*/
	}

	r_travel[i_photon] = r;

	if((rand_interaction > HIInteractionProb) && All.UseDust && (rand_absorption < All.DustAbsorptionProb)){
	  r_travel[i_photon] = 0.0;
	  status[i_photon] = ABSORBED;
	  x_out[i_photon] = x_in[i_photon];
	}
      }
    }
    return 0;
}

void PropagateGetBulkVel(double *BulkVel, double *Pos){
    int i;

    for(i=0;i<3;i++){
	BulkVel[i] = 0.0;
    }

    if(All.ExpandingSphere){
	for(i=0;i<3;i++){
	    BulkVel[i] = (Pos[i]/All.SlabLength)*(All.VmaxSphere*1.0e5); /*in cm/s*/
	}
    }
}

void PropagateGetTemperature(double *Temp, double *Pos){
    *Temp = All.Temperature;
}

void PropagateGetNumberDensity(double *n_HI, double *Pos){
    *n_HI = All.NumberDensityHI;
}


void PropagateLorentzDirChange(double *Dir, double *BulkVel, int sign){
    int i;
    for(i=0;i<3;i++){
	Dir[i] = Dir[i]*(1.0 + sign*(BulkVel[i]/C_LIGHT)); 
    }
    
    /*renormalize, I just want to know the direction change, bust still be normalized*/
    normalize(Dir);
}

void PropagateLorentzFreqChange(double *x, double *Dir, 
			    double *BulkVel, double ThermalVel, int sign){
    double lorentz_factor;
    lorentz_factor = (BulkVel[0]*Dir[0] + BulkVel[1]*Dir[1] + BulkVel[2]*Dir[2])/ThermalVel;
    
    *x  = *x + sign*lorentz_factor;
}

int count_active(int *status, int n_points){
  int i, n_active;
  n_active = 0;
  for(i=0;i<n_points;i++){
    if(status[i]==ACTIVE){
      n_active++;
    }
  }
  return n_active;
}

int PropagatePackage(double *PosX, double *PosY, double *PosZ, 
		     double *DirX, double *DirY, double *DirZ, int *n_scatt, 
		     double *x_in, int *status, int n_points){
    int i;
    double Pos[3];
    double Dir[3];
    double r_travel, x, x_comoving;
    double a, nu_doppler, n_HI, v_thermal, BulkVel[3], temperature;
    int n_iter;
    int stat;
    float last_x;
    int  n_active;
    int n_global_scatt;
    double *x_aux_in;
    double *x_aux_out;
    double *r_travel_aux;
    n_iter=0;
    n_global_scatt=0;

    /* auxiliary variables */
    if(!(x_aux_in = malloc(sizeof(double) * n_points))){
      fprintf(stderr, "Problem with x_aux allocation\n");
      exit(1);
    }

    if(!(x_aux_out = malloc(sizeof(double) * n_points))){
      fprintf(stderr, "Problem with x_aux allocation\n");
      exit(1);
    }

    if(!(r_travel_aux = malloc(sizeof(double) * n_points))){
      fprintf(stderr, "Problem with r_travel_aux allocation\n");
      exit(1);
    }


    /* count the number of active photons */
    n_active = count_active(status, n_points);
#ifdef DEBUG
    fprintf(stdout, "Active photons: %d", n_active);
    fflush(stdout);
#endif



    /* difuse the photon in space and frequency until it gets out */
    while(n_active>0){

      for(i=0;i<n_points;i++){
	/*Make the initialization*/
	Pos[0] = PosX[i];
	Pos[1] = PosY[i];
	Pos[2] = PosZ[i];
	Dir[0] = DirX[i];
	Dir[1] = DirY[i];
	Dir[2] = DirZ[i];
	stat = status[i];
	if(stat==ACTIVE){
	  x_aux_in[i] = x_in[i];    
	  x_aux_out[i] = x_in[i];    
	}else{
	  x_aux_in[i] = 0.0;
	  x_aux_out[i] = 0.0;
	}
	
	/*If the photon is stil inside and active, update its direction of propagation and frequency 
	  to  gas rest frame */
	if(PropagateIsInside(Pos[0],Pos[1],Pos[2])&&(stat==ACTIVE)){
	  /* get the temperature at this point*/
	  PropagateGetTemperature(&temperature, Pos);
	  
	  /* get the number density at this point*/
	  PropagateGetNumberDensity(&n_HI, Pos);
	  
	  /*get the bulk velocity of the fluid at this point*/
	  PropagateGetBulkVel(BulkVel, Pos);
	  
	  /*Get the thermal velocity and doppler broadening*/
	  nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
	  a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
	  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
	  
	  /*change the value of the frequency to one comoving with the fluid*/
	  PropagateLorentzFreqChange(&(x_aux_in[i]), Dir, BulkVel, v_thermal, -1); 
	  
	  /*change the direction of the photon to the fluid frame*/
	  PropagateLorentzDirChange(&(Dir[0]), BulkVel, -1);
	}
	DirX[i] = Dir[0];
	DirY[i] = Dir[1];
	DirZ[i] = Dir[2];
      }
	  
      PropagateGetNumberDensity(&n_HI, Pos);
      PropagateGetTemperature(&temperature, Pos);
      nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
      a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
      v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/

      /*Change the frequency and the Propagation direction, find the displacement*/	
      PropagateStep(x_aux_in, x_aux_out, DirX, DirY, DirZ, r_travel_aux, status, a, n_HI, n_points);	    	
            
      for(i=0;i<n_points;i++){
	Pos[0] = PosX[i];
	Pos[1] = PosY[i];
	Pos[2] = PosZ[i];	
	Dir[0] = DirX[i];
	Dir[1] = DirY[i];
	Dir[2] = DirZ[i];
	stat = status[i];	

	if(PropagateIsInside(Pos[0],Pos[1],Pos[2])&&(stat==ACTIVE)){
	  /* get the temperature at this point*/
	  PropagateGetTemperature(&temperature, Pos);
	  
	  /* get the number density at this point*/
	  PropagateGetNumberDensity(&n_HI, Pos);
	  
	  /*get the bulk velocity of the fluid at this point*/
	  PropagateGetBulkVel(BulkVel, Pos);

	  /*Get the thermal velocity and doppler broadening*/
	  nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
	  a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
	  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
	  
	  /*Change the new direction to the lab frame value*/
	  PropagateLorentzDirChange(Dir, BulkVel, 1);
	  
	  /*Change the frequency comoving to the lab frame value*/
	  PropagateLorentzFreqChange(&(x_aux_out[i]), Dir, BulkVel, v_thermal, 1); 
	  
	  /*Update the position*/
	  Pos[0] += r_travel_aux[i] * Dir[0];	    
	  Pos[1] += r_travel_aux[i] * Dir[1];	    
	  Pos[2] += r_travel_aux[i] * Dir[2];    
	  n_scatt[i]++;

	  /*update the final status of the photon, just to know if it was absorbed, or what*/
	  if(!(PropagateIsInside(Pos[0],Pos[1],Pos[2]))){
	    status[i] = OUT_OF_BOX;
	  }

	  PosX[i] = Pos[0];
	  PosY[i] = Pos[1];
	  PosZ[i] = Pos[2];
	  DirX[i] = Dir[0];
	  DirY[i] = Dir[1];
	  DirZ[i] = Dir[2];	  
	  x_in[i] = x_aux_out[i];
	}

	if(i==0){
	  n_global_scatt++;
	}
      }
      n_active = count_active(status, n_points);
#ifdef DEBUG      
      fprintf(stdout, "Active photons: %d, n_scatt %d \n", n_active, n_global_scatt);
#endif
    }

    free(x_aux_in);
    free(x_aux_out);
    free(r_travel_aux);
    return 0;    
}

    
void PropagateAll(void)
{
    int n_packages;
    int i;
    lyman_RT_photons *Ph;
    int n_iter;
    char FileName[MAX_FILENAME_SIZE];
    int status;

    /*update some geometrical values*/
    PropagateAllSetup();
    
    /*get the number of photons to propagate*/
    n_packages = (int)((All.TotalLuminosity/SizeProc)/All.LuminosityPerPackage);
#ifdef DEBUG
    fprintf(stdout, "%d packages to propagate\n", n_packages);
#endif 

    /*create the packages*/
    Ph = PhotonListCreate(n_packages);
    
    /*initialize the packages*/
    PhotonListInitialize(Ph);
    
    if(All.OutputInitList){
	sprintf(FileName, "%s/%s_in.proc.%d.ascii", All.OutputDir, All.OutputFile, ThisProc);
	SavePhotonListAscii(FileName, Ph);
    }
    
    if(All.OutputFinalList){
	sprintf(FileName, "%s/%s_out.proc.%d.ascii", All.OutputDir, All.OutputFile, ThisProc);
	OpenPhotonListAscii(FileName, Ph);
    }
    
    /*propagate each package*/

    status = PropagatePackage(Ph->PosX, Ph->PosY, Ph->PosZ,
			      Ph->DirX, Ph->DirY, Ph->DirZ,
			      Ph->ScatterHI, Ph->x_out, Ph->Active, 
			      n_packages);
    

    for(i=0;i<n_packages;i++){
	if(All.OutputFinalList){
	  sprintf(FileName, "%s/%s_out.proc.%d.ascii", All.OutputDir, All.OutputFile, ThisProc);
	  AppendPhotonListAscii(FileName, Ph, i);    
	}    
    }
    fprintf(stdout, "finished writing (prc %d)\n", ThisProc);
    /*free the memory*/
    PhotonFree(Ph);
}

void TestRND(void)
{
    double x_in;
    double r_travel;
    double k_in_photon[3];
    double n_HI;
    int status;
    int i;
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;

    sprintf(FileName, "%s/%s_RND.proc.%d", All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "TestRND problem opening file %s\n", FileName);
	exit(0);
    }

    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, 0.0, 0.0);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	r_travel = RandFloatUnit();
	fprintf(out, "%e\n", r_travel);
    }

    fclose(out);
} 

void TestAll(void){
    if(All.TestParallelVel){
	RND_lyman_parallel_vel_test(All.Test_x,All.Test_a);
    }

    if(All.TestParallelVelFast){
	RND_lyman_parallel_vel_fast_test(All.Test_x,All.Test_a);
    }    
    /*
    if(All.TestFirstScatter){
	TestFirstScatter(All.Test_x, All.Test_a);
    }
    */

    if(All.TestRND){
	TestRND();
    }

    if(All.TestPerpVel){
       RND_lyman_perp_vel_test();
    }
}
