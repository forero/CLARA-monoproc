#ifndef PROP_H
#define PROP_H
void PropagateAllSetup(void);
int PropagateStep(double *x, double *k_in_photon, double *r_travel, double a, double n_HI);
int PropagatePackage(double *PosX, double *PosY, double *PosZ, 
		     double *DirX, double *DirY, double *DirZ, 
		     double *x_in, int *status);
void PropagateAll(void);
void TestAll(void);
int PropagateIsInside(double PosX, double PosY, double PosZ);
void PropagateGetBulkVel(double *BulkVel, double *Pos);
void PropagateGetTemperature(double *Temp, double *Pos);
void PropagateGetNumberDensity(double *n_HI, double *Pos);
void TestFirstScatter(double x, double a);
void PropagateLorentzFreqChange(double *x, double *Dir, 
				double *BulkVel, double ThermalVel, int sign);
void PropagateLorentzDirChange(double *Dir, double *BulkVel, int sign);
#endif
