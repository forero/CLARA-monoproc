#ifndef PROP_H
#define PROP_H
void PropagateAllSetup(void);
int PropagatePackage(double *PosX, double *PosY, double *PosZ, 
		     double *DirX, double *DirY, double *DirZ, int *n_scatt, 
		     double *x_in, int *status, int n_points);
void PropagateAll(void);
void TestAll(void);
int PropagateIsInside(double PosX, double PosY, double PosZ);
void PropagateGetBulkVel(double *BulkVel, double *Pos);
void PropagateGetTemperature(double *Temp, double *Pos);
void PropagateGetNumberDensity(double *n_HI, double *Pos);
//void TestFirstScatter(double x, double a);
void PropagateLorentzFreqChange(double *x, double *Dir, 
				double *BulkVel, double ThermalVel, int sign);
void PropagateLorentzDirChange(double *Dir, double *BulkVel, int sign);
int count_active(int *status, int n_points);
int PropagateStep(double *v_parallel, double *v_perp_1, double *v_perp_2, double *x_in, double *x_out, double *k_in_x, double *k_in_y, double *k_in_z, double *r_travel, int *status, double a, double n_HI, int n_points);
#endif
