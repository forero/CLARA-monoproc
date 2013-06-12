#ifndef RND_LYMAN_H
#define RND_LYMAN_H
void RND_spherical(double *vec_x, double *vec_y, double *vec_z);
void RND_pair(double *r_1, double *r_2);
void RND_lyman_perp_vel_test(void);
void RND_lyman_parallel_vel_test(double x, double a);
void RND_lyman_parallel_vel_fast_test(double x, double a);
double RND_lyman_parallel_vel(double x, double a);
double RND_lyman_parallel_vel_fast(double my_x, double a);
void RND_lyman_perp_vel(double *u_1, double *u_2);
void RND_lyman_atom(double *v_parallel, double *v_perp_1, double *v_perp_2, double *DirPhotonX, double *DirPhotonY, double *DirPhotonZ, 
		    double *x, double a, double g_recoil, int n_points);
#endif 

