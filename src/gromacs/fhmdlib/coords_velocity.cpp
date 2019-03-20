#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "data_structures.h"
#include "interpolation.h"
#include "sfunction.h"
#include "macro.h"

double velocity(real lg, rvec v[], int n, int d, double S, rvec f[], dvec f_fh, dvec alpha_term, dvec beta_term, double invro_dt)
{
	double vn;
    vn = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
    
	return vn;
}

double coords(rvec x[], int n, int d, double S, real vn, dvec u_fh, double dt, dvec grad_ro, double invro_dt)
{
    double xprime_n_d;
	xprime_n_d = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;
	
	return xprime_n_d;
}
	