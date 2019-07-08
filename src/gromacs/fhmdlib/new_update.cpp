#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "data_structures.h"
#include "interpolation.h"
#include "sfunction.h"
#include "macro.h"

// new code start
#include <fstream> // to write out to files
#include <iostream> // to write out to files

std::ofstream ofs;

int start_step = 500;

double norm(double vect[3]) {
		return sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
}

bool writing_condition(int step);

class StepWriter {
    const char *name;
    const int width;
    const int stepnum;
    bool started_writing = false;
    char filename[100] = "";
public:
    StepWriter(const char *inname, int width, int step) : name(inname), width(width), stepnum(step) {

//    	char str[100];
//		sprintf(str, "%02d", file_num++);
//
//		strcpy(filename, str);
//		strcat(filename, "_");
		strcat(filename, name);
		strcat(filename, ".csv");

		if(stepnum == 0) {
			write_header();
		}


    	if(writing_condition(step)) {

//			printf("Filename: %s\n Step: %d\n\n", filename, step);
			ofs.open(filename, std::ofstream::out | std::ofstream::app);
			// printf("step IS %d", step);
			ofs << step;
			ofs.close();
    	}

    };
    void write_header();
    void write(double value);
    ~StepWriter();
};

void StepWriter::write_header() {
    ofs.open(filename, std::ofstream::out);
    ofs << "Step";
    for(int i = 0; i < width; i++)
    {
        ofs << ",";
        ofs << i;
    }
    ofs << "\n";
    ofs.close();
}

void StepWriter::write(double value) {
	if(not writing_condition(stepnum)) return;
    ofs.open(filename, std::ofstream::out | std::ofstream::app);
    ofs << ",";
    ofs << value;
    ofs.close();
}

StepWriter::~StepWriter() {
	if(not writing_condition(stepnum)) return;
    ofs.open(filename, std::ofstream::out | std::ofstream::app);
    ofs << "\n";
    ofs.close();
}

bool writing_condition(int step)
{
	if (step >= start_step)
		return true;
	else
		return false;
}

static const bool output_gromacs_parts = true;
static const bool alpha_beta_calculations = true;

// new code end


void fhmd_do_update_md(int start, int nrend,
                       double dt, int nstpcouple,
                       t_grp_tcstat *tcstat,
                       double nh_vxi[],
                       gmx_bool bNEMD, t_grp_acc *gstat, rvec accel[],
                       ivec nFreeze[],
                       real invmass[],
                       unsigned short ptype[], unsigned short cFREEZE[],
                       unsigned short cACC[], unsigned short cTC[],
                       rvec x[], rvec xprime[], rvec v[],
                       rvec f[], matrix M,
                       gmx_bool bNH, gmx_bool bPR, t_commrec *cr, FHMD *fh)
{
    double imass, w_dt;
    int    gf = 0, ga = 0, gt = 0;
    rvec   vrel;
    real   vn, vv, va, vb, vnrel;
    real   lg, vxi = 0, u;
    int    n, d;

    /* FHMD variables */
    FH_arrays   *arr = fh->arr;
    int          ind;
    double       invro_dt;
    double       S = fh->S;
    double       gamma_u, gamma_x;
    int          nbr[8];
    dvec         xi;
    dvec         f_fh, u_fh, alpha_term, beta_term, grad_ro;
    const double g_eps = 1e-10;

    // new code start
    double beta[nrend - start]; // the array [0, 7999] = beta[8000]
    double alpha[nrend - start]; // the array [0, 7999] = beta[8000]

//    printf("Counting %d atoms....\n", nrend - start);


    StepWriter beta_writer = StepWriter("beta_values", nrend-start, fh->step_MD);
    StepWriter alpha_writer = StepWriter("alpha_values", nrend-start, fh->step_MD);
    StepWriter a_coef_writer = StepWriter("a_coef", nrend-start, fh->step_MD);
    StepWriter b_coef_writer = StepWriter("b_coef", nrend-start, fh->step_MD);

	// vn, f[n][d], lg, w_dt,

	StepWriter vn_writer = StepWriter("vn_writer", nrend-start, fh->step_MD);
	StepWriter fnd_writer = StepWriter("fnd_writer", nrend-start, fh->step_MD);
	StepWriter lg_writer = StepWriter("lg_writer", nrend-start, fh->step_MD);
	StepWriter w_dt_writer = StepWriter("w_dt_writer", nrend-start, fh->step_MD);


    StepWriter f_fh_norm = StepWriter("f_fh_norm", nrend-start, fh->step_MD);
    StepWriter alpha_n = StepWriter("alpha_n", nrend-start, fh->step_MD);
    StepWriter aplha_term_norm = StepWriter("aplha_term_norm", nrend-start, fh->step_MD);
    StepWriter beta_n = StepWriter("beta_n", nrend-start, fh->step_MD);
    StepWriter beta_term_norm = StepWriter("beta_term_norm", nrend-start, fh->step_MD);
    StepWriter u_fh_norm = StepWriter("u_fh_norm", nrend-start, fh->step_MD);
    StepWriter grad_ro_writer = StepWriter("grad_ro", nrend-start, fh->step_MD);

    // adding new logic

    for (int ind = 0; ind < fh->Ntot; ind++)
    {
        arr[ind].numerator[0] = 0;
        arr[ind].numerator[1] = 0;
        arr[ind].numerator[2] = 0;
        arr[ind].numerator[3] = 0;
        arr[ind].denominator[0] = 0;
        arr[ind].denominator[1] = 0;
    }


	rvec component;

	// begin to update the statistics
	for (int n = start; n < nrend; n++)
	{

		ind      = fh->ind[n]; // index of cell containing atom n


		// interpolation of terms for cells
		{
			trilinear_find_neighbours(x[n], n, xi, nbr, fh);

			if(fh->scheme == Two_Way)
				trilinear_interpolation(f_fh,   xi, INTERPOLATE(f_fh));
			else
				clear_dvec(f_fh);

			trilinear_interpolation(u_fh,       xi, INTERPOLATE(u_fh));
			trilinear_interpolation(alpha_term, xi, INTERPOLATE(alpha_term));
			trilinear_interpolation(beta_term,  xi, INTERPOLATE(beta_term));
			trilinear_interpolation(grad_ro,    xi, INTERPOLATE(grad_ro));

		}


		if(fh->S_function == moving_sphere)
			S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
		else if(fh->S_function == fixed_sphere)
			S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

		// calculations

		// compute values under that comprise the sum

		// numerator 0
		for (d = 0; d < DIM; d++)
		{
			component[d] = -(1.0/arr[ind].ppm) * f[n][d] * S * (u_fh[d] - v[n][d]);
		}
		arr[ind].numerator[0] += SUM(component);


		// numerator 1
		for (d = 0; d < DIM; d++)
		{
			component[d] = -(1.0/arr[ind].ppm) * v[n][d] * (S*f[n][d] - S*f_fh[d]);
		}
		arr[ind].numerator[1] += SUM(component);

		// numerator 2, without alpha
		for (d = 0; d < DIM; d++)
		{
			component[d] = (1.0/arr[ind].ppm) * (v[n][d]/invmass[n]) * arr[ind].inv_ro * (arr[ind].alpha_term[d]);
		}
		arr[ind].numerator[2] += SUM(component);

		// numerator 3, without alpha
		for (d = 0; d < DIM; d++)
		{
			component[d] = -(1.0/arr[ind].ppm) * f[n][d] * S * (1 - S) * arr[ind].grad_ro[d] * arr[ind].inv_ro;
		}
		arr[ind].numerator[3] += SUM(component);

		double sum_uro_without_p[DIM];

		rvec pr_v;

		for (d = 0; d < DIM; d++)
		{
			sum_uro_without_p[d] = 0;
			for(int k = start; k < nrend; k++)
			{
				if ((k == n) or (ind != fh->ind[k])) continue;
				sum_uro_without_p[d] += v[k][d] * (1.0/invmass[k]) * fh->grid.ivol[ind];
			}

			pr_v[d] = arr[ind].u_fh[d] * arr[ind].ro_fh - sum_uro_without_p[d];

			component[d] = -(1.0/arr[ind].ppm) * (v[n][d]/invmass[n]) * S * (1 - S) * pr_v[d] * arr[ind].inv_ro;
		}

		arr[ind].denominator[1] += SUM(component);
	} // calculations of components


	for (int n = start; n < nrend; n++) // 0, 1, 2, ... 7999
	{
		ind = fh->ind[n]; // index of cell containing atom n

		for (d = 0; d < DIM; d++)
		{
			component[d] = S * (1 - S) * (1.0/invmass[n]) * fh->grid.ivol[ind] / arr[ind].ro_md;
		}

		arr[ind].denominator[0] = SUM(component);


		double a_coef = (arr[ind].numerator[2] + arr[ind].numerator[3])/(arr[ind].denominator[0] + arr[ind].denominator[1]);
		double b_coef = (arr[ind].numerator[0] + arr[ind].numerator[1])/(arr[ind].denominator[0] + arr[ind].denominator[1]);

//		float alph; // = 10;
//
//
//		if (a_coef > 0 and b_coef > 0)
//		{
//			// alpha[n-start] = 50;
//		}
//		else if (a_coef < 0 and b_coef < 0)
//		{
//			alpha[n-start] = 0;
//		}
//		else
//		{
//			alpha[n-start] = -a_coef*b_coef/(a_coef*a_coef + 1);
//		}
//
//
//		beta[n-start] = a_coef * alpha[n-start] + b_coef;

		double alpha = 20;
		beta[n-start] = a_coef * alpha + b_coef;


		beta_writer.write(beta[n-start]);
		// alpha_writer.write(alpha[n-start]);
		alpha_writer.write(alpha);
		a_coef_writer.write(a_coef);
		b_coef_writer.write(b_coef);

	} // calculations of first_bottom_v and beta coef

    // new code end

    if (bNH || bPR)
    {
        /* Update with coupling to extended ensembles, used for
         * Nose-Hoover and Parrinello-Rahman coupling
         * Nose-Hoover uses the reversible leap-frog integrator from
         * Holian et al. Phys Rev E 52(3) : 2338, 1995
         */

        /* FHMD Error */
        printf(MAKE_RED "\nFHMD: ERROR: FH-MD coupling doesn't support Nose-Hoover and Parrinello-Rahman\n" RESET_COLOR "\n");
        exit(11);

    }
    else if (cFREEZE != NULL ||
             nFreeze[0][XX] || nFreeze[0][YY] || nFreeze[0][ZZ] ||
             bNEMD)
    {
        /* Update with Berendsen/v-rescale coupling and freeze or NEMD */

        /* FHMD Error */
        printf(MAKE_RED "\nFHMD: ERROR: FH-MD coupling doesn't support freeze or NEMD\n" RESET_COLOR "\n");
        exit(12);

    }
    else
    {
        /* Plain update with Berendsen/v-rescale coupling */
        for (n = start; n < nrend; n++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell))
            {
                w_dt     = invmass[n]*dt;
                ind      = fh->ind[n];
                invro_dt = arr[ind].inv_ro*dt;

                trilinear_find_neighbours(x[n], n, xi, nbr, fh);

                if(fh->scheme == Two_Way)
                    trilinear_interpolation(f_fh,   xi, INTERPOLATE(f_fh));
                else
                    clear_dvec(f_fh);
                trilinear_interpolation(u_fh,       xi, INTERPOLATE(u_fh));
                trilinear_interpolation(alpha_term, xi, INTERPOLATE(alpha_term));
                trilinear_interpolation(beta_term,  xi, INTERPOLATE(beta_term));
                trilinear_interpolation(grad_ro,    xi, INTERPOLATE(grad_ro));

#ifdef FHMD_DEBUG_INTERPOL
                if(!(n % 10000) && !(fh->step_MD % 100))
                    printf("\nStep %d, atom #%d (%g %g %g): %g %g %g\nneighbour cells: %d %d %d %d %d %d %d %d\nrescaled coordinates: %g %g %g\n",
                            fh->step_MD, n, x[n][0], x[n][1], x[n][2], u_fh[0], u_fh[1], u_fh[2],
                            nbr[0], nbr[1], nbr[2], nbr[3], nbr[4], nbr[5], nbr[6], nbr[7], xi[0], xi[1], xi[2]);
#endif

                if(fh->S_function == moving_sphere)
                    S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
                else if(fh->S_function == fixed_sphere)
                    S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

                if (cTC)
                {
                    gt = cTC[n];
                }
                lg = tcstat[gt].lambda;                             // Thermostat

                /* Local thermostat */
                if(fh->S_berendsen >= 0)
                {
                    if(S > fh->S_berendsen) lg = 1;
                }
                else
                {
                    lg = lg*(1 - pow(S, -fh->S_berendsen)) + pow(S, -fh->S_berendsen);
                }

                for (d = 0; d < DIM; d++)
                {
                     /* Pure MD: */
//                     vn           = lg*v[n][d] + f[n][d]*w_dt;
//                     v[n][d]      = vn;
//                     xprime[n][d] = x[n][d] + vn*dt;

                	// vn, f[n][d], lg, w_dt,


                    if(fh->scheme == One_Way)
                    {
                    	// new code start

                    	/*
                    	 * old:
                       		vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
                    		v[n][d]      = vn;
                    		xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;


                    	 */

                    	if(alpha_beta_calculations)
                    	{
                        	// what is affected if we do alpha and beta calculations?
                        	// - grad_ro
                        	// - alpha_u_grad
                        	// - alpha_term

                    		vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha[n] * alpha_term[d] + beta[n] * S*(1 - S)*beta_term[d])*invro_dt;
                    		v[n][d]      = vn;
                    		xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + alpha[n] * S*(1 - S)*grad_ro[d]*invro_dt;

                    	}
                    	else
                    	{
                    		vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + 	 	   alpha_term[d] +			 S*(1 - S)*beta_term[d])*invro_dt;
							v[n][d]      = vn;
							xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + 			  S*(1 - S)*grad_ro[d]*invro_dt;
                    	}
                    	// new code end

                    }
                    else if(fh->scheme == Two_Way)
                    {
                        gamma_u = fh->gamma_u*S*S*S*S*dt*(fh->stat.std_u_fh[d]*fh->stat.std_u_fh[d]/(fh->std_u*fh->std_u) - 1);
                        gamma_x = fh->gamma_x*S*S*S*S*dt*(fh->stat.std_rho_fh/fh->std_rho - 1);

                        if(fabs(gamma_u) < g_eps) gamma_u = g_eps;
                        if(fabs(gamma_x) < g_eps) gamma_x = g_eps;

                        vn = lg*v[n][d]*exp(-gamma_u) + ((1 - S)*f[n][d]*invmass[n] + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*arr[ind].inv_ro)
                                 *(1 - exp(-gamma_u))/gamma_u*dt;

                        v[n][d] = vn;

                        xprime[n][d] = x[n][d] + (1 - S)*vn*(1 - exp(-gamma_u))/gamma_u*dt +
                                           (S*u_fh[d] + S*(1 - S)*grad_ro[d]*arr[ind].inv_ro)*(1 - exp(-gamma_x))/gamma_x*dt;
                    } // if


                }

                // new code start
                f_fh_norm.write(norm(f_fh));
                alpha_n.write(alpha[n]);
                aplha_term_norm.write(norm(alpha_term));
                beta_n.write(beta[n]);
                beta_term_norm.write(norm(beta_term));
                u_fh_norm.write(norm(u_fh));
                grad_ro_writer.write(norm(grad_ro));

                if(output_gromacs_parts == true) {
					vn_writer.write(norm(v[n]));
					fnd_writer.write(norm(f[n]));
					lg_writer.write(lg);
					w_dt_writer.write(w_dt);
                }

                // new code end
            }
            else
            {
                for (d = 0; d < DIM; d++)
                {
                    v[n][d]        = 0.0;
                    xprime[n][d]   = x[n][d];
                }
            }
        }
    }
}
