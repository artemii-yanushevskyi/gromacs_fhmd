#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "data_structures.h"
#include "interpolation.h"
#include "sfunction.h"
#include "macro.h"

#include <fstream> // to write out to files
#include <iostream> // to write out to files
using namespace std;

bool writing_condition(int step)
{
	if (step % 50 == 0)
		return true;
	else
		return false;
}

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

    double beta[nrend - start]; // the array [0, 7999] = beta[8000]
    double alpha[nrend - start]; // the array [0, 7999] = beta[8000]

    ofstream ofs;

//    printf("Counting %d atoms....\n", nrend - start);



    if (fh->step_MD == 0)
	{
    	int i;

		ofs.open("log.txt", std::ofstream::out);
		ofs << fh->step_MD;
		ofs << "atoms ";
		ofs << start << nrend;
		ofs << "\n start to output on ";
		ofs << "step\n alpha = " << fh->alpha;

		ofs << "\n";
		ofs.close();

		ofs.open("beta_values.csv", std::ofstream::out);
		ofs << "Step";
		for(i = start; i < nrend; i++)
		{
			ofs << ",";
			ofs << i;
		}
		ofs << "\n";
		ofs.close();

		ofs.open("alpha_values.csv", std::ofstream::out);
		ofs << "Step";
		for(i = start; i < nrend; i++)
		{
			ofs << ",";
			ofs << i;
		}
		ofs << "\n";
		ofs.close();

		ofs.open("a_coef.csv", std::ofstream::out);
		ofs << "Step";
		for(i = start; i < nrend; i++)
		{
			ofs << ",";
			ofs << i;
		}
		ofs << "\n";
		ofs.close();

		ofs.open("b_coef.csv", std::ofstream::out);
		ofs << "Step";
		for(i = start; i < nrend; i++)
		{
			ofs << ",";
			ofs << i;
		}
		ofs << "\n";
		ofs.close();

		ofs.open("atom_cell.csv", std::ofstream::out);
		ofs << "step";
		for(i = start; i < nrend; i++)
		{
			ofs << ",";
			ofs << i;
		}
		ofs << "\n";
		ofs.close();

		ofs.open("numerator.csv", std::ofstream::out);
		ofs << "step";
		ofs << "0,1,2,3,4,0,1";
		ofs << "\n";
		ofs.close();

	}

    if (writing_condition(fh->step_MD))
    {
		ofs.open("a_coef.csv", std::ofstream::out | std::ofstream::app);
		ofs << fh->step_MD;
		ofs.close();

		ofs.open("b_coef.csv", std::ofstream::out | std::ofstream::app);
		ofs << fh->step_MD;
		ofs.close();

		ofs.open("atom_cell.csv", std::ofstream::out | std::ofstream::app);
		ofs << fh->step_MD;
		ofs.close();

		ofs.open("beta_values.csv", std::ofstream::out | std::ofstream::app);
		ofs << fh->step_MD;
		ofs.close();

		ofs.open("alpha_values.csv", std::ofstream::out | std::ofstream::app);
		ofs << fh->step_MD;
		ofs.close();

		ofs.open("numerator.csv", std::ofstream::out | std::ofstream::app);
		ofs << fh->step_MD;
		ofs.close();
    }

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

    double mean = 0, min = 1000, max = -1000;

    if(1)
    {

        rvec component;

        // begin to update the statistics
        for (int n = start; n < nrend; n++)
        {

            ind      = fh->ind[n]; // index of cell containing atom n

    	    if (writing_condition(fh->step_MD))
    	    {
        		ofs.open("atom_cell.csv", std::ofstream::out | std::ofstream::app);
        		ofs << ",";
        		ofs << ind;
        		ofs.close();
    	    }

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

        for (int ind = 0; ind < fh->Ntot; ind++)
        {
        	if (fh->step_MD % 500)
        	{
        		ofs.open("numerator.csv", std::ofstream::out | std::ofstream::app);
        		ofs << "," << arr[ind].numerator[0] << "," << arr[ind].numerator[1] << "," << arr[ind].numerator[2] << "," << arr[ind].numerator[3] << "," << "\n";
        		ofs.close();
        	}

        }

        for (int n = start; n < nrend; n++) // 0, 1, 2, ... 7999
        {
            ind = fh->ind[n]; // index of cell containing atom n

            for (d = 0; d < DIM; d++)
            {
                component[d] = S * (1 - S) * (1.0/invmass[n]) * fh->grid.ivol[ind] / arr[ind].ro_md;
            }

            arr[ind].denominator[0] = SUM(component);

        	if (fh->step_MD % 500)
        	{
        		ofs.open("numerator.csv", std::ofstream::out | std::ofstream::app);
        		ofs << "," << arr[ind].denominator[0] << "," << arr[ind].denominator[1] << "\n";
        		ofs.close();
        	}

            double a_coef = (arr[ind].numerator[2] + arr[ind].numerator[3])/(arr[ind].denominator[0] + arr[ind].denominator[1]);
            double b_coef = (arr[ind].numerator[0] + arr[ind].numerator[1])/(arr[ind].denominator[0] + arr[ind].denominator[1]);

        	float alph; // = 10;


        	if (a_coef > 0 and b_coef > 0)
        	{
        		alpha[n-start] = 50;
        	}
        	else if (a_coef < 0 and b_coef < 0)
        	{
    			alpha[n-start] = 0;
        	}
        	else
        	{
        		alpha[n-start] = -a_coef*b_coef/(a_coef*a_coef + 1);
        	}



//        	{
//        		alpha = 50;
//        		beta[n-start] = 10;
//        	}

            beta[n-start] = a_coef * alpha[n-start] + b_coef;

            if (max < beta[n-start]) max = beta[n-start];
            if (min > beta[n-start]) min = beta[n-start];

            mean += beta[n-start]/(nrend-start);

            if(writing_condition(fh->step_MD))
            {
				ofs.open("a_coef.csv", std::ofstream::out | std::ofstream::app);
				ofs << ",";
				ofs << a_coef;
				ofs.close();

				ofs.open("b_coef.csv", std::ofstream::out | std::ofstream::app);
				ofs << ",";
				ofs << b_coef;
				ofs.close();

				ofs.open("beta_values.csv", std::ofstream::out | std::ofstream::app);
				ofs << ",";
				ofs << beta[n-start];
				ofs.close();

				ofs.open("alpha_values.csv", std::ofstream::out | std::ofstream::app);
				ofs << ",";
				ofs << alpha[n-start];
				ofs.close();
            }

        } // calculations of first_bottom_v and beta coef
    } // if (1)



	ofs.open("log.txt", std::ofstream::out | std::ofstream::app);
	ofs << fh->step_MD << "," << min << "," << mean << "," << max;
	ofs << "\n";
	ofs.close();

    if(writing_condition(fh->step_MD))
    {
		int i;

		ofs.open("atom_cell.csv", std::ofstream::out | std::ofstream::app);
		ofs << "\n";
		ofs.close();

		ofs.open("a_coef.csv", std::ofstream::out | std::ofstream::app);
		ofs << "\n";
		ofs.close();

		ofs.open("b_coef.csv", std::ofstream::out | std::ofstream::app);
		ofs << "\n";
		ofs.close();

		ofs.open("beta_values.csv", std::ofstream::out | std::ofstream::app);
		ofs << "\n";
		ofs.close();

		ofs.open("alpha_values.csv", std::ofstream::out | std::ofstream::app);
		ofs << "\n";
		ofs.close();
    }

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
                     /* vn           = lg*v[n][d] + f[n][d]*w_dt; */
                     /* v[n][d]      = vn; */
                     /* xprime[n][d] = x[n][d] + vn*dt; */

                    if(fh->scheme == One_Way)
                    {
                        // vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt;
//                        vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + beta[n-start]*S*(1 - S)*beta_term[d])*invro_dt;
                        vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha[n-start]*alpha_term[d] + beta[n-start]*S*(1 - S)*beta_term[d])*invro_dt;

                        v[n][d]      = vn;
                        xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + S*(1 - S)*grad_ro[d]*invro_dt;
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


                    }
                }
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
