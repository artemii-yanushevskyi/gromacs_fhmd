#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "data_structures.h"
#include "interpolation.h"
#include "sfunction.h"
#include "macro.h"


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

    printf("start %d, end %d\n", start, nrend);

    double beta[nrend - start];

    if(0)
    {
        double alpha = fh->alpha;

        double first_top[DIM];
        double first_top_v;

        double second_top[DIM];
        double second_top_v;

        double first_bottom[DIM];
        double first_bottom_v;

        double second_bottom[DIM];
        double second_bottom_v;

        // beta is not constant
        // we have n - is particle


        for(int i = 0; i < fh->Ntot; i++)
        {
            arr[i].first_top =        0;
            arr[i].second_top =       0;
            arr[i].second_bottom =    0;
        }


        for (int n = start; n < nrend; n++)
        {
            ind = fh->ind[n]; // index of cell containing atom n

            if(fh->S_function == moving_sphere)
                S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
            else if(fh->S_function == fixed_sphere)
                S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

            double ro_fh = arr[ind].ro_fh;
            double ro_md = arr[ind].ro_md;
            double ro_tilde = arr[ind].ro_fh;
            double m = 1/invmass[n];

            double ppm = arr[ind].ppm;

            rvec F_md;
            rvec mv; // momentum
            rvec u_tilde;
            dvec grad_ro_ss_ro_alpha;


            for(d = 0; d < DIM; d++)
            {
                F_md[d] = f[n][d];
                mv[d] = m*v[n][d]; // only v_md ?
                u_tilde[d] = u_fh[d]; //  fh-velocity
                grad_ro_ss_ro_alpha[d] = S*(1 - S)*arr[ind].grad_ro[d] * arr[ind].inv_ro;
            }

            // compute

            // first top
            for (d = 0; d < DIM; d++)
            {
                first_top[d] = - (1/ppm) * F_md[d] * (S*(u_tilde[d] - v[n][d]) + grad_ro_ss_ro_alpha[d]);
            }

            arr[ind].first_top += SUM(first_top);

            // second top

            for (d = 0; d < DIM; d++)
            {
                second_top[d] = - (1/ppm) * mv[d]/m * (S*F_md[d] - m * alpha_term[d] / ro_md);
                // in mv what's v?
            }

            arr[ind].second_top += SUM(second_top);


            double sum_uro_without_p[DIM];

            rvec pr_v;

            // second bottom
            for (d = 0; d < DIM; d++)
            {
                sum_uro_without_p[d] = 0;
                for(int k = start; k < nrend; k++)
                {
                    if ((k == n) || (ind == fh->ind[k])) continue;
                    sum_uro_without_p[d] += v[k][d] * m * fh->grid.ivol[ind];
                }

                pr_v[d] = u_tilde[d]*ro_tilde - sum_uro_without_p[d];

                second_bottom[d] = - 1/ppm * m*S*(1-S)*pr_v[d]/ro_md;
            }
            arr[ind].second_bottom += SUM(second_bottom);
        }

        for (int n = start; n < nrend; n++)
        {
            w_dt     = invmass[n]*dt;
            ind      = fh->ind[n]; // index of cell containing atom n
            invro_dt = arr[ind].inv_ro*dt;

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

            double ro_fh = arr[ind].ro_fh;
            double ro_md = arr[ind].ro_md;
            double ro_tilde = arr[ind].ro_fh;
            double m = 1/invmass[n];

            double ppm = arr[ind].ppm;

            rvec F_md;
            rvec mv; // momentum
            rvec u_tilde;
            dvec grad_ro_ss_ro_alpha;


            for(d = 0; d < DIM; d++)
            {
                F_md[d] = f[n][d];
                mv[d] = m*v[n][d];
                u_tilde[d] = u_fh[d]; //  fh-velocity
                grad_ro_ss_ro_alpha[d] = S*(1 - S)*arr[ind].grad_ro[d] * arr[ind].inv_ro;
            }



            // // first top
            // for(int p = start; p < nrend; p++)
            // {
            //     for (d = 0; d < DIM; d++)
            //     {
            //         u_tilde[d] = u_fh[d];
            //         first_top[d] = - (1/ppm) * F_md[d] * (S*(u_tilde[d] - v[p][d]) + grad_ro_ss_ro_alpha[d]);
            //     }
            //     first_top_v += SUM(first_top);
            // }
            //
            //
            //
            // // second top
            // for(int p = start; p < nrend; p++)
            // {
            //     for (d = 0; d < DIM; d++)
            //     {
            //         second_top[d] = mv[d]/m * (S*F_md[d] - m * alpha_term[d]);
            //     }
            //     second_top_v += SUM(second_top);
            // }


            // first bottom
            for (d = 0; d < DIM; d++)
            {
                first_bottom[d] = S*(1 - S)*(m*fh->grid.ivol[ind])/ro_md; //
                // ro_p is m_p/V
                // Це те саме ро, що в їхніх формулах
                // Просто виділене для однієї частинки
            }

            first_bottom_v = SUM(first_bottom);

            //
            // double sum_uro_without_p[DIM];
            //
            // rvec pr_v;
            //
            //
            // // second bottom
            // for(int p = start; p < nrend; p++)
            // {
            //     for (d = 0; d < DIM; d++)
            //     {
            //         sum_uro_without_p[d] = 0;
            //         for(int k = start; k < nrend; k++)
            //         {
            //             if (k == p) continue;
            //             sum_uro_without_p[d] += v[k][d]*ro_md;
            //             // if (k%1000 == 0) printf("%d\n", k);
            //
            //         }
            //         // printf("w");
            //         pr_v[d] = u_tilde[d]*ro_tilde - sum_uro_without_p[d];
            //
            //         second_bottom[d] = -mv[d]*pr_v[d]/ro_md;
            //     }
            //     second_bottom_v += SUM(second_bottom);
            //
            // }

            if(n % 10000 == 0) printf("%d ", n);
            // if(n > 8000 == 0) printf("/%d/", n);

            // final

            beta[n-start] = (arr[ind].first_top + arr[ind].second_top)/(first_bottom_v + arr[ind].second_bottom);

        }
    }

    for(int i = 0; i < fh->Ntot; i++)
    {
        printf("\n%e\t", arr[i].first_top);
        printf("%e\t", arr[i].second_top);
        printf("%e\t", arr[i].second_bottom);
        printf("%e\t", beta[i + 10000]);
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
                        // vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + S*(1 - S)*beta_term[d])*invro_dt; // *
                        // vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + beta[n-start]*S*(1 - S)*beta_term[d])*invro_dt; // *
                        vn           = lg*v[n][d] + (1 - S)*f[n][d]*w_dt + (S*f_fh[d] + alpha_term[d] + fh->beta*S*(1 - S)*beta_term[d])*invro_dt; // *
                        // beta_term we multiply by beta coef for each particle n (beta_term * beta_p)  beta[n]*...*ber
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
