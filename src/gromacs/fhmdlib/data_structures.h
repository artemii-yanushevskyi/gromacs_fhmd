#ifndef FHMD_DATA_STRUCTURES_H_
#define FHMD_DATA_STRUCTURES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "params.h"
#include "gromacs/topology/topology.h"      /* GROMACS gmx_mtop_t  structure definition */
#include "gromacs/domdec/domdec_struct.h"   /* GROMACS gmx_domdec_t structure definition */
#include "gromacs/mdtypes/commrec.h"        /* GROMACS MPI definitions, 't_commrec', MASTER(), PAR(), etc. */
#include "gromacs/gmxlib/network.h"         /* GROMACS MPI functions, gmx_bcast(), gmx_sumf(), etc. */
#include "gromacs/math/vectypes.h"          /* GROMACS vector types: rvec, dvec, ivec, etc. */
#include "gromacs/math/vec.h"               /* GROMACS vector operations: copy_ivec, dvec_add, etc. */


typedef struct FH_arrays                    /* FH/MD arrays */
{
    double      ro_md, ro_fh;               /* densities */
    double      inv_ro;                     /* inverse density: 1/ro_md */
    dvec        u_md, u_fh;                 /* velocities */
    dvec        uro_md;                     /* momentum */
    dvec        f_fh;                       /* FH force */
    dvec        alpha_term;                 /* alpha term for du/dt equation */
    dvec        beta_term;                  /* beta term for du/dt equation */

    double      delta_ro;                   /* delta rho for 1-way coupling or ro_prime for 2-way */
    dvec        grad_ro;                    /* grad of density */
    matrix      alpha_u_grad;               /* preliminary alpha-term [u-index][grad-index] */

    double      ro_prime, ron_prime;        /* density prime */
    double      ro_star, ron_star;          /* density star */
    dvec        m_prime, mn_prime;          /* m prime */
    dvec        m_star, mn_star;            /* m star */

    double      ros_md, ropr_md;            /* sources of MD density */
    dvec        uros_md, uropr_md;          /* sources of MD momentum */

    double      S;                          /* S parameter in the FH cell centres */
    dvec        Sf;                         /* S parameter in the FH cell faces */

    double      p, pn;                      /* FH pressure */
    dvec        rof, rofn, pf, pfn;         /* FH flux variables */
    matrix      uf, ufn;                    /* FH flux velocities */
    double      ro_fh_n;                    /* FH density (new time layer) */
    dvec        u_fh_n;                     /* FH velocity (new time layer) */
    matrix      rans;                       /* FH random stress */
} FH_arrays;


typedef struct FH_grid          /* Computational grid */
{
    dvec       *c;              /* FH cell centres coordinates */
    dvec       *n;              /* FH cell nodes coordinates */
    dvec       *h;              /* FH cell steps */
    double     *vol;            /* FH cell volume */
    double     *ivol;           /* 1/cellVolume */
} FH_grid;


typedef struct MD_stat          /* Particle statistics */
{
    int         N, n;
    double      invN;

    double     *avg_rho_md_cell, *avg_rho_fh_cell;

    double      davg_rho_md,  davg_rho_fh;
    double      davg_rho2_md, davg_rho2_fh;
    dvec        davg_u_md,    davg_u_fh;
    dvec        davg_u2_md,   davg_u2_fh;

    double      avg_rho_md,   avg_rho_fh;
    double      std_rho_md,   std_rho_fh;
    dvec        avg_u_md,     avg_u_fh;
    dvec        std_u_md,     std_u_fh;
} MD_stat;


typedef struct FHMD
{
    FH_arrays  *arr;            /* FH/MD arrays */
    FH_grid     grid;           /* FH grid */
    MD_stat     stat;           /* Particle statistics */
    int        *ind;            /* FH cell number for each atom */
    ivec       *indv;           /* 3-component FH cell number for each atom (vector) */
    double     *mpi_linear;     /* Linear array to summarise MDFH arrays */

    int         scheme;         /* 0 - Pure MD, 1 - One-way coupling, 2 - Two-way coupling */
    double      S;              /* Parameter S (-1 - fixed sphere, -2 - moving sphere) */
    double      R1;             /* MD sphere radius for variable S, [0..1] */
    double      R2;             /* FH sphere radius for variable S, [0..1] */
    double      Smin;           /* Minimum S for variable S */
    double      Smax;           /* Maximum S for variable S */
    double      R12, R22, RS;   /* Derived variables from R1, R2, Smin, Smax */
    double      alpha;          /* Alpha parameter for dx/dt and du/dt equations, nm^2/ps */
    double      beta;           /* Beta parameter, nm^2/ps or ps^-1 depending on the scheme */
    double      gamma_x;        /* Gamma_x parameter (density fluctuations dissipator), ps^-1 */
    double      gamma_u;        /* Gamma_u parameter (velocity fluctuations dissipator), ps^-1 */
    double      std_rho;        /* Analytical STD of density */
    double      std_u;          /* Analytical STD of velocity */

    ivec        N;              /* Number of FH cells along each direction */
    dvec        box;            /* Box size */
    dvec        box05;          /* Half of box size */
    dvec        protein_com;    /* Protein COM coordinates */
    double      protein_mass;   /* Protein mass */
    int         protein_N;      /* Number of atoms in the protein */
    double      box_volume;     /* Volume of the box, nm^3 */
    double      total_density;  /* Total density of the box, a.m.u./nm^3 */
    int         Ntot;           /* Total number of FH cells */
    int         step_MD;        /* Current MD time step */
    int         Noutput;        /* Write arrays to files every Noutput MD time steps (0 - do not write) */

    int         FH_EOS;         /* EOS: 0 - Liquid Argon, 1 - SPC/E water */
    FHMD_EOS    eos;            /* Equation of state */
    int         FH_step;        /* dt_FH = FH_step * dt_MD */
    int         FH_equil;       /* Number of time steps for the FH model equilibration */
    double      FH_dens;        /* FH mean density */
    double      FH_temp;        /* FH mean temperature */
    double      FH_blend;       /* FH Blending: -1 - dynamic, or define static blending parameter (0.0 = Central Diff., 1.0 = Classic CABARET) */
    double      dt_FH;          /* FH time step */
} FHMD;

#endif /* FHMD_DATA_STRUCTURES_H_ */
