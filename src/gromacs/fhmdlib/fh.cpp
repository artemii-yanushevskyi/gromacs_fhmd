#include "data_structures.h"
#include "macro.h"
#include "fh.h"


void FH_init(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    switch(fh->eos)
    {
    case eos_argon:
        MU     = 0;
        KAPPA  = 0;
        EOS_A  = 0;
        EOS_B  = 0;
        EOS_C  = 0;
        EOS_D  = 0;
        EOS_E  = 0;
        P_INIT = 0;
        SOUND  = 0;
        break;
    case eos_spce:
        MU     = FHMD_EOS_SPCE_MU;
        KAPPA  = FHMD_EOS_SPCE_KAPPA;
        EOS_A  = FHMD_EOS_SPCE_A;
        EOS_B  = FHMD_EOS_SPCE_B;
        EOS_C  = FHMD_EOS_SPCE_C;
        P_INIT = fh->FH_dens*(fh->FH_dens*EOS_A + EOS_B) + EOS_C;
        SOUND  = sqrt(2.0*fh->FH_dens*EOS_A + EOS_B);
        break;
    }

#ifdef FHMD_DEBUG
    printf(MAKE_YELLOW "FHMD DEBUG: MU = %g, KAPPA = %g\n", MU, KAPPA);
    printf(MAKE_YELLOW "FHMD DEBUG: A = %g, B = %g, C = %g\n", EOS_A, EOS_B, EOS_C);
    printf(MAKE_YELLOW "FHMD DEBUG: Initial Pressure = %g, Speed of Sound = %g\n", P_INIT, SOUND*1000.0);
    printf(RESET_COLOR "\n");
#endif

    // Viscosity terms for viscous stress
    VISC1 = 4.0/3.0 + KAPPA/MU;
    VISC2 = KAPPA/MU - 2.0/3.0;

    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_fh   = fh->FH_dens;
        arr[i].ro_fh_n = fh->FH_dens;
        arr[i].p       = P_INIT;
        arr[i].pn      = P_INIT;

        for(int d = 0; d < DIM; d++)
        {
            arr[i].rof[d]  = fh->FH_dens;
            arr[i].rofn[d] = fh->FH_dens;
            arr[i].pf[d]   = P_INIT;
            arr[i].pfn[d]  = P_INIT;
        }
    }

    T        = 0;       // Current time
    STEP     = 0;       // Current time step

    // Reset statistics
    T_AVG    = 0;
    NT_AVG   = 0;
    RHO_AVG  = 0;
    RHO2_AVG = 0;
    P_AVG    = 0;
    N_AVG    = 0;

    for(int d = 0; d < DIM; d++)
    {
        U_AVG[d]  = 0;
        U2_AVG[d] = 0;
    }
}


void FH_predictor(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec   ind;
    double DT = fh->dt_FH;
    dvec   F, PG, TAU, TAURAN;
    matrix TAUL, TAUR;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    // Mass flux
                    F[d] = (arr[R].rof[d]*arr[R].uf[d][d] - arr[L].rof[d]*arr[L].uf[d][d])/fh->grid.h[C][d];
                }

                // Mass conservation
                arr[C].ro_fh_n = arr[C].ro_fh - 0.5*DT*SUM(F);

                for(int d = 0; d < DIM; d++)
                {
                    // Pressure gradient
                    PG[d] = (arr[R].pf[d] - arr[L].pf[d])/fh->grid.h[C][d];
                }

                for(int dim = 0; dim < DIM; dim++)
                {
                    VISCOUS_FLUX(u_fh);     // UNIFORM GRID!

                    // Momentum flux, viscous and random stress
                    for(int d = 0; d < DIM; d++)
                    {
                        F[d]      = (arr[R].rof[d]*arr[R].uf[d][d]*arr[R].uf[dim][d] - arr[L].rof[d]*arr[L].uf[d][d]*arr[L].uf[dim][d])/fh->grid.h[C][d];
                        TAU[d]    = (TAUR[dim][d] - TAUL[dim][d])/fh->grid.h[C][d];
                        TAURAN[d] = (arr[R].rans[dim][d] - arr[L].rans[dim][d])/fh->grid.h[C][d];
                    }

                    // FH force
                    arr[C].f_fh[dim] = -SUM(F) - PG[dim] + SUM(TAU) + SUM(TAURAN);

                    // Momentum conservation
                    arr[C].u_fh_n[dim] = (arr[C].ro_fh*arr[C].u_fh[dim] + 0.5*DT*arr[C].f_fh[dim])/arr[C].ro_fh_n;
                }

                // Pressure
                switch(fh->eos)
                {
                case eos_argon:
                    arr[C].pn = EOS_D*(1./3.*EOS_A*EOS_A*pow(EOS_E*arr[C].ro_fh_n - EOS_B, 3) + EOS_C);
                    break;
                case eos_spce:
                    arr[C].pn = arr[C].ro_fh_n*(arr[C].ro_fh_n*EOS_A + EOS_B) + EOS_C;
                    break;
                }
            }
        }
    }
}


void FH_corrector(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec   ind;
    double DT = fh->dt_FH;
    dvec   F, PG, TAU, TAURAN;
    matrix TAUL, TAUR;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    // Mass flux
                    F[d] = (arr[R].rofn[d]*arr[R].ufn[d][d] - arr[L].rofn[d]*arr[L].ufn[d][d])/fh->grid.h[C][d];
                }

                // Mass conservation
                arr[C].ro_fh = arr[C].ro_fh_n - 0.5*DT*SUM(F);

                for(int d = 0; d < DIM; d++)
                {
                    // Pressure gradient
                    PG[d] = (arr[R].pfn[d] - arr[L].pfn[d])/fh->grid.h[C][d];
                }

                for(int dim = 0; dim < DIM; dim++)
                {
                    VISCOUS_FLUX(u_fh_n);     // UNIFORM GRID!

                    // Momentum flux, viscous and random stress
                    for(int d = 0; d < DIM; d++)
                    {
                        F[d]      = (arr[R].rofn[d]*arr[R].ufn[d][d]*arr[R].ufn[dim][d] - arr[L].rofn[d]*arr[L].ufn[d][d]*arr[L].ufn[dim][d])/fh->grid.h[C][d];
                        TAU[d]    = (TAUR[dim][d] - TAUL[dim][d])/fh->grid.h[C][d];
                        TAURAN[d] = (arr[R].rans[dim][d] - arr[L].rans[dim][d])/fh->grid.h[C][d];
                    }

                    // FH force
                    arr[C].f_fh[dim] = -SUM(F) - PG[dim] + SUM(TAU) + SUM(TAURAN);

                    // Momentum conservation
                    arr[C].u_fh[dim] = (arr[C].ro_fh_n*arr[C].u_fh_n[dim] + 0.5*DT*arr[C].f_fh[dim])/arr[C].ro_fh;
                }

                // Pressure
                switch(fh->eos)
                {
                case eos_argon:
                    arr[C].p = EOS_D*(1./3.*EOS_A*EOS_A*pow(EOS_E*arr[C].ro_fh - EOS_B, 3) + EOS_C);
                    break;
                case eos_spce:
                    arr[C].p = arr[C].ro_fh*(arr[C].ro_fh*EOS_A + EOS_B) + EOS_C;
                    break;
                }
            }
        }
    }
}


void FH_char(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    const double RO1 = 1.0/fh->FH_dens;
    const double DT  = fh->dt_FH;

    ivec ind;
    int  dn;

    double RCLN, RCL, QCLN, QCL, VCLN, VCL, VL, VPL, WCLN, WCL, WL, WPL, RL, RPL, QL, QPL;

    /*
    double L1L,L2L,L3L,L1R,L2R,L3R;
    double RCL,RCLN,RCR,RCRN,RL,RPL,RR,RPR,RPNL,RPNR,RPN;
    double QCL,QCLN,QCR,QCRN,QL,QPL,QR,QPR,QPNL,QPNR,QPN;
    double VCL,VCLN,VCR,VCRN,VL,VPL,VR,VPR,VPNL,VPNR,VPN;
    double WCL,WCLN,WCR,WCRN,WL,WPL,WR,WPR,WPNL,WPNR,WPN;
    double QQ1,QQ2,QQ3,QQ4,RMIN,RMAX,QMIN,QMAX,VMIN,VMAX,WMIN,WMAX;
    double RO1,P1,P2,P3;
    double RPN1,QPN1,VPN1,WPN1,RPN2,QPN2,VPN2,WPN2;
*/

    for(int d = 0; d < DIM; d++)
    {
        for(int k = 0; k < NZ; k++)
        {
            for(int j = 0; j < NY; j++)
            {
                for(int i = 0; i < NX; i++)
                {
                    ASSIGN_IND(ind, i, j, k);

                    RCLN = arr[L0].u_fh_n[d] + SOUND*log(arr[L0].ro_fh_n*RO1);
                    RCL  = arr[L0].u_fh[d]   + SOUND*log(arr[L0].ro_fh*RO1);

                    QCLN = arr[L0].u_fh_n[d] - SOUND*log(arr[L0].ro_fh_n*RO1);
                    QCL  = arr[L0].u_fh[d]   - SOUND*log(arr[L0].ro_fh*RO1);

                    dn = d + 1; if(dn >= DIM) dn = 0;

                    VCLN = arr[L0].u_fh_n[dn];
                    VCL  = arr[L0].u_fh[dn];
                    VL   = arr[L0].uf[dn][d];
                    VPL  = arr[L1].uf[dn][d];

                    dn++; if(dn >= DIM) dn = 0;

                    WCLN = arr[L0].u_fh_n[dn];
                    WCL  = arr[L0].u_fh[dn];
                    WL   = arr[L0].uf[dn][d];
                    WPL  = arr[L1].uf[dn][d];

                    RL   = arr[L0].uf[d][d] + SOUND*log(arr[L0].rof[d]*RO1);
                    RPL  = arr[L1].uf[d][d] + SOUND*log(arr[L1].rof[d]*RO1);
                    QL   = arr[L0].uf[d][d] - SOUND*log(arr[L0].rof[d]*RO1);
                    QPL  = arr[L1].uf[d][d] - SOUND*log(arr[L1].rof[d]*RO1);




                }
            }
        }
    }

    for(k=1; k<=N3-1; k++)
    {
        for(j=1; j<=N2-1; j++)
        {


            for(i=1; i<=N1-1; i++)
            {
            // --------------- LOCAL RIEMANN INVARIANTS ----------------------------------------
            // =================================================================================
                P1 = SOUND*log(rocn[i  ][j][k]*RO1);
                P2 = SOUND*log(roc [i  ][j][k]*RO1);
                P3 = SOUND*log(rox [i+1][j][k]*RO1);
            // =================================================================================
                RCRN = uxcn[i  ][j][k]+P1;      QCRN = uxcn[i  ][j][k]-P1;      VCRN = uycn[i  ][j][k];     WCRN = uzcn[i  ][j][k];
                RCR  = uxc [i  ][j][k]+P2;      QCR  = uxc [i  ][j][k]-P2;      VCR  = uyc [i  ][j][k];     WCR  = uzc [i  ][j][k];
            // =================================================================================
                RR   = uxx [i+1][j][k]+P3;      QR   = uxx [i+1][j][k]-P3;      VR   = uyx [i+1][j][k];     WR   = uzx [i+1][j][k];
            // =================================================================================
                RPR  = RPL; //uxx [i  ][j  ]+SOUND*log(rox [i  ][j  ]*RO1)
                QPR  = QPL; //uxx [i  ][j  ]-SOUND*log(rox [i  ][j  ]*RO1)
                VPR  = VPL;
                WPR  = WPL;
            // =================================================================================

            // --------------- CHARACTERISTIC SPEEDS -------------------------------------------
                L1L = uxcn[i-1][j][k]+SOUND;  L2L = uxcn[i-1][j][k]-SOUND;   L3L = uxcn[i-1][j][k];
                L1R = uxcn[i  ][j][k]+SOUND;  L2R = uxcn[i  ][j][k]-SOUND;   L3R = uxcn[i  ][j][k];
            // ---------------------------------------------------------------------------------

            // --------------- LEFT CELL -------------------------------------------------------
            // ----- EXTRAPOLATION FOR THE LEFT GOING WAVES ------------------------------------
                RPNL = 2.0*RCLN-RPL;      QPNL = 2.0*QCLN-QPL;      VPNL = 2.0*VCLN-VPL;    WPNL = 2.0*WCLN-WPL;
            // --------------- RIGHT CELL ----------------------------------------------------
            // ----- EXTRAPOLATION FOR THE RIGHT GOING WAVES -----------------------------------
                RPNR = 2.0*RCRN-RPR;      QPNR = 2.0*QCRN-QPR;      VPNR = 2.0*VCRN-VPR;    WPNR = 2.0*WCRN-WPR;
            // ---------------------------------------------------------------------------------

//      QQ1 = 2.0*(RCLN-RCL) + DT*L1L*(RPL-RL)*HXI;
//      QQ2 = 2.0*(QCLN-QCL) + DT*L2L*(QPL-QL)*HXI;
//      QQ3 = 2.0*(VCLN-VCL) + DT*L3L*(VPL-VL)*HXI;
//      QQ4 = 2.0*(WCLN-WCL) + DT*L3L*(WPL-WL)*HXI;

        QQ1 = 2.0*(RCLN-RCL) + DT*L1L*(RPL-RL)/hc[i];
        QQ2 = 2.0*(QCLN-QCL) + DT*L2L*(QPL-QL)/hc[i];
        QQ3 = 2.0*(VCLN-VCL) + DT*L3L*(VPL-VL)/hc[i];
        QQ4 = 2.0*(WCLN-WCL) + DT*L3L*(WPL-WL)/hc[i];

        RMAX = DMAX1(RL,RCLN,RPL)+QQ1;        RMIN = DMIN1(RL,RCLN,RPL)+QQ1;
        QMAX = DMAX1(QL,QCLN,QPL)+QQ2;        QMIN = DMIN1(QL,QCLN,QPL)+QQ2;
        VMAX = DMAX1(VL,VCLN,VPL)+QQ3;        VMIN = DMIN1(VL,VCLN,VPL)+QQ3;
        WMAX = DMAX1(WL,WCLN,WPL)+QQ4;        WMIN = DMIN1(WL,WCLN,WPL)+QQ4;

        if(RPNL>RMAX) RPNL = RMAX;        if(RPNL<RMIN) RPNL = RMIN;
        if(QPNL>QMAX) QPNL = QMAX;        if(QPNL<QMIN) QPNL = QMIN;
        if(VPNL>VMAX) VPNL = VMAX;        if(VPNL<VMIN) VPNL = VMIN;
        if(WPNL>WMAX) WPNL = WMAX;        if(WPNL<WMIN) WPNL = WMIN;

//      QQ1 = 2.0*(RCRN-RCR) + DT*L1R*(RR-RPR)*HXI;
//      QQ2 = 2.0*(QCRN-QCR) + DT*L2R*(QR-QPR)*HXI;
//      QQ3 = 2.0*(VCRN-VCR) + DT*L3R*(VR-VPR)*HXI;
//      QQ4 = 2.0*(WCRN-WCR) + DT*L3R*(WR-WPR)*HXI;

        QQ1 = 2.0*(RCRN-RCR) + DT*L1R*(RR-RPR)/hc[i];
        QQ2 = 2.0*(QCRN-QCR) + DT*L2R*(QR-QPR)/hc[i];
        QQ3 = 2.0*(VCRN-VCR) + DT*L3R*(VR-VPR)/hc[i];
        QQ4 = 2.0*(WCRN-WCR) + DT*L3R*(WR-WPR)/hc[i];

        RMAX = DMAX1(RR,RCRN,RPR)+QQ1;   RMIN = DMIN1(RR,RCRN,RPR)+QQ1;
        QMAX = DMAX1(QR,QCRN,QPR)+QQ2;   QMIN = DMIN1(QR,QCRN,QPR)+QQ2;
        VMAX = DMAX1(VR,VCRN,VPR)+QQ3;   VMIN = DMIN1(VR,VCRN,VPR)+QQ3;
        WMAX = DMAX1(WR,WCRN,WPR)+QQ4;   WMIN = DMIN1(WR,WCRN,WPR)+QQ4;

        if(RPNR>RMAX) RPNR = RMAX;    if(RPNR<RMIN) RPNR = RMIN;
        if(QPNR>QMAX) QPNR = QMAX;    if(QPNR<QMIN) QPNR = QMIN;
        if(VPNR>VMAX) VPNR = VMAX;    if(VPNR<VMIN) VPNR = VMIN;
        if(WPNR>WMAX) WPNR = WMAX;    if(WPNR<WMIN) WPNR = WMIN;

        if((L1L+L1R)>=0) { RPN1 = RPNL; } else { RPN1 = RPNR; }
        if((L2L+L2R)>=0) { QPN1 = QPNL; } else { QPN1 = QPNR; }
        if((L3L+L3R)>=0) { VPN1 = VPNL;   WPN1 = WPNL; } else { VPN1 = VPNR;   WPN1 = WPNR; }

    // ---------------------------------------------------------------------------------
        RPN2 = 0.5*(RPNL+RPNR);    QPN2 = 0.5*(QPNL+QPNR);  VPN2 = 0.5*(VPNL+VPNR); WPN2 = 0.5*(WPNL+WPNR);
    // ---------------------------------------------------------------------------------
//      QQ1 = 0.5*(2.0*(RCRN-RCR+RCLN-RCL) +DT*0.5*(L1L+L1R)*(RR-RPR+RPL-RL)*HXI);
//      QQ2 = 0.5*(2.0*(QCRN-QCR+QCLN-QCL) +DT*0.5*(L2L+L2R)*(QR-QPR+QPL-QL)*HXI);
//      QQ3 = 0.5*(2.0*(VCRN-VCR+VCLN-VCL) +DT*0.5*(L3L+L3R)*(VR-VPR+VPL-VL)*HXI);
//      QQ4 = 0.5*(2.0*(WCRN-WCR+WCLN-WCL) +DT*0.5*(L3L+L3R)*(WR-WPR+WPL-WL)*HXI);

        QQ1 = 0.5*(2.0*(RCRN-RCR+RCLN-RCL) +DT*0.5*(L1L+L1R)*(RR-RPR+RPL-RL)/hc[i]);
        QQ2 = 0.5*(2.0*(QCRN-QCR+QCLN-QCL) +DT*0.5*(L2L+L2R)*(QR-QPR+QPL-QL)/hc[i]);
        QQ3 = 0.5*(2.0*(VCRN-VCR+VCLN-VCL) +DT*0.5*(L3L+L3R)*(VR-VPR+VPL-VL)/hc[i]);
        QQ4 = 0.5*(2.0*(WCRN-WCR+WCLN-WCL) +DT*0.5*(L3L+L3R)*(WR-WPR+WPL-WL)/hc[i]);

        RMAX = DMAX1(0.5*(RL+RR),0.5*(RCLN+RCRN),0.5*(RPL+RPR))+QQ1;
        RMIN = DMIN1(0.5*(RL+RR),0.5*(RCLN+RCRN),0.5*(RPL+RPR))+QQ1;
        QMAX = DMAX1(0.5*(QL+QR),0.5*(QCLN+QCRN),0.5*(QPL+QPR))+QQ2;
        QMIN = DMIN1(0.5*(QL+QR),0.5*(QCLN+QCRN),0.5*(QPL+QPR))+QQ2;
        VMAX = DMAX1(0.5*(VL+VR),0.5*(VCLN+VCRN),0.5*(VPL+VPR))+QQ3;
        VMIN = DMIN1(0.5*(VL+VR),0.5*(VCLN+VCRN),0.5*(VPL+VPR))+QQ3;
        WMAX = DMAX1(0.5*(WL+WR),0.5*(WCLN+WCRN),0.5*(WPL+WPR))+QQ4;
        WMIN = DMIN1(0.5*(WL+WR),0.5*(WCLN+WCRN),0.5*(WPL+WPR))+QQ4;

        if(RPN2>RMAX) {RPN2 = RMAX;}    if(RPN2<RMIN) {RPN2 = RMIN;}
        if(QPN2>QMAX) {QPN2 = QMAX;}    if(QPN2<QMIN) {QPN2 = QMIN;}
        if(VPN2>VMAX) {VPN2 = VMAX;}    if(VPN2<VMIN) {VPN2 = VMIN;}
        if(WPN2>WMAX) {WPN2 = WMAX;}    if(WPN2<WMIN) {WPN2 = WMIN;}

    // ------------- blend -------------------------------------------------------------
        RPN = Fblend*RPN1+(1.0-Fblend)*RPN2;
        QPN = Fblend*QPN1+(1.0-Fblend)*QPN2;
        VPN = Fblend*VPN1+(1.0-Fblend)*VPN2;
        WPN = Fblend*WPN1+(1.0-Fblend)*WPN2;
    // ------------- blend -------------------------------------------------------------

            // ================================================================================]
            // ------------------ START THE APProxIMATE RIEMANN SOLVER -------------------------
            // ========= fluxes ================================================================
                roxn[i  ][j][k] = RO_INIT*exp(0.5*(RPN-QPN)/SOUND);
                uxxn[i  ][j][k] = 0.5*(RPN+QPN);
                uyxn[i  ][j][k] = VPN;
                uzxn[i  ][j][k] = WPN;

                if(fh_EOS == hmd_EOS_ARGON)
                    pxn[i  ][j][k] = D*(1./3.*A*A*pow(E*roxn[i  ][j][k]-B,3)+C);
                else if(fh_EOS == hmd_EOS_SPC)
                    pxn[i  ][j][k] = roxn[i  ][j][k]*(roxn[i  ][j][k]*A+B)+C;
                else {
                    printf("\nERROR in FHLJModel3D - unknown EOS\n\n");
                    exit(3);
                }

            // ================================================================================]

            // ---------------- For next time step -----------------------------------------------
                RCLN= RCRN;      QCLN= QCRN;        VCLN= VCRN;     WCLN= WCRN;
                RCL = RCR;       QCL = QCR;         VCL = VCR;      WCL = WCR;
                 RL = RPL;        QL = QPL;          VL = VPL;       WL = WPL;
                RPL = RR;        QPL = QR;          VPL = VR;       WPL = WR;
            // ===================================================================================
            }
        }
    }

}


void compute_random_stress(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    double FLUCT, TCELL, TRACE;
    double ranv[7];
    double Tcurr, Tf;

    T_INST = 0;

    for(int i = 0; i < fh->Ntot; i++)
    {
        // Fluctuation-Dissipation factor
        FLUCT = sqrt(2.0*FHMD_kB*fh->FH_temp)/sqrt(fh->dt_FH*fh->grid.vol[i]);

        TCELL   = arr[i].ro_fh*(arr[i].u_fh[0]*arr[i].u_fh[0] + arr[i].u_fh[1]*arr[i].u_fh[1] + arr[i].u_fh[2]*arr[i].u_fh[2])*fh->grid.vol[i]/(3.0*FHMD_kB);
        T_INST += TCELL;

        for(int k = 0; k < 7; k++)
            ranv[k] = DRNOR();

        TRACE = (ranv[1] + ranv[2] + ranv[3])/3.0;

        arr[i].rans[0][0] = FLUCT*(sqrt(2.0*MU)*(ranv[1] - TRACE) + sqrt(KAPPA)*ranv[0]);
        arr[i].rans[1][1] = FLUCT*(sqrt(2.0*MU)*(ranv[2] - TRACE) + sqrt(KAPPA)*ranv[0]);
        arr[i].rans[2][2] = FLUCT*(sqrt(2.0*MU)*(ranv[3] - TRACE) + sqrt(KAPPA)*ranv[0]);

        arr[i].rans[0][1] = FLUCT*(sqrt(MU)*ranv[4]);
        arr[i].rans[0][2] = FLUCT*(sqrt(MU)*ranv[5]);
        arr[i].rans[1][2] = FLUCT*(sqrt(MU)*ranv[6]);

        arr[i].rans[1][0] = arr[i].rans[0][1];
        arr[i].rans[2][0] = arr[i].rans[0][2];
        arr[i].rans[2][1] = arr[i].rans[1][2];
    }

    // Compute current temperature and compute dynamic blending factor if enabled
    if(FHMD_BLENDING < 0)
    {
        blend = 0.0;
        Tcurr = T_INST/(double)(fh->Ntot);

        if(Tcurr > fh->FH_temp)
        {
            Tf = Tcurr/fh->FH_temp;
            blend = 1.0*(Tf - 1.0);
            if(blend > 1.0) blend = 1.0;
        }
    }
}


void swap_var(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    for(int i = 0; i < fh->Ntot; i++)
    {
        for(int d = 0; d < DIM; d++)
        {
            arr[i].rof[d] = arr[i].rofn[d];
            arr[i].pf[d]  = arr[i].pfn[d];

            for(int d1 = 0; d1 < DIM; d1++)
            {
                arr[i].uf[d][d1] = arr[i].ufn[d][d1];
            }
        }
    }
}


void collect_statistics(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    T_AVG += T_INST/(double)(fh->Ntot);
    NT_AVG++;

    for(int i = 0; i < fh->Ntot; i++)
    {
        RHO_AVG  += arr[i].ro_fh;
        RHO2_AVG += (arr[i].ro_fh - fh->FH_dens)*(arr[i].ro_fh - fh->FH_dens);
        P_AVG    += arr[i].p;
        N_AVG++;

        for(int d = 0; d < DIM; d++)
        {
            U_AVG[d]  += arr[i].u_fh[d];
            U2_AVG[d] += arr[i].u_fh[d]*arr[i].u_fh[d];
        }
    }
}


void update_statistics(FHMD *fh)
{
    const double DN_AVG   = (double)(N_AVG);

    avg_rho = RHO_AVG/DN_AVG;
    avg_p   = P_AVG/DN_AVG;
    std_rho = sqrt(fabs(RHO2_AVG/DN_AVG));

    for(int d = 0; d < DIM; d++)
    {
        avg_u[d]     = U_AVG[d]/DN_AVG;
        std_u[d]     = sqrt(fabs(U2_AVG[d] - avg_u[d]*avg_u[d]));
        if(std_rho > 0)
            avg_sound[d] = fh->FH_dens*std_u[d]/std_rho;
        else
            avg_sound[d] = 0;
    }

    sound = (avg_sound[0] + avg_sound[1] + avg_sound[2])/3.;
    avg_T = T_AVG/(double)(NT_AVG);
}


void FH_do_single_timestep(FHMD *fh)
{
    compute_random_stress(fh);

    FH_predictor(fh);
    FH_char(fh);
    FH_corrector(fh);

    swap_var(fh);
}


void FH_equilibrate(FHMD *fh)
{
    const double STD_RHO  = sqrt(fh->FH_temp*fh->FH_dens*FHMD_kB/(fh->box_volume/(double)(fh->Ntot)))/SOUND;
    const double STD_U    = sqrt(fh->FH_temp*FHMD_kB/(fh->FH_dens*fh->box_volume/(double)(fh->Ntot)));
    const int    N_OUTPUT = 1000;

    printf(MAKE_BLUE "FHMD: FH equilibration in process...\n\n");
    printf("%8s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %6s\n",
           "Step", "STD rho", "STD Ux", "STD Uy", "STD Uz", "C_s, m/s", "T, K", "<T>, K", "<rho>", "<P>", "<Ux>", "<Uy>", "<Uz>", "blend");
    printf("---------------------------------------------------------------------------------------------------------------------------------------\n");
    printf(MAKE_LIGHT_BLUE "->       %9.4f %9.5f %9.5f %9.5f %9.3f %9.4f %9.4f %9.4f %9.4f      <- Theoretical Values",
           STD_RHO, STD_U, STD_U, STD_U, SOUND*1000.0, fh->FH_temp, fh->FH_temp, fh->FH_dens, P_INIT);
    printf(MAKE_BLUE "\n");

    while(STEP <= fh->FH_equil)
    {
        collect_statistics(fh);

        if(!(STEP % N_OUTPUT))
        {
            update_statistics(fh);

            printf("\r%8d %9.4f %9.5f %9.5f %9.5f %9.3f %9.4f %9.4f %9.4f %9.4f %9.5f %9.5f %9.5f %6.4f",
                   STEP, std_rho, std_u[0], std_u[1], std_u[2], sound*1000.0, T_INST/(double)(fh->Ntot), avg_T,
                   avg_rho, avg_p, avg_u[0], avg_u[1], avg_u[2], blend);

            fflush(stdout);
        }

        FH_do_single_timestep(fh);

        T += fh->dt_FH;
        STEP++;
    }

    printf("\n---------------------------------------------------------------------------------------------------------------------------------------\n");
    printf(RESET_COLOR "\n");
    fflush(stdout);
}


/*
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void CHARX()
{
// ---------------------------------------------------------------------------------------
    int i,j,k;
    double L1L,L2L,L3L,L1R,L2R,L3R;
    double RCL,RCLN,RCR,RCRN,RL,RPL,RR,RPR,RPNL,RPNR,RPN;
    double QCL,QCLN,QCR,QCRN,QL,QPL,QR,QPR,QPNL,QPNR,QPN;
    double VCL,VCLN,VCR,VCRN,VL,VPL,VR,VPR,VPNL,VPNR,VPN;
    double WCL,WCLN,WCR,WCRN,WL,WPL,WR,WPR,WPNL,WPNR,WPN;
    double QQ1,QQ2,QQ3,QQ4,RMIN,RMAX,QMIN,QMAX,VMIN,VMAX,WMIN,WMAX;
    double RO1,P1,P2,P3;
    double RPN1,QPN1,VPN1,WPN1,RPN2,QPN2,VPN2,WPN2;

    RO1 = 1.0/RO_INIT;

    for(k=1; k<=N3-1; k++)
    {
        for(j=1; j<=N2-1; j++)
        {
            RCLN = uxcn[0][j][k]+SOUND*log(rocn[0][j][k]*RO1);        RCL  = uxc [0][j][k]+SOUND*log(roc [0][j][k]*RO1);
            QCLN = uxcn[0][j][k]-SOUND*log(rocn[0][j][k]*RO1);        QCL  = uxc [0][j][k]-SOUND*log(roc [0][j][k]*RO1);
            VCLN = uycn[0][j][k];                                  VCL  = uyc [0][j][k];
            WCLN = uzcn[0][j][k];                                  WCL  = uzc [0][j][k];

            RL   = uxx [0][j][k]+SOUND*log(rox [0][j][k]*RO1);        RPL  = uxx [1][j][k]+SOUND*log(rox [1][j][k]*RO1);
            QL   = uxx [0][j][k]-SOUND*log(rox [0][j][k]*RO1);        QPL  = uxx [1][j][k]-SOUND*log(rox [1][j][k]*RO1);
            VL   = uyx [0][j][k];                                  VPL  = uyx [1][j][k];
            WL   = uzx [0][j][k];                                  WPL  = uzx [1][j][k];

            for(i=1; i<=N1-1; i++)
            {
            // --------------- LOCAL RIEMANN INVARIANTS ----------------------------------------
            // =================================================================================
                P1 = SOUND*log(rocn[i  ][j][k]*RO1);
                P2 = SOUND*log(roc [i  ][j][k]*RO1);
                P3 = SOUND*log(rox [i+1][j][k]*RO1);
            // =================================================================================
                RCRN = uxcn[i  ][j][k]+P1;      QCRN = uxcn[i  ][j][k]-P1;      VCRN = uycn[i  ][j][k];     WCRN = uzcn[i  ][j][k];
                RCR  = uxc [i  ][j][k]+P2;      QCR  = uxc [i  ][j][k]-P2;      VCR  = uyc [i  ][j][k];     WCR  = uzc [i  ][j][k];
            // =================================================================================
                RR   = uxx [i+1][j][k]+P3;      QR   = uxx [i+1][j][k]-P3;      VR   = uyx [i+1][j][k];     WR   = uzx [i+1][j][k];
            // =================================================================================
                RPR  = RPL; //uxx [i  ][j  ]+SOUND*log(rox [i  ][j  ]*RO1)
                QPR  = QPL; //uxx [i  ][j  ]-SOUND*log(rox [i  ][j  ]*RO1)
                VPR  = VPL;
                WPR  = WPL;
            // =================================================================================

            // --------------- CHARACTERISTIC SPEEDS -------------------------------------------
                L1L = uxcn[i-1][j][k]+SOUND;  L2L = uxcn[i-1][j][k]-SOUND;   L3L = uxcn[i-1][j][k];
                L1R = uxcn[i  ][j][k]+SOUND;  L2R = uxcn[i  ][j][k]-SOUND;   L3R = uxcn[i  ][j][k];
            // ---------------------------------------------------------------------------------

            // --------------- LEFT CELL -------------------------------------------------------
            // ----- EXTRAPOLATION FOR THE LEFT GOING WAVES ------------------------------------
                RPNL = 2.0*RCLN-RPL;      QPNL = 2.0*QCLN-QPL;      VPNL = 2.0*VCLN-VPL;    WPNL = 2.0*WCLN-WPL;
            // --------------- RIGHT CELL ----------------------------------------------------
            // ----- EXTRAPOLATION FOR THE RIGHT GOING WAVES -----------------------------------
                RPNR = 2.0*RCRN-RPR;      QPNR = 2.0*QCRN-QPR;      VPNR = 2.0*VCRN-VPR;    WPNR = 2.0*WCRN-WPR;
            // ---------------------------------------------------------------------------------

//      QQ1 = 2.0*(RCLN-RCL) + DT*L1L*(RPL-RL)*HXI;
//      QQ2 = 2.0*(QCLN-QCL) + DT*L2L*(QPL-QL)*HXI;
//      QQ3 = 2.0*(VCLN-VCL) + DT*L3L*(VPL-VL)*HXI;
//      QQ4 = 2.0*(WCLN-WCL) + DT*L3L*(WPL-WL)*HXI;

        QQ1 = 2.0*(RCLN-RCL) + DT*L1L*(RPL-RL)/hc[i];
        QQ2 = 2.0*(QCLN-QCL) + DT*L2L*(QPL-QL)/hc[i];
        QQ3 = 2.0*(VCLN-VCL) + DT*L3L*(VPL-VL)/hc[i];
        QQ4 = 2.0*(WCLN-WCL) + DT*L3L*(WPL-WL)/hc[i];

        RMAX = DMAX1(RL,RCLN,RPL)+QQ1;        RMIN = DMIN1(RL,RCLN,RPL)+QQ1;
        QMAX = DMAX1(QL,QCLN,QPL)+QQ2;        QMIN = DMIN1(QL,QCLN,QPL)+QQ2;
        VMAX = DMAX1(VL,VCLN,VPL)+QQ3;        VMIN = DMIN1(VL,VCLN,VPL)+QQ3;
        WMAX = DMAX1(WL,WCLN,WPL)+QQ4;        WMIN = DMIN1(WL,WCLN,WPL)+QQ4;

        if(RPNL>RMAX) RPNL = RMAX;        if(RPNL<RMIN) RPNL = RMIN;
        if(QPNL>QMAX) QPNL = QMAX;        if(QPNL<QMIN) QPNL = QMIN;
        if(VPNL>VMAX) VPNL = VMAX;        if(VPNL<VMIN) VPNL = VMIN;
        if(WPNL>WMAX) WPNL = WMAX;        if(WPNL<WMIN) WPNL = WMIN;

//      QQ1 = 2.0*(RCRN-RCR) + DT*L1R*(RR-RPR)*HXI;
//      QQ2 = 2.0*(QCRN-QCR) + DT*L2R*(QR-QPR)*HXI;
//      QQ3 = 2.0*(VCRN-VCR) + DT*L3R*(VR-VPR)*HXI;
//      QQ4 = 2.0*(WCRN-WCR) + DT*L3R*(WR-WPR)*HXI;

        QQ1 = 2.0*(RCRN-RCR) + DT*L1R*(RR-RPR)/hc[i];
        QQ2 = 2.0*(QCRN-QCR) + DT*L2R*(QR-QPR)/hc[i];
        QQ3 = 2.0*(VCRN-VCR) + DT*L3R*(VR-VPR)/hc[i];
        QQ4 = 2.0*(WCRN-WCR) + DT*L3R*(WR-WPR)/hc[i];

        RMAX = DMAX1(RR,RCRN,RPR)+QQ1;   RMIN = DMIN1(RR,RCRN,RPR)+QQ1;
        QMAX = DMAX1(QR,QCRN,QPR)+QQ2;   QMIN = DMIN1(QR,QCRN,QPR)+QQ2;
        VMAX = DMAX1(VR,VCRN,VPR)+QQ3;   VMIN = DMIN1(VR,VCRN,VPR)+QQ3;
        WMAX = DMAX1(WR,WCRN,WPR)+QQ4;   WMIN = DMIN1(WR,WCRN,WPR)+QQ4;

        if(RPNR>RMAX) RPNR = RMAX;    if(RPNR<RMIN) RPNR = RMIN;
        if(QPNR>QMAX) QPNR = QMAX;    if(QPNR<QMIN) QPNR = QMIN;
        if(VPNR>VMAX) VPNR = VMAX;    if(VPNR<VMIN) VPNR = VMIN;
        if(WPNR>WMAX) WPNR = WMAX;    if(WPNR<WMIN) WPNR = WMIN;

        if((L1L+L1R)>=0) { RPN1 = RPNL; } else { RPN1 = RPNR; }
        if((L2L+L2R)>=0) { QPN1 = QPNL; } else { QPN1 = QPNR; }
        if((L3L+L3R)>=0) { VPN1 = VPNL;   WPN1 = WPNL; } else { VPN1 = VPNR;   WPN1 = WPNR; }

    // ---------------------------------------------------------------------------------
        RPN2 = 0.5*(RPNL+RPNR);    QPN2 = 0.5*(QPNL+QPNR);  VPN2 = 0.5*(VPNL+VPNR); WPN2 = 0.5*(WPNL+WPNR);
    // ---------------------------------------------------------------------------------
//      QQ1 = 0.5*(2.0*(RCRN-RCR+RCLN-RCL) +DT*0.5*(L1L+L1R)*(RR-RPR+RPL-RL)*HXI);
//      QQ2 = 0.5*(2.0*(QCRN-QCR+QCLN-QCL) +DT*0.5*(L2L+L2R)*(QR-QPR+QPL-QL)*HXI);
//      QQ3 = 0.5*(2.0*(VCRN-VCR+VCLN-VCL) +DT*0.5*(L3L+L3R)*(VR-VPR+VPL-VL)*HXI);
//      QQ4 = 0.5*(2.0*(WCRN-WCR+WCLN-WCL) +DT*0.5*(L3L+L3R)*(WR-WPR+WPL-WL)*HXI);

        QQ1 = 0.5*(2.0*(RCRN-RCR+RCLN-RCL) +DT*0.5*(L1L+L1R)*(RR-RPR+RPL-RL)/hc[i]);
        QQ2 = 0.5*(2.0*(QCRN-QCR+QCLN-QCL) +DT*0.5*(L2L+L2R)*(QR-QPR+QPL-QL)/hc[i]);
        QQ3 = 0.5*(2.0*(VCRN-VCR+VCLN-VCL) +DT*0.5*(L3L+L3R)*(VR-VPR+VPL-VL)/hc[i]);
        QQ4 = 0.5*(2.0*(WCRN-WCR+WCLN-WCL) +DT*0.5*(L3L+L3R)*(WR-WPR+WPL-WL)/hc[i]);

        RMAX = DMAX1(0.5*(RL+RR),0.5*(RCLN+RCRN),0.5*(RPL+RPR))+QQ1;
        RMIN = DMIN1(0.5*(RL+RR),0.5*(RCLN+RCRN),0.5*(RPL+RPR))+QQ1;
        QMAX = DMAX1(0.5*(QL+QR),0.5*(QCLN+QCRN),0.5*(QPL+QPR))+QQ2;
        QMIN = DMIN1(0.5*(QL+QR),0.5*(QCLN+QCRN),0.5*(QPL+QPR))+QQ2;
        VMAX = DMAX1(0.5*(VL+VR),0.5*(VCLN+VCRN),0.5*(VPL+VPR))+QQ3;
        VMIN = DMIN1(0.5*(VL+VR),0.5*(VCLN+VCRN),0.5*(VPL+VPR))+QQ3;
        WMAX = DMAX1(0.5*(WL+WR),0.5*(WCLN+WCRN),0.5*(WPL+WPR))+QQ4;
        WMIN = DMIN1(0.5*(WL+WR),0.5*(WCLN+WCRN),0.5*(WPL+WPR))+QQ4;

        if(RPN2>RMAX) {RPN2 = RMAX;}    if(RPN2<RMIN) {RPN2 = RMIN;}
        if(QPN2>QMAX) {QPN2 = QMAX;}    if(QPN2<QMIN) {QPN2 = QMIN;}
        if(VPN2>VMAX) {VPN2 = VMAX;}    if(VPN2<VMIN) {VPN2 = VMIN;}
        if(WPN2>WMAX) {WPN2 = WMAX;}    if(WPN2<WMIN) {WPN2 = WMIN;}

    // ------------- blend -------------------------------------------------------------
        RPN = Fblend*RPN1+(1.0-Fblend)*RPN2;
        QPN = Fblend*QPN1+(1.0-Fblend)*QPN2;
        VPN = Fblend*VPN1+(1.0-Fblend)*VPN2;
        WPN = Fblend*WPN1+(1.0-Fblend)*WPN2;
    // ------------- blend -------------------------------------------------------------

            // ================================================================================]
            // ------------------ START THE APProxIMATE RIEMANN SOLVER -------------------------
            // ========= fluxes ================================================================
                roxn[i  ][j][k] = RO_INIT*exp(0.5*(RPN-QPN)/SOUND);
                uxxn[i  ][j][k] = 0.5*(RPN+QPN);
                uyxn[i  ][j][k] = VPN;
                uzxn[i  ][j][k] = WPN;

                if(fh_EOS == hmd_EOS_ARGON)
                    pxn[i  ][j][k] = D*(1./3.*A*A*pow(E*roxn[i  ][j][k]-B,3)+C);
                else if(fh_EOS == hmd_EOS_SPC)
                    pxn[i  ][j][k] = roxn[i  ][j][k]*(roxn[i  ][j][k]*A+B)+C;
                else {
                    printf("\nERROR in FHLJModel3D - unknown EOS\n\n");
                    exit(3);
                }

            // ================================================================================]

            // ---------------- For next time step -----------------------------------------------
                RCLN= RCRN;      QCLN= QCRN;        VCLN= VCRN;     WCLN= WCRN;
                RCL = RCR;       QCL = QCR;         VCL = VCR;      WCL = WCR;
                 RL = RPL;        QL = QPL;          VL = VPL;       WL = WPL;
                RPL = RR;        QPL = QR;          VPL = VR;       WPL = WR;
            // ===================================================================================
            }
        }
    }
// ---------------------------------------------------------------------------------------
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void CHARY()
{
// ---------------------------------------------------------------------------------------
    int i,j,k;
    double L1B,L2B,L3B,L1T,L2T,L3T;
    double RCB,RCBN,RCT,RCTN,RB,RPB,RT,RPT,RPNB,RPNT,RPN;
    double QCB,QCBN,QCT,QCTN,QB,QPB,QT,QPT,QPNB,QPNT,QPN;
    double VCB,VCBN,VCT,VCTN,VB,VPB,VT,VPT,VPNB,VPNT,VPN;
    double WCB,WCBN,WCT,WCTN,WB,WPB,WT,WPT,WPNB,WPNT,WPN;
    double QQ1,QQ2,QQ3,QQ4,RMIN,RMAX,QMIN,QMAX,VMIN,VMAX,WMIN,WMAX;
    double RO1,P1,P2,P3;
    double RPN1,QPN1,VPN1,WPN1,RPN2,QPN2,VPN2,WPN2;

    RO1 = 1.0/RO_INIT;

    for(k=1; k<=N3-1; k++)
    {
        for(i=1; i<=N1-1; i++)
        {
            RCBN = uycn[i][0][k]+SOUND*log(rocn[i][0][k]*RO1);        RCB  = uyc [i][0][k]+SOUND*log(roc [i][0][k]*RO1);
            QCBN = uycn[i][0][k]-SOUND*log(rocn[i][0][k]*RO1);        QCB  = uyc [i][0][k]-SOUND*log(roc [i][0][k]*RO1);
            VCBN = uxcn[i][0][k];                                  VCB  = uxc [i][0][k];
            WCBN = uzcn[i][0][k];                                  WCB  = uzc [i][0][k];

            RB   = uyy [i][0][k]+SOUND*log(roy [i][0][k]*RO1);        RPB  = uyy [i][1][k]+SOUND*log(roy [i][1][k]*RO1);
            QB   = uyy [i][0][k]-SOUND*log(roy [i][0][k]*RO1);        QPB  = uyy [i][1][k]-SOUND*log(roy [i][1][k]*RO1);
            VB   = uxy [i][0][k];                                  VPB  = uxy [i][1][k];
            WB   = uzy [i][0][k];                                  WPB  = uzy [i][1][k];

            for(j=1; j<=N2-1; j++)
            {
            // --------------- LOCAL RIEMANN INVARIANTS ----------------------------------------
            // =================================================================================
                P1 = SOUND*log(rocn[i][j  ][k]*RO1);
                P2 = SOUND*log(roc [i][j  ][k]*RO1);
                P3 = SOUND*log(roy [i][j+1][k]*RO1);
            // =================================================================================
                RCTN = uycn[i][j  ][k]+P1;  QCTN = uycn[i][j  ][k]-P1;  VCTN = uxcn[i][j  ][k]; WCTN = uzcn[i][j  ][k];
                RCT  = uyc [i][j  ][k]+P2;  QCT  = uyc [i][j  ][k]-P2;  VCT  = uxc [i][j  ][k]; WCT  = uzc [i][j  ][k];
            // =================================================================================
                RT   = uyy [i][j+1][k]+P3;  QT   = uyy [i][j+1][k]-P3;  VT   = uxy [i][j+1][k]; WT   = uzy [i][j+1][k];
            // =================================================================================
                RPT  = RPB; //uyy [i  ][j  ]+SOUND*log(roy [i  ][j  ]*RO1)
                QPT  = QPB; //uyy [i  ][j  ]-SOUND*log(roy [i  ][j  ]*RO1)
                VPT  = VPB;
                WPT  = WPB;
            // =================================================================================

            // --------------- CHARACTERISTIC SPEEDS -------------------------------------------
                L1B = uycn[i][j-1][k]+SOUND;  L2B = uycn[i][j-1][k]-SOUND;   L3B = uycn[i][j-1][k];
                L1T = uycn[i][j  ][k]+SOUND;  L2T = uycn[i][j  ][k]-SOUND;   L3T = uycn[i][j  ][k];
            // ---------------------------------------------------------------------------------

            // --------------- LEFT CELL -------------------------------------------------------
            // ----- EXTRAPOLATION FOR THE LEFT GOING WAVES ------------------------------------
                RPNB = 2.0*RCBN-RPB;      QPNB = 2.0*QCBN-QPB;      VPNB = 2.0*VCBN-VPB;    WPNB = 2.0*WCBN-WPB;
            // --------------- RIGHT CELL ----------------------------------------------------
            // ----- EXTRAPOLATION FOR THE RIGHT GOING WAVES -----------------------------------
                RPNT = 2.0*RCTN-RPT;      QPNT = 2.0*QCTN-QPT;      VPNT = 2.0*VCTN-VPT;    WPNT = 2.0*WCTN-WPT;
            // ---------------------------------------------------------------------------------
                QQ1 = 2.0*(RCBN-RCB) + DT*L1B*(RPB-RB)*HYI;
        QQ2 = 2.0*(QCBN-QCB) + DT*L2B*(QPB-QB)*HYI;
        QQ3 = 2.0*(VCBN-VCB) + DT*L3B*(VPB-VB)*HYI;
        QQ4 = 2.0*(WCBN-WCB) + DT*L3B*(WPB-WB)*HYI;

        RMAX = DMAX1(RB,RCBN,RPB)+QQ1;        RMIN = DMIN1(RB,RCBN,RPB)+QQ1;
        QMAX = DMAX1(QB,QCBN,QPB)+QQ2;        QMIN = DMIN1(QB,QCBN,QPB)+QQ2;
        VMAX = DMAX1(VB,VCBN,VPB)+QQ3;        VMIN = DMIN1(VB,VCBN,VPB)+QQ3;
        WMAX = DMAX1(WB,WCBN,WPB)+QQ4;        WMIN = DMIN1(WB,WCBN,WPB)+QQ4;

        if(RPNB>RMAX) RPNB = RMAX;        if(RPNB<RMIN) RPNB = RMIN;
        if(QPNB>QMAX) QPNB = QMAX;        if(QPNB<QMIN) QPNB = QMIN;
        if(VPNB>VMAX) VPNB = VMAX;        if(VPNB<VMIN) VPNB = VMIN;
        if(WPNB>WMAX) WPNB = WMAX;        if(WPNB<WMIN) WPNB = WMIN;

        QQ1 = 2.0*(RCTN-RCT) + DT*L1T*(RT-RPT)*HYI;
        QQ2 = 2.0*(QCTN-QCT) + DT*L2T*(QT-QPT)*HYI;
        QQ3 = 2.0*(VCTN-VCT) + DT*L3T*(VT-VPT)*HYI;
        QQ4 = 2.0*(WCTN-WCT) + DT*L3T*(WT-WPT)*HYI;

        RMAX = DMAX1(RT,RCTN,RPT)+QQ1;   RMIN = DMIN1(RT,RCTN,RPT)+QQ1;
        QMAX = DMAX1(QT,QCTN,QPT)+QQ2;   QMIN = DMIN1(QT,QCTN,QPT)+QQ2;
        VMAX = DMAX1(VT,VCTN,VPT)+QQ3;   VMIN = DMIN1(VT,VCTN,VPT)+QQ3;
        WMAX = DMAX1(WT,WCTN,WPT)+QQ4;   WMIN = DMIN1(WT,WCTN,WPT)+QQ4;

        if(RPNT>RMAX) RPNT = RMAX;    if(RPNT<RMIN) RPNT = RMIN;
        if(QPNT>QMAX) QPNT = QMAX;    if(QPNT<QMIN) QPNT = QMIN;
        if(VPNT>VMAX) VPNT = VMAX;    if(VPNT<VMIN) VPNT = VMIN;
        if(WPNT>WMAX) WPNT = WMAX;    if(WPNT<WMIN) WPNT = WMIN;

        if((L1B+L1T)>=0) { RPN1 = RPNB; } else { RPN1 = RPNT; }
        if((L2B+L2T)>=0) { QPN1 = QPNB; } else { QPN1 = QPNT; }
        if((L3B+L3T)>=0) { VPN1 = VPNB;   WPN1 = WPNB; } else { VPN1 = VPNT;   WPN1 = WPNT; }

    // ---------------------------------------------------------------------------------
        RPN2 = 0.5*(RPNB+RPNT);    QPN2 = 0.5*(QPNB+QPNT);    VPN2 = 0.5*(VPNB+VPNT);   WPN2 = 0.5*(WPNB+WPNT);
    // ---------------------------------------------------------------------------------
        QQ1 = 0.5*(2.0*(RCTN-RCT+RCBN-RCB) +DT*0.5*(L1B+L1T)*(RT-RPT+RPB-RB)*HYI);
        QQ2 = 0.5*(2.0*(QCTN-QCT+QCBN-QCB) +DT*0.5*(L2B+L2T)*(QT-QPT+QPB-QB)*HYI);
        QQ3 = 0.5*(2.0*(VCTN-VCT+VCBN-VCB) +DT*0.5*(L3B+L3T)*(VT-VPT+VPB-VB)*HYI);
        QQ4 = 0.5*(2.0*(WCTN-WCT+WCBN-WCB) +DT*0.5*(L3B+L3T)*(WT-WPT+WPB-WB)*HYI);

        RMAX = DMAX1(0.5*(RB+RT),0.5*(RCBN+RCTN),0.5*(RPB+RPT))+QQ1;
        RMIN = DMIN1(0.5*(RB+RT),0.5*(RCBN+RCTN),0.5*(RPB+RPT))+QQ1;
        QMAX = DMAX1(0.5*(QB+QT),0.5*(QCBN+QCTN),0.5*(QPB+QPT))+QQ2;
        QMIN = DMIN1(0.5*(QB+QT),0.5*(QCBN+QCTN),0.5*(QPB+QPT))+QQ2;
        VMAX = DMAX1(0.5*(VB+VT),0.5*(VCBN+VCTN),0.5*(VPB+VPT))+QQ3;
        VMIN = DMIN1(0.5*(VB+VT),0.5*(VCBN+VCTN),0.5*(VPB+VPT))+QQ3;
        WMAX = DMAX1(0.5*(WB+WT),0.5*(WCBN+WCTN),0.5*(WPB+WPT))+QQ4;
        WMIN = DMIN1(0.5*(WB+WT),0.5*(WCBN+WCTN),0.5*(WPB+WPT))+QQ4;

        if(RPN2>RMAX) RPN2 = RMAX;    if(RPN2<RMIN) RPN2 = RMIN;
        if(QPN2>QMAX) QPN2 = QMAX;    if(QPN2<QMIN) QPN2 = QMIN;
        if(VPN2>VMAX) VPN2 = VMAX;    if(VPN2<VMIN) VPN2 = VMIN;
        if(WPN2>WMAX) WPN2 = WMAX;    if(WPN2<WMIN) WPN2 = WMIN;

    // ------------- blend -------------------------------------------------------------
        RPN = Fblend*RPN1+(1.0-Fblend)*RPN2;
        QPN = Fblend*QPN1+(1.0-Fblend)*QPN2;
        VPN = Fblend*VPN1+(1.0-Fblend)*VPN2;
        WPN = Fblend*WPN1+(1.0-Fblend)*WPN2;
    // ------------- blend -------------------------------------------------------------

            // ================================================================================]
            // ------------------ START THE APProxIMATE RIEMANN SOLVER -------------------------
            // ========= fluxes ================================================================
                royn[i][j  ][k] = RO_INIT*exp(0.5*(RPN-QPN)/SOUND);
                uxyn[i][j  ][k] = VPN;
                uzyn[i][j  ][k] = WPN;
                uyyn[i][j  ][k] = 0.5*(RPN+QPN);

                    if(fh_EOS == hmd_EOS_ARGON)
                        pyn[i][j  ][k] = D*(1./3.*A*A*pow(E*royn[i][j  ][k]-B,3)+C);
                    else if(fh_EOS == hmd_EOS_SPC)
                        pyn[i][j  ][k] = royn[i][j  ][k]*(royn[i][j  ][k]*A+B)+C;
                    else {
                        printf("\nERROR in FHLJModel3D - unknown EOS\n\n");
                        exit(3);
                    }

            // ================================================================================]

            // ---------------- For next time step -----------------------------------------------
                RCBN= RCTN;      QCBN= QCTN;        VCBN= VCTN;     WCBN= WCTN;
                RCB = RCT;       QCB = QCT;         VCB = VCT;      WCB = WCT;
                 RB = RPB;        QB = QPB;          VB = VPB;      WB = WPB;
                RPB = RT;        QPB = QT;          VPB = VT;       WPB = WT;
            // ===================================================================================
            }
        }
    }
// ---------------------------------------------------------------------------------------
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void CHARZ()
{
// ---------------------------------------------------------------------------------------
    int i,j,k;
    double L1D,L2D,L3D,L1U,L2U,L3U;
    double RCD,RCDN,RCU,RCUN,RD,RPD,RU,RPU,RPND,RPNU,RPN;
    double QCD,QCDN,QCU,QCUN,QD,QPD,QU,QPU,QPND,QPNU,QPN;
    double VCD,VCDN,VCU,VCUN,VD,VPD,VU,VPU,VPND,VPNU,VPN;
    double WCD,WCDN,WCU,WCUN,WD,WPD,WU,WPU,WPND,WPNU,WPN;
    double QQ1,QQ2,QQ3,QQ4,RMIN,RMAX,QMIN,QMAX,VMIN,VMAX,WMIN,WMAX;
    double RO1,P1,P2,P3;
    double RPN1,QPN1,VPN1,WPN1,RPN2,QPN2,VPN2,WPN2;

    RO1 = 1.0/RO_INIT;

    for(j=1; j<=N2-1; j++)
    {
        for(i=1; i<=N1-1; i++)
        {
            RCDN = uzcn[i][j][0]+SOUND*log(rocn[i][j][0]*RO1);        RCD  = uzc [i][j][0]+SOUND*log(roc [i][j][0]*RO1);
            QCDN = uzcn[i][j][0]-SOUND*log(rocn[i][j][0]*RO1);        QCD  = uzc [i][j][0]-SOUND*log(roc [i][j][0]*RO1);
            VCDN = uxcn[i][j][0];                                  VCD  = uxc [i][j][0];
            WCDN = uycn[i][j][0];                                  WCD  = uyc [i][j][0];

            RD   = uzz [i][j][0]+SOUND*log(roz [i][j][0]*RO1);        RPD  = uzz [i][j][1]+SOUND*log(roz [i][j][1]*RO1);
            QD   = uzz [i][j][0]-SOUND*log(roz [i][j][0]*RO1);        QPD  = uzz [i][j][1]-SOUND*log(roz [i][j][1]*RO1);
            VD   = uxz [i][j][0];                                  VPD  = uxz [i][j][1];
            WD   = uyz [i][j][0];                                  WPD  = uyz [i][j][1];

            for(k=1; k<=N3-1; k++)
            {
            // --------------- LOCAL RIEMANN INVARIANTS ----------------------------------------
            // =================================================================================
                P1 = SOUND*log(rocn[i][j][k  ]*RO1);
                P2 = SOUND*log(roc [i][j][k  ]*RO1);
                P3 = SOUND*log(roz [i][j][k+1]*RO1);
            // =================================================================================
                RCUN = uzcn[i][j][k  ]+P1;  QCUN = uzcn[i][j][k  ]-P1;  VCUN = uxcn[i][j][k  ]; WCUN = uycn[i][j][k  ];
                RCU  = uzc [i][j][k  ]+P2;  QCU  = uzc [i][j][k  ]-P2;  VCU  = uxc [i][j][k  ]; WCU  = uyc [i][j][k  ];
            // =================================================================================
                RU   = uzz [i][j][k+1]+P3;  QU   = uzz [i][j][k+1]-P3;  VU   = uxz [i][j][k+1]; WU   = uyz [i][j][k+1];
            // =================================================================================
                RPU  = RPD; //uyy [i  ][j  ]+SOUND*log(roy [i  ][j  ]*RO1)
                QPU  = QPD; //uyy [i  ][j  ]-SOUND*log(roy [i  ][j  ]*RO1)
                VPU  = VPD;
                WPU  = WPD;
            // =================================================================================

            // --------------- CHARACTERISTIC SPEEDS -------------------------------------------
                L1D = uzcn[i][j][k-1]+SOUND;  L2D = uzcn[i][j][k-1]-SOUND;   L3D = uzcn[i][j][k-1];
                L1U = uzcn[i][j][k  ]+SOUND;  L2U = uzcn[i][j][k  ]-SOUND;   L3U = uzcn[i][j][k  ];
            // ---------------------------------------------------------------------------------

            // --------------- LEFT CELL -------------------------------------------------------
            // ----- EXTRAPOLATION FOR THE LEFT GOING WAVES ------------------------------------
                RPND = 2.0*RCDN-RPD;      QPND = 2.0*QCDN-QPD;      VPND = 2.0*VCDN-VPD;    WPND = 2.0*WCDN-WPD;
            // --------------- RIGHU CELL ----------------------------------------------------
            // ----- EXTRAPOLATION FOR THE RIGHT GOING WAVES -----------------------------------
                RPNU = 2.0*RCUN-RPU;      QPNU = 2.0*QCUN-QPU;      VPNU = 2.0*VCUN-VPU;    WPNU = 2.0*WCUN-WPU;
            // ---------------------------------------------------------------------------------
                QQ1 = 2.0*(RCDN-RCD) +DT*L1D*(RPD-RD)*HZI;
        QQ2 = 2.0*(QCDN-QCD) +DT*L2D*(QPD-QD)*HZI;
        QQ3 = 2.0*(VCDN-VCD) +DT*L3D*(VPD-VD)*HZI;
        QQ4 = 2.0*(WCDN-WCD) +DT*L3D*(WPD-WD)*HZI;

        RMAX = DMAX1(RD,RCDN,RPD)+QQ1;        RMIN = DMIN1(RD,RCDN,RPD)+QQ1;
        QMAX = DMAX1(QD,QCDN,QPD)+QQ2;        QMIN = DMIN1(QD,QCDN,QPD)+QQ2;
        VMAX = DMAX1(VD,VCDN,VPD)+QQ3;        VMIN = DMIN1(VD,VCDN,VPD)+QQ3;
        WMAX = DMAX1(WD,WCDN,WPD)+QQ4;        WMIN = DMIN1(WD,WCDN,WPD)+QQ4;

        if(RPND>RMAX) RPND = RMAX;        if(RPND<RMIN) RPND = RMIN;
        if(QPND>QMAX) QPND = QMAX;        if(QPND<QMIN) QPND = QMIN;
        if(VPND>VMAX) VPND = VMAX;        if(VPND<VMIN) VPND = VMIN;
        if(WPND>WMAX) WPND = WMAX;        if(WPND<WMIN) WPND = WMIN;

        QQ1 = 2.0*(RCUN-RCU) + DT*L1U*(RU-RPU)*HZI;
        QQ2 = 2.0*(QCUN-QCU) + DT*L2U*(QU-QPU)*HZI;
        QQ3 = 2.0*(VCUN-VCU) + DT*L3U*(VU-VPU)*HZI;
        QQ4 = 2.0*(WCUN-WCU) + DT*L3U*(WU-WPU)*HZI;

        RMAX = DMAX1(RU,RCUN,RPU)+QQ1;   RMIN = DMIN1(RU,RCUN,RPU)+QQ1;
        QMAX = DMAX1(QU,QCUN,QPU)+QQ2;   QMIN = DMIN1(QU,QCUN,QPU)+QQ2;
        VMAX = DMAX1(VU,VCUN,VPU)+QQ3;   VMIN = DMIN1(VU,VCUN,VPU)+QQ3;
        WMAX = DMAX1(WU,WCUN,WPU)+QQ4;   WMIN = DMIN1(WU,WCUN,WPU)+QQ4;

        if(RPNU>RMAX) RPNU = RMAX;    if(RPNU<RMIN) RPNU = RMIN;
        if(QPNU>QMAX) QPNU = QMAX;    if(QPNU<QMIN) QPNU = QMIN;
        if(VPNU>VMAX) VPNU = VMAX;    if(VPNU<VMIN) VPNU = VMIN;
        if(WPNU>WMAX) WPNU = WMAX;    if(WPNU<WMIN) WPNU = WMIN;

        if((L1D+L1U)>=0) { RPN1 = RPND; } else { RPN1 = RPNU; }
        if((L2D+L2U)>=0) { QPN1 = QPND; } else { QPN1 = QPNU; }
        if((L3D+L3U)>=0) { VPN1 = VPND;   WPN1 = WPND; } else { VPN1 = VPNU;   WPN1 = WPNU; }

    // ---------------------------------------------------------------------------------
        RPN2 = 0.5*(RPND+RPNU);    QPN2 = 0.5*(QPND+QPNU);    VPN2 = 0.5*(VPND+VPNU);   WPN2 = 0.5*(WPND+WPNU);
    // ---------------------------------------------------------------------------------
        QQ1 = 0.5*(2.0*(RCUN-RCU+RCDN-RCD) + DT*0.5*(L1D+L1U)*(RU-RPU+RPD-RD)*HZI);
        QQ2 = 0.5*(2.0*(QCUN-QCU+QCDN-QCD) + DT*0.5*(L2D+L2U)*(QU-QPU+QPD-QD)*HZI);
        QQ3 = 0.5*(2.0*(VCUN-VCU+VCDN-VCD) + DT*0.5*(L3D+L3U)*(VU-VPU+VPD-VD)*HZI);
        QQ4 = 0.5*(2.0*(WCUN-WCU+WCDN-WCD) + DT*0.5*(L3D+L3U)*(WU-WPU+WPD-WD)*HZI);

        RMAX = DMAX1(0.5*(RD+RU),0.5*(RCDN+RCUN),0.5*(RPD+RPU))+QQ1;
        RMIN = DMIN1(0.5*(RD+RU),0.5*(RCDN+RCUN),0.5*(RPD+RPU))+QQ1;
        QMAX = DMAX1(0.5*(QD+QU),0.5*(QCDN+QCUN),0.5*(QPD+QPU))+QQ2;
        QMIN = DMIN1(0.5*(QD+QU),0.5*(QCDN+QCUN),0.5*(QPD+QPU))+QQ2;
        VMAX = DMAX1(0.5*(VD+VU),0.5*(VCDN+VCUN),0.5*(VPD+VPU))+QQ3;
        VMIN = DMIN1(0.5*(VD+VU),0.5*(VCDN+VCUN),0.5*(VPD+VPU))+QQ3;
        WMAX = DMAX1(0.5*(WD+WU),0.5*(WCDN+WCUN),0.5*(WPD+WPU))+QQ4;
        WMIN = DMIN1(0.5*(WD+WU),0.5*(WCDN+WCUN),0.5*(WPD+WPU))+QQ4;

        if(RPN2>RMAX) RPN2 = RMAX;    if(RPN2<RMIN) RPN2 = RMIN;
        if(QPN2>QMAX) QPN2 = QMAX;    if(QPN2<QMIN) QPN2 = QMIN;
        if(VPN2>VMAX) VPN2 = VMAX;    if(VPN2<VMIN) VPN2 = VMIN;
        if(WPN2>WMAX) WPN2 = WMAX;    if(WPN2<WMIN) WPN2 = WMIN;

    // ------------- blend -------------------------------------------------------------
        RPN = Fblend*RPN1+(1.0-Fblend)*RPN2;
        QPN = Fblend*QPN1+(1.0-Fblend)*QPN2;
        VPN = Fblend*VPN1+(1.0-Fblend)*VPN2;
        WPN = Fblend*WPN1+(1.0-Fblend)*WPN2;
    // ------------- blend -------------------------------------------------------------

            // ================================================================================]
            // ------------------ SUARU UHE APProxIMAUE RIEMANN SOLVER -------------------------
            // ========= fluxes ================================================================
                rozn[i][j][k  ] = RO_INIT*exp(0.5*(RPN-QPN)/SOUND);
                uxzn[i][j][k  ] = VPN;
                uyzn[i][j][k  ] = WPN;
                uzzn[i][j][k  ] = 0.5*(RPN+QPN);

                    if(fh_EOS == hmd_EOS_ARGON)
                        pzn[i][j][k  ] = D*(1./3.*A*A*pow(E*rozn[i][j  ][k]-B,3)+C);
                    else if(fh_EOS == hmd_EOS_SPC)
                        pzn[i][j][k  ] = rozn[i][j  ][k]*(rozn[i][j  ][k]*A+B)+C;
                    else {
                        printf("\nERROR in FHLJModel3D - unknown EOS\n\n");
                        exit(3);
                    }

            // ================================================================================]

            // ---------------- For next time step -----------------------------------------------
                RCDN= RCUN;      QCDN= QCUN;        VCDN= VCUN;     WCDN= WCUN;
                RCD = RCU;       QCD = QCU;         VCD = VCU;      WCD = WCU;
                 RD = RPD;        QD = QPD;          VD = VPD;      WD = WPD;
                RPD = RU;        QPD = QU;          VPD = VU;       WPD = WU;
            // ===================================================================================
            }
        }
    }
// ---------------------------------------------------------------------------------------
}
*/




void define_FH_grid(t_commrec *cr, FHMD *fh)
{
    dvec h0;
    ivec ind;

    for(int d = 0; d < DIM; d++)
        h0[d] = fh->box[d]/(double)(fh->N[d]);

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    fh->grid.h[C][d] = h0[d];
                    fh->grid.n[C][d] = (double)(ind[d])*h0[d];
                    fh->grid.c[C][d] = fh->grid.n[C][d] + 0.5*h0[d];
                }
            }
        }
    }

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->grid.vol[i]  = fh->grid.h[i][0]*fh->grid.h[i][1]*fh->grid.h[i][2];
        fh->grid.ivol[i] = 1.0/fh->grid.vol[i];
    }

#ifdef FHMD_DEBUG_GRID
    if(MASTER(cr)) {
        printf(MAKE_YELLOW "FHMD DEBUG: Grid X (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[0]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(i,0,0,fh->N)][0], fh->grid.c[I3(i,0,0,fh->N)][0], fh->grid.h[I3(i,0,0,fh->N)][0]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[0], fh->grid.n[I3(fh->N[0],0,0,fh->N)][0]);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Y (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[1]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(0,i,0,fh->N)][1], fh->grid.c[I3(0,i,0,fh->N)][1], fh->grid.h[I3(0,i,0,fh->N)][1]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[1], fh->grid.n[I3(0,fh->N[1],0,fh->N)][1]);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Z (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[2]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(0,0,i,fh->N)][2], fh->grid.c[I3(0,0,i,fh->N)][2], fh->grid.h[I3(0,0,i,fh->N)][2]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[2], fh->grid.n[I3(0,0,fh->N[2],fh->N)][2]);
        printf("==============================================\n");
        printf(RESET_COLOR "\n");
    }
#endif
}
