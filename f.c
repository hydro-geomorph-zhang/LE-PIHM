/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : f.c   (PIHM v.2.0)                                            *
 * Function    : Model Kernel: Building ODE system for each physical process   *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * Developer of PIHM v.2.0:  Mukesh Kumar (muk139@psu.edu)		       *
 * Developer of PIHM v.1.0:  Yizhong Qu	(quyizhong@gmail.com)		       *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * NOTE: f.c has gone a massive revamp (essentially rewritten) since PIHM v.1.0*
 *                                                                             *
 *...........MODIFICATIONS/ADDITIONS incorporated in f.c (PIHM v.2.0)...........*
 * a) Surface Flow: 							       *
 *	--> Correction of diffusion wave approximation (calculation of dh/ds)  *
 *              i. Calculation of dh/ds performed using planar slope connecting*
 *                 neighboring centroids				       *
 *              ii.Reflection of elements at boundaries and rivers for dh/ds   *
 *		   calculation
 *	--> Correction of kinematic wave approximation (dh/ds calculation based*
 *	    on elevation only instead of head				       *
 *	--> Correction of gradient for cases with steep change in topography   *
 * b) Subsurface Flow:							       *
 *	--> Addition of macropore phenomena				       *
 *	--> Addition of rectangular cell beneath a river element	       *
 *	--> Implementation of two layered subsurface model(sat/unsat) based on *
 *	Richard's eqn							       *
 *	--> Incorporation of Vertical and Horizontal Anisotropy                *
 *	--> Use of geologic data					       *
 * c) River Flow:							       *
 *	--> Correction of kinematic and diff. wave approximation of SV eqn     *
 *	--> Incorporation of flexible river shapes			       *
 *	--> Separate incorporation of leakage and lateral flow		       *
 *	--> Correction of bank overland flow for extreme cases		       *
 *	--> Addition of aquifer cells below river elements		       *
 * c) Surface/Subsurface Coupling:					       *
 *	--> Implementation of First order coupling through (in/ex)filtration   *
 *		based on head continuity 				       *
 * d) Evaporation:							       *
 *	--> Incorporation of ET from ground/subsurface/vegetation	       *
 *	--> Incorporation of landcover properties for calculation of each ET   *
 *	    component							       *
 * e) Computational:							       *
 *	--> Use of temporary state variables in calculation. Note: Never change*
 *		core state variables					       *
 * f) Miscellaneous (other advantages realtive to PIHM1.0): No maximum         *
 *    constraint on gw level. Accordingly, no numerical constraints on subsur- *
 *    face flux terms.Faster Implementation. Led to first large scale model    *
 *    application.
 *-----------------------------------------------------------------------------*
 *									       *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact 				       *
 *	--> Mukesh Kumar (muk139@psu.edu)				       *
 *	--> Prof. Chris Duffy (cxd11@psu.edu)				       *
 * This code is free for research purpose only.		       		       *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * REFERENCES:                                                                 *
 * PIHM2.0:                                                                    *
 *      a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *      Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *      b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *      a Mesoscale Watershed", Advances in Water Resources (submitted)        *
 * PIHM1.0:                                                                    *
 *      a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *      simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *      b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *      for multiprocess watershed simulation". Water Resources Research       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "nvector_serial.h"
#include "sundials_types.h"
#include "pihm.h"
#define multF	2
#define MINpsi	-70
#define EPS 0.05
#define THRESH 0.0
#define UNIT_C 1440		/* Note 60*24 for calculation of yDot in
				 * m/min units while forcing is in m/day. */
#define GRAV 9.8*60*60		/* Note the dependence on physical units */

#define C_air 1004.0
#define Lv (2.503*pow(10,6))
#define SIGMA (5.67*pow(10,-8)*60)
#define R_dry 287.04
#define R_v 461.5
#define TIME_CONVERT 525600

// #define Niu 0.000001004 //Unit: (m^2/s)

realtype        Interpolation(TSD * Data, realtype t);

// realtype
// returnVal(realtype rArea, realtype rPerem, realtype eqWid, realtype ap_Bool)
// {
	// if (ap_Bool == 1) {
		// return rArea;
	// } else if (ap_Bool == 2) {
		// return rPerem;
	// } else {
		// return eqWid;
	// }
// }

// void
// OverlandFlow(realtype ** flux, int loci, int locj, realtype avg_y, realtype grad_y, realtype avg_sf, realtype crossA, realtype avg_rough)
// {
	// flux[loci][locj] = crossA * pow(avg_y, 2.0 / 3.0) * grad_y / (sqrt(fabs(avg_sf)) * avg_rough);
	// //flux[loci][locj] = (grad_y > 0 ? 1 : -1) * crossA * pow(avg_y, 2.0 / 3.0) * sqrt(fabs(grad_y)) / (avg_rough);
// }

/* Note: above formulation doesnt' satisfies upwind/downwind scheme. */
realtype
avgY(realtype diff, realtype yi, realtype yinabr)
{
	if (diff > 0) {
		if (yi > 1 * EPS / 100) {
			//return 0.5 * (yi + yinabr);
			//return ((yinabr > yi) ? 0 : 1.0 * yi);
			/*
			 * Note the if-else TRUE case can be possible only
			 * for Kinematic case
			 */
			return 1.0 * yi;
		} else {
			return 0;
		}
	} else {
		if (yinabr > 1 * EPS / 100) {
			//return 0.5 * (yi + yinabr);
			//return ((yi > yinabr) ? 0 : 1.0 * yinabr);
			/*
			 * Note the if-else TRUE case can be possible only
			 * for Kinematic case
			 */
			return 1.0 * yinabr;
		} else {
			return 0;
		}
	}
}
realtype
effKV(realtype ksatFunc, realtype gradY, realtype macKV, realtype KV, realtype areaF)
{
	if (ksatFunc >= 0.98) {
		return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
	} else {
		if (fabs(gradY) * ksatFunc * KV <= 1 * KV * ksatFunc) {
			return KV * ksatFunc;
		} else {
			if (fabs(gradY) * ksatFunc * KV < (macKV * areaF + KV * (1 - areaF) * ksatFunc)) {
				return (macKV * areaF * ksatFunc + KV * (1 - areaF) * ksatFunc);
			} else {
				return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
			}
		}
	}
}
realtype
effKH(int mp, realtype tmpY, realtype aqDepth, realtype MacD, realtype MacKsatH, realtype areaF, realtype ksatH)
{
	if (mp == 1) {
		if (tmpY > aqDepth - MacD) {
			if (tmpY > aqDepth) {
				return (MacKsatH * MacD * areaF + ksatH * (aqDepth - MacD * areaF)) / aqDepth;
			} else {
				return (MacKsatH * (tmpY - (aqDepth - MacD)) * areaF + ksatH * (aqDepth - MacD + (tmpY - (aqDepth - MacD)) * (1 - areaF))) / tmpY;
			}
		} else {
			return ksatH;
		}
	} else {
		return ksatH;
	}
}
realtype* initiation(Model_Data MD,realtype *DY,realtype *Y, realtype t)
{
	int             i, j, k, inabr;
	for (i = 0; i < MD->NumEle; i++) 
	{
		//MD->DummyY[i] = (Y[i] >= 0) ? Y[i] : 0;
		for(j=0;j<3;j++)
		{
			MD->FluxSub[i][j] = 0;
			MD->FluxSurf[i][j] = 0;
			MD->Seepage[i][j] = 0;
			if (MD->Keyword_LE_PIHM == 1)
			{
				MD->BedloadFlux[i][j] = 0;
				MD->SoilFluxSub[i][j] = 0;
			}
		}
		MD->EleViR[i] = 0;
		MD->Recharge[i] = 0;
		MD->EleET[i][1] = 0;
		MD->EleET[i][2] = 0;
		MD->EleET[i][3] = 0;
		MD->EleET[i][4] = 0;
		MD->EleET[i][5] = 0;
		if (MD->Keyword_LE_PIHM == 1)
		{
			MD->Soilproduction[i] = 0;
			MD->Uplift[i] = 0;
		}
		
		
		MD->DummyY[i] = Y[i];
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			MD->DummyY[i + MD->NumEle] = Y[i + MD->NumEle];
			MD->DummyY[i + 2*MD->NumEle] = Y[i + 2*MD->NumEle];
			if (MD->Keyword_LE_PIHM == 1)
			{
				MD->DummyY[i + 3*MD->NumEle] = Y[i + 3*MD->NumEle];
				MD->DummyY[i + 4*MD->NumEle] = Y[i + 4*MD->NumEle];
			}
		}
		else
		{
			if (MD->Keyword_LE_PIHM == 1)
			{
				MD->DummyY[i + MD->NumEle] = Y[i + MD->NumEle];
				MD->DummyY[i + 2*MD->NumEle] = Y[i + 2*MD->NumEle];
			}
		}		
		
		DY[i] = 0;
		MD->Total_Surfwater_out[i] = 0;
		MD->Total_Unsatwater_out[i] = 0;
		MD->Total_Satwater_out[i] = 0;
		MD->infil_mode[i] = 0;
		MD->Ele[i].Tan_slope = 0.00000000;
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			DY[i + MD->NumEle] = 0;
			DY[i + 2*MD->NumEle] = 0;
			if (MD->Keyword_LE_PIHM == 1){
				MD->Total_soil_in[i] = 0;
				MD->Total_soil_out[i] = 0;
				DY[i + 3*MD->NumEle] = 0;
				DY[i + 4*MD->NumEle] = 0;
			}
		}
		else
		{
			if (MD->Keyword_LE_PIHM == 1)
			{
				MD->Total_soil_in[i] = 0;
				MD->Total_soil_out[i] = 0;
				DY[i + MD->NumEle] = 0;
				DY[i + 2*MD->NumEle] = 0;
			}
		}
		
		
	}
	
	for (i = 0; i < MD->NumEle; i++) 
	{
		if ((MD->SurfMode == 2) && (i < MD->NumEle)) 
		{
			for (j = 0; j < 3; j++) 
			{
				inabr = MD->Ele[i].nabr[j] - 1;
				MD->Ele[i].surfH[j] = (inabr >=0) ? (MD->Ele[inabr].zmax + MD->DummyY[inabr]) : ((MD->Ele[i].BC[j] != 1) ? (MD->Ele[i].zmax + MD->DummyY[i]) : Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
				if (MD->Keyword_LE_PIHM == 1)
				{
					if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
					{
						MD->Ele[i].surfH[j] = (inabr >=0) ? (MD->DummyY[inabr+3*MD->NumEle] + MD->DummyY[inabr]) : ((MD->Ele[i].BC[j] != 1) ? (MD->DummyY[i + 3*MD->NumEle] + MD->DummyY[i]) : Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
					}
					else
					{
						MD->Ele[i].surfH[j] = (inabr >=0) ? (MD->DummyY[inabr+MD->NumEle] + MD->DummyY[inabr]) : ((MD->Ele[i].BC[j] != 1) ? (MD->DummyY[i + MD->NumEle] + MD->DummyY[i]) : Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
					}
				}
			}
			MD->Ele[i].dhBYdx = -1 * (MD->Ele[i].surfY[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfX[2] * (MD->Ele[i].surfY[1] - MD->Ele[i].surfY[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfY[0] - MD->Ele[i].surfY[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfY[2] - MD->Ele[i].surfY[1]));
			MD->Ele[i].dhBYdy = -1 * (MD->Ele[i].surfX[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfY[2] * (MD->Ele[i].surfX[1] - MD->Ele[i].surfX[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfX[0] - MD->Ele[i].surfX[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfX[2] - MD->Ele[i].surfX[1]));
		}
	}
	return DY;
}


void HydroSubLateral(Model_Data MD, int i,int j,int inabr)
{
	//int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	
	if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		if (inabr>=0) 
		{
			//----------------Set up aquifer depth----------------------//
			AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
			if (MD->Keyword_LE_PIHM == 1)
			{
				AquiferDepth = (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[i + 4*MD->NumEle]);
			}
			
			nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
			if (MD->Keyword_LE_PIHM == 1)
			{
				nabrAqDepth = (MD->DummyY[inabr+3*MD->NumEle] - MD->DummyY[inabr+4*MD->NumEle]);
			}					
			
			if (AquiferDepth < MD->Ele[i].macD)
			{
				if (AquiferDepth<0)
				{
					MD->Ele[i].macD = 0;
					AquiferDepth = 0;
				}
				else
				{
					MD->Ele[i].macD = AquiferDepth;
				}
			}
			
			if (nabrAqDepth < MD->Ele[inabr].macD)
			{
				if (nabrAqDepth<0)
				{
					MD->Ele[inabr].macD = 0;
					nabrAqDepth = 0;
				}
				else
				{
					MD->Ele[inabr].macD = nabrAqDepth;
				}
			}
			//----------------Finish Setting up aquifer depth----------------------//
			
			//----------------Calculation----------------------//
			if (AquiferDepth>0 && nabrAqDepth>0)
			{
				Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[inabr + 2 * MD->NumEle] + MD->Ele[inabr].zmin);
				if (MD->Keyword_LE_PIHM == 1)
				{
					Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->DummyY[i + 3*MD->NumEle]-AquiferDepth) - (MD->DummyY[inabr + 2 * MD->NumEle] + MD->DummyY[inabr+3*MD->NumEle]-nabrAqDepth);
				}					
												
				Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], MD->DummyY[inabr + 2 * MD->NumEle]);
				Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
				Grad_Y_Sub = Dif_Y_Sub / Distance;
				/* take care of macropore effect */
				effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
				effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
				Avg_Ksat = 0.5 * (effK + effKnabr);
				/* groundwater flow modeled by Darcy's law */
				MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
				MD->Seepage[i][j] = 0;
			}
			else if (AquiferDepth>0 && nabrAqDepth==0)
			{
				Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[inabr] + MD->Ele[inabr].zmin);
				if (MD->Keyword_LE_PIHM == 1)
				{
					Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->DummyY[i + 4*MD->NumEle]) - (MD->DummyY[inabr] + MD->DummyY[inabr+4*MD->NumEle]);
				}																
				Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], MD->DummyY[inabr]);
				Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
				Grad_Y_Sub = Dif_Y_Sub / Distance;
				/* take care of macropore effect */
				effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
				Avg_Ksat = effK;
				/* groundwater flow modeled by Darcy's law */
				MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
				MD->Seepage[i][j] = 0;
			}
			else if (AquiferDepth==0 && nabrAqDepth>0)
			{
				Dif_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmin) - (MD->DummyY[inabr + 2 * MD->NumEle] + MD->Ele[inabr].zmin);
				if (MD->Keyword_LE_PIHM == 1)
				{
					Dif_Y_Sub = (MD->DummyY[i] + MD->DummyY[i + 4*MD->NumEle]) - (MD->DummyY[inabr + 2 * MD->NumEle] + MD->DummyY[inabr+4*MD->NumEle]);
				}					
												
				Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i], MD->DummyY[inabr + 2 * MD->NumEle]);
				Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
				Grad_Y_Sub = Dif_Y_Sub / Distance;
				/* take care of macropore effect */
				effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
				Avg_Ksat = effKnabr;
				/* groundwater flow modeled by Darcy's law */
				MD->FluxSub[i][j] = 0;
				MD->Seepage[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
			}
			else if (AquiferDepth==0 && nabrAqDepth==0)
			{
				MD->FluxSub[i][j] = 0;
				MD->Seepage[i][j] = 0;
			}
			//----------------Finish Calculation----------------------//
		}
		else
		{
			MD->FluxSub[i][j] = 0;
			MD->Seepage[i][j] = 0;
		}
		
		if (MD->FluxSub[i][j]>0)
			MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] + MD->FluxSub[i][j]/MD->Ele[i].area;
		if (MD->Seepage[i][j]>0)
			MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->Seepage[i][j]/MD->Ele[i].area;
	}
	else
	{
		MD->FluxSub[i][j] = 0;
		MD->Seepage[i][j] = 0;
	}
}

void HydroSurfLateral(Model_Data MD, int i,int j,int inabr)
{
	//int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	
	//----------------Calculation----------------------//
	if (inabr>=0) 
	{
		Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->Ele[i].zmax - MD->Ele[inabr].zmax) : (MD->DummyY[i] + MD->Ele[i].zmax) - (MD->DummyY[inabr] + MD->Ele[inabr].zmax);
		if (MD->Keyword_LE_PIHM == 1)
		{
			if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
			{
				Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[inabr+3*MD->NumEle]) : (MD->DummyY[i] + MD->DummyY[i + 3*MD->NumEle]) - (MD->DummyY[inabr] + MD->DummyY[inabr+3*MD->NumEle]);
			}
			else
			{
				Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->DummyY[i + MD->NumEle] - MD->DummyY[inabr+MD->NumEle]) : (MD->DummyY[i] + MD->DummyY[i + MD->NumEle]) - (MD->DummyY[inabr] + MD->DummyY[inabr+MD->NumEle]);
			}
		}
		Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], MD->DummyY[inabr]);
		Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
		Grad_Y_Surf = Dif_Y_Surf / Distance;
		Avg_Sf = 0.5*(sqrt(pow(MD->Ele[i].dhBYdx, 2) + pow(MD->Ele[i].dhBYdy, 2))+sqrt(pow(MD->Ele[inabr].dhBYdx, 2) + pow(MD->Ele[inabr].dhBYdy, 2)));//?? Xuan Weighting needed
		//updated in 2.2
		Avg_Sf = (MD->SurfMode == 1) ? (Grad_Y_Surf > 0 ? Grad_Y_Surf : EPS / pow(10.0, 6)) : (Avg_Sf > EPS / pow(10.0, 6)) ? Avg_Sf : EPS / pow(10.0, 6);
		/* Weighting needed */
		Avg_Rough = 0.5 * (MD->Ele[i].Rough + MD->Ele[inabr].Rough);
		CrossA = Avg_Y_Surf * MD->Ele[i].edge[j];
		MD->FluxSurf[i][j] = CrossA * pow(Avg_Y_Surf, 2.0 / 3.0) * Grad_Y_Surf / (sqrt(fabs(Avg_Sf)) * Avg_Rough);
	}
	else
	{
		if (MD->Ele[i].BC[j] == 0) 
		{
			MD->FluxSurf[i][j] = 0;
		} 
		else if (MD->Ele[i].BC[j] == 3)
		{
			Dif_Y_Surf = MD->DummyY[i];
			if (MD->Keyword_LE_PIHM == 1)
			{
				if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
				{
					Dif_Y_Surf = (MD->DummyY[i]+MD->DummyY[i + 3*MD->NumEle] - MD->zmax_init[i]);
				}
				else
				{
					Dif_Y_Surf = (MD->DummyY[i]+MD->DummyY[i + MD->NumEle] - MD->zmax_init[i]);
				}
			}
			Avg_Y_Surf = MD->DummyY[i];
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
			Grad_Y_Surf = Dif_Y_Surf/Distance;
			Avg_Rough = MD->Ele[i].Rough;
			CrossA = MD->DummyY[i]*MD->Ele[i].edge[j];        
			MD->FluxSurf[i][j] = sqrt((Grad_Y_Surf>0)?Grad_Y_Surf:0)*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough; 														
		}				
	}
	//----------------Finish Calculation----------------------//
	if (MD->FluxSurf[i][j]>0)
	{
		MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->FluxSurf[i][j]/MD->Ele[i].area;
	}
	
}

void SedDiffu(Model_Data MD, int i,int j,int inabr)
{
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	
	
	//----------------Set up aquifer depth----------------------//
	
	AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
	if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		AquiferDepth = (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[i + 4*MD->NumEle]);
	}
	else
	{
		AquiferDepth = (MD->DummyY[i + MD->NumEle] - MD->DummyY[i + 2*MD->NumEle]);
	}
	//----------------Finish Setting up aquifer depth----------------------//
	
	//----------------Calculation----------------------//
	if (AquiferDepth>0)
	{
		if (inabr>=0) 
		{
			if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
			{
				Dif_Elev = (MD->DummyY[i + 3 * MD->NumEle] - MD->DummyY[inabr + 3 * MD->NumEle]);
			}
			else
			{
				Dif_Elev = (MD->DummyY[i + MD->NumEle] - MD->DummyY[inabr + MD->NumEle]);
			}
			Grad_Elev = Dif_Elev/Distance;
			if (MD->Ele[i].Tan_slope<fabs(Grad_Elev))
			{
				MD->Ele[i].Tan_slope = fabs(Grad_Elev);
			}
			MD->SoilFluxSub[i][j] = 0.5*(MD->Ele[i].DReg + MD->Ele[inabr].DReg)*Grad_Elev*MD->Ele[i].edge[j]*MD->MSF/525600;
		}
		else
		{
			if (MD->Ele[i].BC[j] == 0) 
			{
				MD->SoilFluxSub[i][j] = 0;
			} 
			else if (MD->Ele[i].BC[j] == 3)
			{
				if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
				{
					Dif_Elev = (MD->DummyY[i + 3 * MD->NumEle] - MD->zmax_init[i]);
				}
				else
				{
					Dif_Elev = (MD->DummyY[i + MD->NumEle] - MD->zmax_init[i]);
				}
				Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));			
				Grad_Elev = Dif_Elev/Distance;
				if (MD->Ele[i].Tan_slope<fabs(Grad_Elev))
				{
					MD->Ele[i].Tan_slope = fabs(Grad_Elev);
				}
				MD->SoilFluxSub[i][j] = (MD->Ele[i].DReg)*((Grad_Elev>0)?Grad_Elev:0)*MD->Ele[i].edge[j]*MD->MSF/525600;														
			}
		}
	}
	else
	{
		MD->SoilFluxSub[i][j] = 0;
	}
	
	if (MD->SoilFluxSub[i][j]>0)
	{
		MD->Total_soil_out[i] = MD->Total_soil_out[i] + MD->SoilFluxSub[i][j]/MD->Ele[i].area;
	}
}
void ET(Model_Data MD, int i,realtype t)
{
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	realtype		ETa;
	realtype		tempET;
	
	
	if (MD->Keyword_LE_PIHM_HYDRO == 1)
	{
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			//----------------Set up aquifer depth----------------------//
			AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
			if (MD->Keyword_LE_PIHM == 1)
			{
				AquiferDepth = (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[i + 4*MD->NumEle]);
			}
			Deficit = AquiferDepth - MD->DummyY[i + 2 * MD->NumEle];
			
			if (AquiferDepth < MD->Ele[i].RzD)
			{
				if (AquiferDepth<0)
				{
					MD->Ele[i].RzD = 0;
					AquiferDepth = 0;
				}
				else
				{
					MD->Ele[i].RzD = AquiferDepth;
				}
			}
			//----------------Finish Setting up aquifer depth----------------------//
			
			//----------------Calculation----------------------//
			Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn - 1], t);
			//G = Interpolation(&MD->TSD_G[MD->Ele[i].G - 1], t);
			G = 0.1 * Rn;
			T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp - 1], t);
			Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel - 1], t);
			RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity - 1], t);
			VP = 611.2 * exp(17.67 * T / (T + 243.5)) * RH;
			P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
			if (MD->Keyword_LE_PIHM == 1)
			{
				P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->DummyY[i + 3 * MD->NumEle]) / 293, 5.26);
			}
			qv = 0.622 * VP / P;
			qv_sat = 0.622 * (VP / RH) / P;
			//P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
			//Delta = 2503 * pow(10, 3) * exp(17.27 * T / (T + 237.3)) / (pow(237.3 + T, 2));
			//Gamma = P * 1.0035 * 0.92 / (0.622 * 2441);
			LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC - 1], t);
			/*
			* zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
			* cnpy_h =
			* zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));
			* if(LAI<2.85)	{ rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
			* } else { rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h)); }
			*/
			rl = Interpolation(&MD->TSD_RL[MD->Ele[i].LC - 1], t);
			r_a = 12 * 4.72 * log(MD->Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;

			Gamma = 4 * 0.7 * SIGMA * UNIT_C * R_dry / C_air * pow(T + 273.15, 4) / (P / r_a) + 1;
			Delta = Lv * Lv * 0.622 / R_v / C_air / pow(T + 273.15, 2) * qv_sat;
			ETp = (Rn * Delta + Gamma * (1.2 * Lv * (qv_sat - qv) / r_a)) / (1000.0 * Lv * (Delta + Gamma));
			//MD->EleET[i][2] = MD->pcCal.Et2 * (1 - MD->Ele[i].VegFrac) * (Rn * (1 - MD->Ele[i].Albedo) * Delta + (1.2 * 1003.5 * ((VP / RH) - VP) / r_a)) / (1000.0 * 2441000.0 * (Delta + Gamma));
			// BHATT: MAJOR BUG = AQUIFER DEPTH NOT CALCULATED EARLIER
			if (Deficit <= MD->Ele[i].RzD) 
			{
				elemSatn = 1.0;
			} 
			else 
			{
				elemSatn = ((MD->DummyY[i + MD->NumEle] / (AquiferDepth - MD->DummyY[i + 2 * MD->NumEle])) > 1) ? 1 : ((MD->DummyY[i + MD->NumEle] / (AquiferDepth - MD->DummyY[i + 2 * MD->NumEle])) < 0) ? 0 : 0.5 * (1 - cos(3.14 * (MD->DummyY[i + MD->NumEle] / (AquiferDepth - MD->DummyY[i + 2 * MD->NumEle]))));
			}
			ThetaRef = 0.7 * MD->Soil[(MD->Ele[i].soil - 1)].ThetaS;
			ThetaW = 1.05 * MD->Soil[(MD->Ele[i].soil - 1)].ThetaR;
			beta_s = (elemSatn * MD->Ele[i].Porosity + MD->Soil[(MD->Ele[i].soil - 1)].ThetaR - ThetaW) / (ThetaRef - ThetaW);
			beta_s = (beta_s < 0.0001) ? 0.0001 : (beta_s > 1 ? 1 : beta_s);
			ETa = MD->pcCal.Et2 * (1 - MD->Ele[i].VegFrac) * beta_s * ETp;
			if (LAI > 0.0) 
			{
				Rmax = 5000.0 / (60 * UNIT_C);	/* Unit day_per_m */
				f_r = 1.1 * 1.5 * Rn / (MD->Ele[i].Rs_ref * LAI);
				f_r = f_r < 0 ? 0 : f_r;
				alpha_r = (1 + f_r) / (f_r + (MD->Ele[i].Rmin / Rmax));
				alpha_r = alpha_r > 10000 ? 10000 : alpha_r;
				eta_s = 1 - 0.0016 * (pow((24.85 - T), 2));
				eta_s = eta_s < 0.0001 ? 0.0001 : eta_s;
				gamma_s = 1 / (1 + 0.00025 * (VP / RH - VP));
				gamma_s = (gamma_s < 0.01) ? 0.01 : gamma_s;
				r_s = ((MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s)) > Rmax) ? Rmax : (MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s));
				P_c = (1 + Delta / Gamma) / (1 + r_s / r_a + Delta / Gamma);
				tempET = MD->pcCal.Et1 * MD->Ele[i].VegFrac * P_c * (1 - pow(((MD->EleIS[i] + MD->EleSnowCanopy[i] < 0) ? 0 : (MD->EleIS[i] + MD->EleSnowCanopy[i])) / (MD->EleISmax[i] + MD->EleISsnowmax[i]), 1.0 / 2.0)) * ETp;
				tempET = tempET < 0 ? 0 : tempET;
			} 
			else 
			{
				tempET = 0.0;
			}
			
			if (AquiferDepth > 0)
			{
				if (MD->DummyY[i]>ETa)
				{
					MD->EleET[i][1] = ETa;
					MD->EleET[i][2] = 0;
					MD->EleET[i][3] = 0;
				}
				else
				{
					MD->EleET[i][1] = MD->DummyY[i];
					if (Deficit>MD->Ele[i].infD)
					{
						if (MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity>ETa-MD->DummyY[i])
						{
							MD->EleET[i][2] = ETa-MD->DummyY[i];
							MD->EleET[i][3] = 0;
						}
						else
						{
							MD->EleET[i][2] = MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity;
							MD->EleET[i][3] = 0;
						}
					}
					else
					{
						if (MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity>ETa-MD->DummyY[i])
						{
							MD->EleET[i][2] = 0;
							MD->EleET[i][3] = ETa-MD->DummyY[i];
						}
						else
						{
							MD->EleET[i][2] = 0;
							MD->EleET[i][3] = MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity;
						}
					}
				}
				
			
				if (Deficit>MD->Ele[i].RzD)
				{
					if (MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity>tempET)
					{
						MD->EleET[i][4] = tempET;
						MD->EleET[i][5] = 0;
					}
					else
					{
						MD->EleET[i][4] = MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity;
						MD->EleET[i][5] = 0;
					}
				}
				else
				{
					if (MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity>tempET)
					{
						MD->EleET[i][4] = 0;
						MD->EleET[i][5] = tempET;
					}
					else
					{
						MD->EleET[i][4] = 0;
						MD->EleET[i][5] = MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity;
					}
				}
			}
			else
			{
				if (MD->DummyY[i]>ETa)
				{
					MD->EleET[i][1] = ETa;
					MD->EleET[i][2] = 0;
					MD->EleET[i][3] = 0;
					MD->EleET[i][4] = 0;
					MD->EleET[i][5] = 0;
				}
				else
				{
					MD->EleET[i][1] = MD->DummyY[i];
					MD->EleET[i][2] = 0;
					MD->EleET[i][3] = 0;
					MD->EleET[i][4] = 0;
					MD->EleET[i][5] = 0;
				}
			}
		}
		else
		{
			MD->EleET[i][1] = 0;
			MD->EleET[i][2] = 0;
			MD->EleET[i][3] = 0;
			MD->EleET[i][4] = 0;
			MD->EleET[i][5] = 0;
		}
	}
	else
	{
		MD->EleET[i][1] = 0;
		MD->EleET[i][2] = 0;
		MD->EleET[i][3] = 0;
		MD->EleET[i][4] = 0;
		MD->EleET[i][5] = 0;
	}
	MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->EleET[i][1];
	MD->Total_Unsatwater_out[i] = MD->Total_Unsatwater_out[i] + MD->EleET[i][2] + MD->EleET[i][4];
	MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] + MD->EleET[i][3] + MD->EleET[i][5];
}



void Infiltration(Model_Data MD, int i)
{
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	realtype		ETa;
	realtype		temp;
	
	
	if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		//----------------Set up aquifer depth----------------------//
		AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
		if (MD->Keyword_LE_PIHM == 1)
		{
			AquiferDepth = (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[i + 4*MD->NumEle]);
		}
		Deficit = AquiferDepth - MD->DummyY[i + 2 * MD->NumEle];
			
		if (AquiferDepth < MD->Ele[i].infD)
		{
			if (AquiferDepth<0)
			{
				MD->Ele[i].infD = 0;
				AquiferDepth = 0;
			}
			else
			{
				MD->Ele[i].infD = AquiferDepth;
			}
		}
		//----------------Finish Setting up aquifer depth----------------------//
		
		//----------------Calculation----------------------//
		if (AquiferDepth>0)
		{
			if (Deficit>MD->Ele[i].infD)
			{
				if (MD->DummyY[i + MD->NumEle]<Deficit)
				{
					MD->infil_mode[i] = 1;
					elemSatn = MD->DummyY[i + MD->NumEle]/Deficit;
					elemSatn = (elemSatn < multF * EPS) ? multF * EPS : elemSatn;
					Avg_Y_Sub = (-(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha) < MINpsi) ? MINpsi : -(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha);
								
					TotalY_Ele = Avg_Y_Sub + MD->Ele[i].zmax - Deficit;
					Grad_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmax - TotalY_Ele) / Deficit;
						
					Grad_Y_Sub = ((MD->DummyY[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
					satKfunc = pow(elemSatn, 0.5) * pow(-1 + pow(1 - pow(elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
					effK = (MD->Ele[i].Macropore == 1) ? effKV(satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].infKsatV, MD->Ele[i].hAreaF) : MD->Ele[i].infKsatV;
					if (Avg_Y_Sub < 0)
					{
						MD->EleViR[i] = effK * Grad_Y_Sub;
					}
					else 
					{
						MD->EleViR[i] = 0;
					}
								
					effK = (MD->Ele[i].Macropore == 1) ? ((MD->DummyY[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].macD) ? effK : MD->Ele[i].KsatV * satKfunc) : MD->Ele[i].KsatV * satKfunc;
					MD->Recharge[i] = (elemSatn == 0.0) ? 0 : (Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->DummyY[i + 2 * MD->NumEle] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow(-1 + pow(elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow(Deficit + MD->DummyY[i + 2 * MD->NumEle], 2));
					MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->DummyY[i + MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
					MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->DummyY[i + 2 * MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
				}
				else if(MD->DummyY[i + MD->NumEle]>=Deficit)
				{
					MD->infil_mode[i] = 1;
					elemSatn = 1;
					MD->EleViR[i] = 0;
								
					effK = (MD->Ele[i].Macropore == 1) ? ((MD->DummyY[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].macD) ? effK : MD->Ele[i].KsatV * satKfunc) : MD->Ele[i].KsatV * satKfunc;
					MD->Recharge[i] = (elemSatn == 0.0) ? 0 : (Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->DummyY[i + 2 * MD->NumEle] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow(-1 + pow(elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow(Deficit + MD->DummyY[i + 2 * MD->NumEle], 2));
					MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->DummyY[i + MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
					MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->DummyY[i + 2 * MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
				}	
			}
			else
			{
				MD->infil_mode[i] = -1;
				Grad_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmax - (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin)) / AquiferDepth;
				Grad_Y_Sub = ((MD->DummyY[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
				elemSatn = 1.0;
				satKfunc = pow(elemSatn, 0.5) * pow(-1 + pow(1 - pow(elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
				effK = (MD->Ele[i].Macropore == 1) ? effKV(satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].infKsatV, MD->Ele[i].hAreaF) : MD->Ele[i].infKsatV;
							
				MD->EleViR[i] = effK * Grad_Y_Sub;
				MD->Recharge[i] = MD->EleViR[i];	
			}
		}
		else
		{
			MD->EleViR[i] = 0;
			MD->Recharge[i] = 0;
		}
	}
	else
	{
		MD->EleViR[i] = 0;
		MD->Recharge[i] = 0;
	}
	//----------------Calculation----------------------//
	if (MD->infil_mode[i] == -1)
	{
		if (MD->EleViR[i]>0)
		{
			MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->EleViR[i];
		}
		
		if (MD->EleViR[i]<0)
		{
			MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] - MD->EleViR[i];
		}
	}
	
	if (MD->infil_mode[i] == 1)
	{
		if (MD->EleViR[i]>0)
		{
			MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->EleViR[i];
		}
		
		if (MD->Recharge[i]>0)
		{
			MD->Total_Unsatwater_out[i] = MD->Total_Unsatwater_out[i] + MD->Recharge[i];
		}
		
		if (MD->Recharge[i]<0)
		{
			MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] - MD->Recharge[i];
		}
	}
	
	
	
}

void Adjust_water(Model_Data MD)
{
	int	i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	
	for (i = 0; i < MD->NumEle; i++) 
	{
		//-----------Surface water------------------//
		if (MD->Total_Surfwater_out[i]>0&&MD->Total_Surfwater_out[i]/UNIT_C>MD->DummyY[i]&&MD->DummyY[i]>=0)
		{
			for (j=0;j<3;j++)
			{
				inabr = MD->Ele[i].nabr[j] - 1; 
				
				if (MD->FluxSurf[i][j]>0)
				{
					MD->FluxSurf[i][j] = ((MD->FluxSurf[i][j]/MD->Ele[i].area)/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C*MD->Ele[i].area;
					if (inabr>=0)
						MD->FluxSurf[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSurf[i][j];
				}
				if (MD->Seepage[i][j]>0)
				{
					MD->Seepage[i][j] = ((MD->Seepage[i][j]/MD->Ele[i].area)/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C*MD->Ele[i].area;
					if (inabr>=0)
						MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->Seepage[i][j];
				}	
			}
			if (MD->EleET[i][1]>0)
			{
				MD->EleET[i][1] = ((MD->EleET[i][1])/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C;
			}
			if (MD->EleViR[i]>0)
			{
				MD->EleViR[i] = ((MD->EleViR[i])/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C;
				if (MD->infil_mode[i]==-1)
				{
					MD->Recharge[i] = MD->EleViR[i];
				}
			}
			
		}
		else if(MD->DummyY[i]<0)
		{
			for (j=0;j<3;j++)
			{
				inabr = MD->Ele[i].nabr[j] - 1; 
				if (MD->FluxSurf[i][j]>0)
				{
					MD->FluxSurf[i][j] = 0;
					if (inabr>=0)
						MD->FluxSurf[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSurf[i][j];
				}
				if (MD->Seepage[i][j]>0)
				{
					MD->Seepage[i][j] = 0;
					if (inabr>=0)
						MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->Seepage[i][j];
				}
			}
			if (MD->EleET[i][1]>0)
			{
				MD->EleET[i][1] = 0;
			}
			if (MD->EleViR[i]>0)
			{
				MD->EleViR[i] = 0;
				if (MD->infil_mode[i]==-1)
				{
					MD->Recharge[i] = MD->EleViR[i];
				}
			}
		}
		//-----------Finish surface water------------------//
		

		if (MD->Keyword_LE_PIHM_SUBHYDRO==1)
		{
			//-----------Unsat water------------------//
			if (MD->Total_Unsatwater_out[i]>0 && MD->Total_Unsatwater_out[i]/UNIT_C>MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity && MD->DummyY[i + MD->NumEle]>=0)
			{
				if (MD->EleET[i][2]>0)
				{
					MD->EleET[i][2] = ((MD->EleET[i][2])/MD->Total_Unsatwater_out[i])*MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				if (MD->EleET[i][4]>0)
				{
					MD->EleET[i][4] = ((MD->EleET[i][4])/MD->Total_Unsatwater_out[i])*MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				if (MD->Recharge[i]>0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = ((MD->Recharge[i])/MD->Total_Unsatwater_out[i])*MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
						
			}
			else if (MD->DummyY[i + MD->NumEle]<0)
			{
				if (MD->EleET[i][2]>0)
				{
					MD->EleET[i][2] = 0;
				}
				if (MD->EleET[i][4]>0)
				{
					MD->EleET[i][4] = 0;
				}
				if (MD->Recharge[i]>0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = 0;
				}
			}
			//-----------Finish unsat water------------------//
			
			//-----------Saturated water------------------//
			if (MD->Total_Satwater_out[i]>0 && MD->Total_Satwater_out[i]/UNIT_C>MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity && MD->DummyY[i + 2*MD->NumEle]>=0)
			{
				for (j=0;j<3;j++)
				{
					inabr = MD->Ele[i].nabr[j] - 1; 
					
					if (MD->FluxSub[i][j]>0)
					{
						MD->FluxSub[i][j] = ((MD->FluxSub[i][j]/MD->Ele[i].area)/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C*MD->Ele[i].area;
						if (inabr>=0)
							MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSub[i][j];
					}	
				}				
				
				if (MD->EleET[i][3]>0)
				{
					MD->EleET[i][3] = ((MD->EleET[i][3])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				if (MD->EleET[i][5]>0)
				{
					MD->EleET[i][5] = ((MD->EleET[i][5])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				
				if (MD->EleViR[i]<0)
				{
					MD->EleViR[i] = ((MD->EleViR[i])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
					MD->Recharge[i] = MD->EleViR[i];
				}
				if (MD->Recharge[i]<0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = ((MD->Recharge[i])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}						
			}
			else if (MD->DummyY[i + 2*MD->NumEle]<0)
			{
				for (j=0;j<3;j++)
				{
					inabr = MD->Ele[i].nabr[j] - 1; 
					
					if (MD->FluxSub[i][j]>0)
					{
						MD->FluxSub[i][j] = 0;
						if (inabr>=0)
							MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSub[i][j];
					}	
				}				
				
				if (MD->EleET[i][3]>0)
				{
					MD->EleET[i][3] = 0;
				}
				if (MD->EleET[i][5]>0)
				{
					MD->EleET[i][5] = 0;
				}
				
				if (MD->EleViR[i]<0)
				{
					MD->EleViR[i] = 0;
					MD->Recharge[i] = MD->EleViR[i];
				}
				if (MD->Recharge[i]<0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = 0;
				}
			}
			//-----------Finish saturated water------------------//
		}		
	}
}

void SedAdvec(Model_Data MD, int i,int j,int inabr)
{
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	realtype		Niu = 0.000001004 //Unit: (m^2/s)
	if (inabr < 0)
	{
		Avg_Y_Surf = avgY(MD->FluxSurf[i][j], MD->DummyY[i], MD->DummyY[i]);
		Dif_Y_Surf = MD->DummyY[i];
		Grad_Y_Surf = Dif_Y_Surf/Distance;
		Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
	}
	else
	{
		Avg_Y_Surf = avgY(MD->FluxSurf[i][j], MD->DummyY[i], MD->DummyY[inabr]);
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[inabr+3*MD->NumEle]) : (MD->DummyY[i] + MD->DummyY[i + 3*MD->NumEle]) - (MD->DummyY[inabr] + MD->DummyY[inabr+3*MD->NumEle]);
		}
		else
		{
			Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->DummyY[i + MD->NumEle] - MD->DummyY[inabr+MD->NumEle]) : (MD->DummyY[i] + MD->DummyY[i + MD->NumEle]) - (MD->DummyY[inabr] + MD->DummyY[inabr+MD->NumEle]);
		}
		Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
	}
	CrossA = Avg_Y_Surf * MD->Ele[i].edge[j];
	Grad_Y_Surf = Dif_Y_Surf / Distance;
	if (Avg_Y_Surf>0)
	{
		// Tau0 = 1000*Cf*pow(fabs(MD->FluxSurf[i][j])/(CrossA*60*60*24),2);//shear stress, unit: kg/m*s^2
		Tau0 = 1000*9.8*((Avg_Y_Surf>0)?Avg_Y_Surf:0)*fabs(Grad_Y_Surf);
		U_star = pow(Tau0/1000,0.5);//U star, unit: m/s
		Renolds = U_star*MD->Ele[i].Diameter/Niu; //Renolds number, dimensionless				
		if (MD->SedMode == 0)
		{
			Critical_ShieldsStress = 0;
		}
		else
		{
			Critical_ShieldsStress = 0.5*(0.22*pow(Renolds,(-0.6))+0.06*exp(-17.77*pow(Renolds,(-0.6))));
		}									
		ShieldsStress = Tau0/((MD->Ele[i].RhoReg-1000)*9.8*MD->Ele[i].Diameter);
						
		if (ShieldsStress>Critical_ShieldsStress)
		{
			q_star = 3.97*pow(ShieldsStress-Critical_ShieldsStress,1.5);
			q_b=q_star * pow((MD->Ele[i].RhoReg/1000-1)*9.8*MD->Ele[i].Diameter,0.5)*MD->Ele[i].Diameter;//unit: m^2/min
			MD->BedloadFlux[i][j] = q_b*MD->Ele[i].edge[j]*60*MD->MSF;
			//MD->BedloadFlux[i][j]=0;
			if (MD->FluxSurf[i][j]<0)
			{
				MD->BedloadFlux[i][j]=-MD->BedloadFlux[i][j];
				//MD->BedloadFlux[i][j]=0; 
			}
		}
		else
		{						
			MD->BedloadFlux[i][j] = 0;
		}					
	}
	else
	{
		MD->BedloadFlux[i][j] = 0;
	}
	if (MD->BedloadFlux[i][j]>0)
	{
		MD->Total_soil_out[i] = MD->Total_soil_out[i] + MD->BedloadFlux[i][j]/MD->Ele[i].area;
	}
}

void Adjust_sed(Model_Data MD)
{
	int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	
	is_outlet = 0;
	
	for (i = 0; i < MD->NumEle; i++) 
	{
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			AquiferDepth = (MD->DummyY[i + 3*MD->NumEle] - MD->DummyY[i + 4*MD->NumEle]);
			if (AquiferDepth < 0)
				AquiferDepth = 0;
		}
		else
		{
			AquiferDepth = (MD->DummyY[i + MD->NumEle] - MD->DummyY[i + 2*MD->NumEle]);
			if (AquiferDepth < 0)
				AquiferDepth = 0;
		}
		
		for (q=0;q<MD->NumOutlet;q++)
		{
			if (MD->Outlet_location[q][0]-1==i)
			{
				is_outlet = 1;
			}
		}
		
		for (j = 0; j < 3; j++) 
		{
			inabr = MD->Ele[i].nabr[j] - 1;
			if (is_outlet==0 && MD->Total_soil_out[i]>0 && (MD->Total_soil_out[i]) > (AquiferDepth))
			{
				if (MD->BedloadFlux[i][j]>0)
				{							
					MD->BedloadFlux[i][j] = ((MD->BedloadFlux[i][j]/MD->Ele[i].area)/MD->Total_soil_out[i])*(AquiferDepth)*MD->Ele[i].area;
					if (inabr>=0)
					{
						MD->BedloadFlux[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->BedloadFlux[i][j];
					} 	
				}
						
				if (MD->SoilFluxSub[i][j]>0)
				{
					MD->SoilFluxSub[i][j] = ((MD->SoilFluxSub[i][j]/MD->Ele[i].area)/MD->Total_soil_out[i])*(AquiferDepth )*MD->Ele[i].area;
					if (inabr>=0)
					{
						MD->SoilFluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->SoilFluxSub[i][j];
					} 	
				}									
			}
		}
		is_outlet = 0;
	}
	for (q=0;q<MD->NumOutlet;q++)
	{
		out_i = MD->Outlet_location[q][0]-1;
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			AquiferDepth = MD->DummyY[out_i + 3*MD->NumEle]-MD->zmax_init[out_i];
		}
		else
		{
			AquiferDepth = MD->DummyY[out_i + MD->NumEle]-MD->zmax_init[out_i];
		}
		
		if (MD->Total_soil_out[out_i]>0 && MD->Total_soil_out[out_i]> AquiferDepth && AquiferDepth>0)
		{
			for (j = 0; j < 3; j++) 
			{
				if (MD->BedloadFlux[out_i][j]>0)
				{
					MD->BedloadFlux[out_i][j] = (MD->BedloadFlux[out_i][j]/MD->Ele[out_i].area)/(MD->Total_soil_out[out_i]) * (AquiferDepth)*MD->Ele[out_i].area;
				}
				if (MD->SoilFluxSub[out_i][j]>0)
				{
					MD->SoilFluxSub[out_i][j] = (MD->SoilFluxSub[out_i][j]/MD->Ele[out_i].area)/(MD->Total_soil_out[out_i]) * (AquiferDepth)*MD->Ele[out_i].area;
				}
			}
		}
		else if (AquiferDepth<=0)
		{
			for (j = 0; j < 3; j++) 
			{
				if (MD->BedloadFlux[out_i][j]>0)
				{
					MD->BedloadFlux[out_i][j] = 0;
				}
				if (MD->SoilFluxSub[out_i][j]>0)
				{
					MD->SoilFluxSub[out_i][j] = 0;
				}
			}
		}
	}
}

realtype* Calculate_DY(Model_Data MD,realtype *DY)
{
	int             i,j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	
	for (i=0;i<MD->NumEle;i++)
	{
		for (j = 0; j < 3; j++) 
		{			
			DY[i] = DY[i] - MD->FluxSurf[i][j] / MD->Ele[i].area - MD->Seepage[i][j] / MD->Ele[i].area;
			if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
			{
				DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->FluxSub[i][j] / MD->Ele[i].area;
				if (MD->Keyword_LE_PIHM == 1)
				{
					DY[i + 3 * MD->NumEle] = DY[i + 3 * MD->NumEle] - MD->SoilFluxSub[i][j]/MD->Ele[i].area - MD->BedloadFlux[i][j]/MD->Ele[i].area;
				}
			}
			else
			{
				if (MD->Keyword_LE_PIHM == 1)
				{
					DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->SoilFluxSub[i][j]/MD->Ele[i].area - MD->BedloadFlux[i][j]/MD->Ele[i].area;
				}
			}
		}
		DY[i] = DY[i] + MD->EleNetPrep[i] - MD->EleViR[i];
		DY[i] = DY[i] - MD->EleET[i][1];
		DY[i] = DY[i] / (UNIT_C);
		if (MD->Keyword_LE_PIHM_SUBHYDRO == 1)
		{
			DY[i + MD->NumEle] = DY[i + MD->NumEle] + MD->EleViR[i] - MD->Recharge[i];
			DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] + MD->Recharge[i];
			DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->EleET[i][2];
			DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->EleET[i][3];
			DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->EleET[i][4];
			DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->EleET[i][5];
			
			DY[i + MD->NumEle] =  DY[i + MD->NumEle] / (MD->Ele[i].Porosity * UNIT_C);
			DY[i + 2 * MD->NumEle] =  DY[i + 2 * MD->NumEle] / (MD->Ele[i].Porosity * UNIT_C);
			
			if (MD->Keyword_LE_PIHM == 1)
			{
				DY[i + 3 * MD->NumEle] = DY[i + 3 * MD->NumEle] + (MD->Uplift[i] + (MD->Ele[i].RhoBed/MD->Ele[i].RhoReg-1)*MD->Soilproduction[i]);
				DY[i + 4 * MD->NumEle] = DY[i + 4 * MD->NumEle] + (MD->Uplift[i] - MD->Soilproduction[i]);		
			}
		}
		else
		{
			if (MD->Keyword_LE_PIHM == 1)
			{
				DY[i + MD->NumEle] = DY[i + MD->NumEle] + (MD->Uplift[i] + (MD->Ele[i].RhoBed/MD->Ele[i].RhoReg-1)*MD->Soilproduction[i]);
				DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] + (MD->Uplift[i] - MD->Soilproduction[i]);		
			}
		}
	}
	return DY;
}
void SoilProduction(Model_Data MD, int i)
{
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev; 
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	
	if (MD->Keyword_LE_PIHM_SUBHYDRO==1)
	{
		AquiferDepth = MD->DummyY[i+3*MD->NumEle]-MD->DummyY[i+4*MD->NumEle];
		if (AquiferDepth<0)
			AquiferDepth = 0;
	}
	else
	{
		AquiferDepth = MD->DummyY[i+MD->NumEle]-MD->DummyY[i+2*MD->NumEle];
		if (AquiferDepth<0)
			AquiferDepth = 0;
	}
	Cos_slope = 1/sqrt(1 + pow(MD->Ele[i].Tan_slope,2));
	Sec_slope = sqrt(1 + pow(MD->Ele[i].Tan_slope,2));
	MD->Soilproduction[i] = MD->Ele[i].P0*exp(-MD->Ele[i].CoefP0*(AquiferDepth*Cos_slope))*Sec_slope*MD->MSF/525600;	
}

void BedUplift(Model_Data MD, int i)
{
	MD->Uplift[i] = MD->Ele[i].Uplift*MD->MSF/525600;
}

int
f(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
{
	int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	realtype        Delta, Gamma;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	realtype        ThetaRef, ThetaW;
	realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	realtype        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	realtype		Dif_Elev, Grad_Elev;
	realtype 		c1_square,c2_square,L,b;
	realtype       *Y, *DY;
	realtype		time_convert = 60*60*24/UNIT_C;
	realtype		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	realtype		Cf = 0.003;
	realtype		Cos_slope, Sec_slope;
	Model_Data      MD;
	Y = NV_DATA_S(CV_Y);
	DY = NV_DATA_S(CV_Ydot);
	MD = (Model_Data) DS;

	
	
	/*********initialize some parameters******************/
	DY = initiation(MD,DY,Y,t); //All DY[*] = 0
	/********************************************************/

	/*********Calculations among different Processes******************/
	for (i = 0; i < MD->NumEle; i++) 
	{		
		for (j = 0; j < 3; j++) 
		{
			inabr = MD->Ele[i].nabr[j] - 1;
			HydroSubLateral(MD,i,j,inabr);
			HydroSurfLateral(MD,i,j,inabr);
			if (MD->Keyword_LE_PIHM == 1)
				SedDiffu(MD,i,j,inabr);
		}
		Infiltration(MD,i);
		ET(MD,i,t);
		
		if (MD->Keyword_LE_PIHM == 1)
		{
			SoilProduction(MD, i);
			BedUplift(MD, i);
		}
	}
	Adjust_water(MD);
	
	if (MD->Keyword_LE_PIHM == 1)
	{
		/*********************Calculate sediment flux***********************************************/	
		for (i = 0; i < MD->NumEle; i++) 
		{
			for (j = 0; j < 3; j++)
			{	
				inabr = MD->Ele[i].nabr[j] - 1;
				SedAdvec(MD,i,j,inabr);		
			}
		}
		/****************************End of Bedload Sediment transport******************************/
		Adjust_sed(MD);
	}	
		
	//////////////////////calculate new state variables////////////////////
	DY = Calculate_DY(MD,DY);
			
	return 0;
}

realtype
Interpolation(TSD * Data, realtype t)
{
	int             i, success;
	realtype        result;
	i = Data->iCounter;
	success = 0;
	t = t / (UNIT_C);
	while (i < Data->length && t > Data->TS[i][0]) {
		i++;
	}
	if (i == 0) {
		/* t is smaller than the 1st node */
		result = Data->TS[i][1];
	} else if (i >= Data->length) {
		result = Data->TS[i - 1][1];
	} else {
		result = ((Data->TS[i][0] - t) * Data->TS[i - 1][1] + (t - Data->TS[i - 1][0]) * Data->TS[i][1]) / (Data->TS[i][0] - Data->TS[i - 1][0]);
		success = 1;
	}
	if (success == 0) {
		/*
    	printf("\nWarning:  Extrapolation is used ...\n");
    	*/
	}
	return result;
}