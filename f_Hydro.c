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

#define Niu 0.000001004 //Unit: (m^2/s)

int
f_Hydro(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
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
	Model_Data      MD;
	Y = NV_DATA_S(CV_Y);
	DY = NV_DATA_S(CV_Ydot);
	MD = (Model_Data) DS;

	
	
	//*********initialize some parameters******************
	for (i = 0; i < MD->NumEle; i++) {
		//MD->DummyY_Hydro[i] = (Y[i] >= 0) ? Y[i] : 0;
		MD->DummyY_Hydro[i] = Y[i];
		MD->DummyY_Hydro[i+MD->NumEle] = Y[i+MD->NumEle];
		MD->DummyY_Hydro[i+2*MD->NumEle] = Y[i+2*MD->NumEle];
		
		DY[i] = 0;
		DY[i+MD->NumEle] = 0;
		DY[i+2*MD->NumEle] = 0;
		// MD->Total_water_in[i] = MD->EleNetPrep[i];
		// MD->Total_water_out[i] = 0;
		MD->infil_mode[i] = 0;
		
		
		// #ifdef LE_PIHM_SED
			
		// #endif
		
	}
	is_outlet = 0;
	
	for (i = 0; i < MD->NumEle; i++) 
	{
		if ((MD->SurfMode == 2) && (i < MD->NumEle)) 
		{
			for (j = 0; j < 3; j++) 
			{
				// BHATT: MAJOR BUG DUMMYY OF NABR MAY BE NOT INITIALIZED
				MD->Ele[i].surfH[j] = (MD->Ele[i].nabr[j] > 0) ? (MD->Ele[MD->Ele[i].nabr[j] - 1].zmax + MD->DummyY_Hydro[MD->Ele[i].nabr[j] - 1]) : ((MD->Ele[i].BC[j] != 1) ? (MD->grdelev[i] + MD->DummyY_Hydro[i]) : Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
			}
			MD->Ele[i].dhBYdx = -1 * (MD->Ele[i].surfY[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfX[2] * (MD->Ele[i].surfY[1] - MD->Ele[i].surfY[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfY[0] - MD->Ele[i].surfY[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfY[2] - MD->Ele[i].surfY[1]));
			MD->Ele[i].dhBYdy = -1 * (MD->Ele[i].surfX[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfY[2] * (MD->Ele[i].surfX[1] - MD->Ele[i].surfX[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfX[0] - MD->Ele[i].surfX[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfX[2] - MD->Ele[i].surfX[1]));
		}
	}
	//********************************************************

	
	
	/* Lateral Flux Calculation between Triangular elements Follows  */
	for (i = 0; i < MD->NumEle; i++) 
	{		
		AquiferDepth = (MD->grdelev[i] - MD->bedelev[i]);
		#ifndef BEDROCK
			AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
		#endif
	
		if (AquiferDepth < MD->Ele[i].macD)
			MD->Ele[i].macD = AquiferDepth;
		
		for (j = 0; j < 3; j++) 
		{
			if (MD->Ele[i].nabr[j] > 0) 
			{
				/***************************************************************************/
				/*
				 * Subsurface Lateral Flux Calculation
				 * between Triangular elements Follows
				 */
				/***************************************************************************/
				
				inabr = MD->Ele[i].nabr[j] - 1;
				nabrAqDepth = (MD->grdelev[inabr] - MD->bedelev[inabr]);
				#ifndef BEDROCK
					nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
				#endif
				if (nabrAqDepth < MD->Ele[inabr].macD)
					MD->Ele[inabr].macD = nabrAqDepth;
				if (AquiferDepth>0.1&&nabrAqDepth>0.1)
				{
					Dif_Y_Sub = (MD->DummyY_Hydro[i + 2 * MD->NumEle] + MD->grdelev[i]-AquiferDepth) - (MD->DummyY_Hydro[inabr + 2 * MD->NumEle] + MD->grdelev[inabr]-nabrAqDepth);
					
					
					Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY_Hydro[i + 2 * MD->NumEle], MD->DummyY_Hydro[inabr + 2 * MD->NumEle]);
					Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
					Grad_Y_Sub = Dif_Y_Sub / Distance;
					/* take care of macropore effect */
					effK = effKH(MD->Ele[i].Macropore, MD->DummyY_Hydro[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
					effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY_Hydro[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
					/*
					 * It should be weighted average. However,
					* there is an ambiguity about distance used
					 */
					Avg_Ksat = 0.5 * (effK + effKnabr);
					/* groundwater flow modeled by Darcy's law */
					MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
				}
				else
				{
					MD->FluxSub[i][j] = 0;
				}
				
				/***************************************************************************/
				/*
				 * Surface Lateral Flux Calculation between
				 * Triangular elements Follows
				 */
				/***************************************************************************/
				Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->grdelev[i] - MD->grdelev[inabr]) : (MD->DummyY_Hydro[i] + MD->grdelev[i]) - (MD->DummyY_Hydro[inabr] + MD->grdelev[inabr]);
					Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY_Hydro[i], MD->DummyY_Hydro[inabr]);
					Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
					Grad_Y_Surf = Dif_Y_Surf / Distance;
					Avg_Sf = 0.5*(sqrt(pow(MD->Ele[i].dhBYdx, 2) + pow(MD->Ele[i].dhBYdy, 2))+sqrt(pow(MD->Ele[inabr].dhBYdx, 2) + pow(MD->Ele[inabr].dhBYdy, 2)));//?? Xuan Weighting needed
					//updated in 2.2
					Avg_Sf = (MD->SurfMode == 1) ? (Grad_Y_Surf > 0 ? Grad_Y_Surf : EPS / pow(10.0, 6)) : (Avg_Sf > EPS / pow(10.0, 6)) ? Avg_Sf : EPS / pow(10.0, 6);
					/* Weighting needed */
					Avg_Rough = 0.5 * (MD->Ele[i].Rough + MD->Ele[inabr].Rough);
					CrossA = Avg_Y_Surf * MD->Ele[i].edge[j];
					MD->FluxSurf[i][j] = CrossA * pow(Avg_Y_Surf, 2.0 / 3.0) * Grad_Y_Surf / (sqrt(fabs(Avg_Sf)) * Avg_Rough);
					//OverlandFlow(MD->FluxSurf, i, j, Avg_Y_Surf, Grad_Y_Surf, Avg_Sf, CrossA, Avg_Rough);
					if (MD->FluxSurf[i][j]>0)
					{
						// MD->Total_water_out[i] = MD->Total_water_out[i] + MD->FluxSurf[i][j]/MD->Ele[i].area;
					}
					else
					{
						// MD->Total_water_in[i] = MD->Total_water_in[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
					}
				
				// #ifdef LE_PIHM_SED
					
				// #endif
			}
			/************************************************/
			/* Boundary Condition Flux Calculations Follows */
			/************************************************/
			else 
			{
				/*
				 * No flow (natural) boundary condition is
				 * default
				 */
				if (MD->Ele[i].BC[j] == 0) 
				{
					MD->FluxSurf[i][j] = 0;
					MD->FluxSub[i][j] = 0;
				} 
				else if (MD->Ele[i].BC[j] == 1) 
				{	/* Note: ideally
									 * different boundary
									 * conditions need to be
									 * incorporated	for surf
									 * and subsurf
									 * respectively */
					/*
					 * Note: the formulation assumes only
					 * dirichlet TS right now
					 */
					MD->SoilFluxSurf[i][j] = 0;
					MD->SoilFluxSub[i][j] = 0;
					MD->FluxSurf[i][j] = 0;	/* Note the assumption
								 * here is no flow for
								 * surface */
					Dif_Y_Sub = (MD->DummyY_Hydro[i + 2 * MD->NumEle] + MD->bedelev[i]) - Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
					Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY_Hydro[i + 2 * MD->NumEle], (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->bedelev[i]));
					//Avg_Y_Sub = (MD->DummyY_Hydro[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->bedelev[i])) / 2;
					/*
					 * Minimum Distance from circumcenter
					 * to the edge of the triangle on
					 * which BDD. condition is defined
					 */
					Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
					effK = effKH(MD->Ele[i].Macropore, MD->DummyY_Hydro[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
					Avg_Ksat = effK;
					Grad_Y_Sub = Dif_Y_Sub / Distance;
					MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
				} 
				else if (MD->Ele[i].BC[j] == 3)
				{
						Dif_Y_Surf = (MD->DummyY_Hydro[i]+MD->grdelev[i] - MD->zmax_init[i]);
						Avg_Y_Surf = MD->DummyY_Hydro[i];
						Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
						Grad_Y_Surf = Dif_Y_Surf/Distance;
						Avg_Rough = MD->Ele[i].Rough;
						CrossA = MD->DummyY_Hydro[i]*MD->Ele[i].edge[j];        
						MD->FluxSurf[i][j] = sqrt((Grad_Y_Surf>0)?Grad_Y_Surf:0)*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough; 														
						MD->FluxSub[i][j] = 0;
						if (MD->FluxSurf[i][j]>0)
						{
							// MD->Total_water_out[i] = MD->Total_water_out[i] + MD->FluxSurf[i][j]/MD->Ele[i].area;
						}
						else
						{
							// MD->Total_water_in[i] = MD->Total_water_in[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
						}																					
				}
				else 
				{	/* Neumann BC (Note:
					* MD->Ele[i].BC[j] value
					* have to be = 2+(index of
					* neumann boundary TS) */
					MD->FluxSurf[i][j] = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
					MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC[j]) - 1], t);
					// MD->SoilFluxSurf[i][j] = 0;
					// MD->SoilFluxSub[i][j] = 0;					
				}
			}
		}
							
		// #ifdef LE_PIHM_SED
								
		// #endif
				
		/**************************************************************************************************/
		/*
		 * Evaporation Module: [2] is ET from OVLF/SUBF, [1] is
		 * Transpiration, [0] is ET loss from canopy
		 */
		/**************************************************************************************************/
		
		#ifdef LE_PIHM_HYDRO
			/* Physical Unit Dependent. Change this */
			Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn - 1], t);
			//G = Interpolation(&MD->TSD_G[MD->Ele[i].G - 1], t);
			G = 0.1 * Rn;
			T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp - 1], t);
			Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel - 1], t);
			RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity - 1], t);
			VP = 611.2 * exp(17.67 * T / (T + 243.5)) * RH;
			P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->grdelev[i]) / 293, 5.26);
			qv = 0.622 * VP / P;
			qv_sat = 0.622 * (VP / RH) / P;
			//P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->grdelev[i]) / 293, 5.26);
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
			if (AquiferDepth - MD->DummyY_Hydro[i + 2 * MD->NumEle] < MD->Ele[i].RzD) {
				elemSatn = 1.0;
			} else {
				elemSatn = ((MD->DummyY_Hydro[i + MD->NumEle] / (AquiferDepth - MD->DummyY_Hydro[i + 2 * MD->NumEle])) > 1) ? 1 : ((MD->DummyY_Hydro[i + MD->NumEle] / (AquiferDepth - MD->DummyY_Hydro[i + 2 * MD->NumEle])) < 0) ? 0 : 0.5 * (1 - cos(3.14 * (MD->DummyY_Hydro[i + MD->NumEle] / (AquiferDepth - MD->DummyY_Hydro[i + 2 * MD->NumEle]))));
			}
			ThetaRef = 0.7 * MD->Soil[(MD->Ele[i].soil - 1)].ThetaS;
			ThetaW = 1.05 * MD->Soil[(MD->Ele[i].soil - 1)].ThetaR;
			beta_s = (elemSatn * MD->Ele[i].Porosity + MD->Soil[(MD->Ele[i].soil - 1)].ThetaR - ThetaW) / (ThetaRef - ThetaW);
			beta_s = (beta_s < 0.0001) ? 0.0001 : (beta_s > 1 ? 1 : beta_s);
			MD->EleET[i][2] = MD->pcCal.Et2 * (1 - MD->Ele[i].VegFrac) * beta_s * ETp;
			MD->EleET[i][2] = MD->EleET[i][2] < 0 ? 0 : MD->EleET[i][2];
			MD->EleET[i][2] = (MD->DummyY_Hydro[i] < EPS / 100) ? elemSatn * MD->EleET[i][2] : MD->EleET[i][2];
			MD->EleET[i][2] = (MD->DummyY_Hydro[i] < EPS / 100) ? (elemSatn <= multF * EPS ? 0 : MD->EleET[i][2]) : MD->EleET[i][2];
			if (LAI > 0.0) {
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
				MD->EleET[i][1] = MD->pcCal.Et1 * MD->Ele[i].VegFrac * P_c * (1 - pow(((MD->EleIS[i] + MD->EleSnowCanopy[i] < 0) ? 0 : (MD->EleIS[i] + MD->EleSnowCanopy[i])) / (MD->EleISmax[i] + MD->EleISsnowmax[i]), 1.0 / 2.0)) * ETp;
				MD->EleET[i][1] = MD->EleET[i][1] < 0 ? 0 : MD->EleET[i][1];
				//? ? BHATT
					MD->EleET[i][1] = ((MD->DummyY_Hydro[i + 2 * MD->NumEle] < (AquiferDepth - MD->Ele[i].RzD)) && MD->DummyY_Hydro[i + MD->NumEle] <= 0) ? 0 : MD->EleET[i][1];
				//? ? BHATT
			} else {
				MD->EleET[i][1] = 0.0;
			}
		#else
			MD->EleET[i][2] = 0;
			MD->EleET[i][1] = 0;
		#endif
		
		#ifdef LE_PIHM
			#ifndef LE_PIHM_HYDRO
				MD->EleET[i][2] = 0;
				MD->EleET[i][1] = 0;
			#endif
		#endif
		
		/*
		 * Note: Assumption is OVL flow depth less than EPS/100 is
		 * immobile water
		 */		
		Deficit = AquiferDepth - MD->DummyY_Hydro[i + 2 * MD->NumEle];
		if (AquiferDepth>0.1)
		{
			if (MD->DummyY_Hydro[i+MD->NumEle]<Deficit&&Deficit>MD->Ele[i].infD)
			{
				MD->infil_mode[i] = 1;
				elemSatn = MD->DummyY_Hydro[i+MD->NumEle]/Deficit;
				elemSatn = (elemSatn < multF * EPS) ? multF * EPS : elemSatn;
				Avg_Y_Sub = (-(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha) < MINpsi) ? MINpsi : -(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha);
				
				TotalY_Ele = Avg_Y_Sub + MD->grdelev[i] - Deficit;
				Grad_Y_Sub = (MD->DummyY_Hydro[i] + MD->grdelev[i] - TotalY_Ele) / Deficit;
				Grad_Y_Sub = ((MD->DummyY_Hydro[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
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
				
				effK = (MD->Ele[i].Macropore == 1) ? ((MD->DummyY_Hydro[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].macD) ? effK : MD->Ele[i].KsatV * satKfunc) : MD->Ele[i].KsatV * satKfunc;
				MD->Recharge[i] = (elemSatn == 0.0) ? 0 : (Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->DummyY_Hydro[i + 2 * MD->NumEle] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow(-1 + pow(elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow(Deficit + MD->DummyY_Hydro[i + 2 * MD->NumEle], 2));
				MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->DummyY_Hydro[i + MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
				MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->DummyY_Hydro[i + 2 * MD->NumEle] <= 0) ? 0 : MD->Recharge[i];							
			}
			else
			{
				MD->infil_mode[i] = -1;
				Grad_Y_Sub = (MD->DummyY_Hydro[i] + MD->grdelev[i] - (MD->DummyY_Hydro[i + 2 * MD->NumEle] +  MD->grdelev[i] - AquiferDepth)) / AquiferDepth;
				
				Grad_Y_Sub = ((MD->DummyY_Hydro[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
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
	
///////////////////Ajust surface water flux////////////////////
	for (i = 0; i < MD->NumEle; i++)  
	{
		if ((MD->EleViR[i] - MD->EleNetPrep[i])/UNIT_C>MD->DummyY_Hydro[i])
		{
			MD->EleViR[i] = MD->EleNetPrep[i]+MD->DummyY_Hydro[i]*UNIT_C;
		}
		if (MD->infil_mode[i]==-1)
		{
			MD->Recharge[i] = MD->EleViR[i];
		}
		
		for (j = 0; j < 3; j++)
		{	
			inabr = MD->Ele[i].nabr[j] - 1;
			// if (MD->Total_water_out[i]>0&&(MD->Total_water_out[i])/UNIT_C > (MD->DummyY_Hydro[i]+(MD->EleNetPrep[i]-MD->EleViR[i])/UNIT_C))
			{
				if (MD->FluxSurf[i][j]>0)
				{
					// MD->FluxSurf[i][j] = ((MD->FluxSurf[i][j]/MD->Ele[i].area)/MD->Total_water_out[i])*(MD->DummyY_Hydro[i] * UNIT_C + MD->EleNetPrep[i] - MD->EleViR[i])*MD->Ele[i].area;
					if (inabr>=0)
					{
						MD->FluxSurf[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSurf[i][j];
					}
				}
			}
						
		}
		
	}
	
//////////////////////calculate new state variables////////////////////
	for (i = 0; i < MD->NumEle; i++) 
	{
		
		AquiferDepth = MD->grdelev[i] - MD->bedelev[i];
		Deficit = AquiferDepth - MD->DummyY_Hydro[i + 2 * MD->NumEle];
		
		for (j = 0; j < 3; j++) 
		{			
			DY[i] = DY[i] - MD->FluxSurf[i][j] / MD->Ele[i].area;
			DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->FluxSub[i][j] / MD->Ele[i].area;
		}
		
		DY[i] = DY[i] + MD->EleNetPrep[i] - MD->EleViR[i];
		if (Deficit>0)
		{
			DY[i + MD->NumEle] = DY[i + MD->NumEle] + MD->EleViR[i] - MD->Recharge[i];
		}
		DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] + MD->Recharge[i];
		
		if (MD->DummyY_Hydro[i] < EPS / 100)
		{
			if (Deficit>0)
			{
				DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->EleET[i][2];
			}
			else
			{
				DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->EleET[i][2];
			}
		}
		else
		{
			DY[i] = DY[i] - MD->EleET[i][2];
		}
		
		if (MD->DummyY_Hydro[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].RzD) 
		{
			DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->EleET[i][1];
		} 
		else 
		{
			DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->EleET[i][1];
		}		
				 
		DY[i] = DY[i] / (UNIT_C);
		DY[i + MD->NumEle] =  DY[i + MD->NumEle] / (MD->Ele[i].Porosity * UNIT_C);
		DY[i + 2 * MD->NumEle] =  DY[i + 2 * MD->NumEle] / (MD->Ele[i].Porosity * UNIT_C);   
		// if ((DY[i]<-1)||(DY[i + MD->NumEle])<-1||(DY[i + 2 * MD->NumEle]<-1))
		// {
			// printf("%d\n", i);
		// }
		
		// if ((DY[i]>1)||(DY[i + MD->NumEle])>1||(DY[i + 2 * MD->NumEle]>1))
		// {
			// printf("%d\n", i);
		// }
	}
	return 0;
}