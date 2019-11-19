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
f_LE(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
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
		//MD->SurfDepth[i] = (Y[i] >= 0) ? Y[i] : 0;
		
		#ifdef LE_PIHM
			MD->DummyY_LE[i] = Y[i];
			MD->DummyY_LE[i+MD->NumEle] = Y[i+MD->NumEle];
		#endif
		
		// #ifdef LE_PIHM_SED
			
		// #endif
		
		
		// MD->Total_water_in[i] = MD->EleNetPrep[i];
		// MD->Total_water_out[i] = 0;
		MD->infil_mode[i] = 0;
		#ifdef LE_PIHM
			MD->Total_soil_in[i] = 0;
			MD->Total_soil_out[i] = 0;
			DY[i] = 0;
			DY[i+MD->NumEle] = 0;
		#endif
		
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
				#ifdef LE_PIHM
					MD->Ele[i].surfH[j] = (MD->Ele[i].nabr[j] > 0) ? (MD->SurfDepth[MD->Ele[i].nabr[j] - 1] + MD->DummyY_LE[MD->Ele[i].nabr[j] - 1]) : ((MD->Ele[i].BC[j] != 1) ? (MD->DummyY_LE[i] + MD->SurfDepth[i]) : Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
				#endif
			}
			MD->Ele[i].dhBYdx = -1 * (MD->Ele[i].surfY[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfX[2] * (MD->Ele[i].surfY[1] - MD->Ele[i].surfY[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfY[0] - MD->Ele[i].surfY[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfY[2] - MD->Ele[i].surfY[1]));
			MD->Ele[i].dhBYdy = -1 * (MD->Ele[i].surfX[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfY[2] * (MD->Ele[i].surfX[1] - MD->Ele[i].surfX[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfX[0] - MD->Ele[i].surfX[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfX[2] - MD->Ele[i].surfX[1]));
		}
	}
	//********************************************************

	
	
	/* Lateral Flux Calculation between Triangular elements Follows  */
	for (i = 0; i < MD->NumEle; i++) 
	{		
		#ifdef LE_PIHM
			AquiferDepth = (MD->DummyY_LE[i] - MD->DummyY_LE[i+MD->NumEle]);
			#ifndef BEDROCK
				AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
			#endif
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
				#ifdef LE_PIHM
					nabrAqDepth = (MD->DummyY_LE[inabr] - MD->DummyY_LE[inabr+MD->NumEle]);
					#ifndef BEDROCK
						nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
					#endif
				#endif
				if (nabrAqDepth < MD->Ele[inabr].macD)
					MD->Ele[inabr].macD = nabrAqDepth;
				if (AquiferDepth>0.1&&nabrAqDepth>0.1)
				{
					#ifdef LE_PIHM
						Dif_Y_Sub = (MD->GW[i] + MD->DummyY_LE[i]-AquiferDepth) - (MD->GW[inabr] + MD->DummyY_LE[inabr]-nabrAqDepth);
					#endif
					
					Avg_Y_Sub = avgY(Dif_Y_Sub, MD->GW[i], MD->GW[inabr]);
					Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
					Grad_Y_Sub = Dif_Y_Sub / Distance;
					/* take care of macropore effect */
					effK = effKH(MD->Ele[i].Macropore, MD->GW[i], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
					effKnabr = effKH(MD->Ele[inabr].Macropore, MD->GW[inabr], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
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
				#ifdef LE_PIHM
					Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->DummyY_LE[i] - MD->DummyY_LE[inabr]) : (MD->SurfDepth[i] + MD->DummyY_LE[i]) - (MD->SurfDepth[inabr] + MD->DummyY_LE[inabr]);
				#endif
					Avg_Y_Surf = avgY(Dif_Y_Surf, MD->SurfDepth[i], MD->SurfDepth[inabr]);
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
				#ifdef LE_PIHM
					/*********************Regolith lateral flux movement************************/				
					Dif_Elev = (MD->DummyY_LE[i] - MD->DummyY_LE[inabr]);
					Grad_Elev = Dif_Elev/Distance;
					MD->SoilFluxSub[i][j] = 0.5*(MD->Ele[i].DReg + MD->Ele[inabr].DReg)*Grad_Elev*MD->Ele[i].edge[j]*MD->MSF/525600;
						
					if (MD->SoilFluxSub[i][j]>0)
					{
						MD->Total_soil_out[i] = MD->Total_soil_out[i] + MD->SoilFluxSub[i][j]/MD->Ele[i].area;
					}
					/***************************************************************************/
				#endif
				
				
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
					#ifdef LE_PIHM
						MD->SoilFluxSurf[i][j] = 0;
						MD->SoilFluxSub[i][j] = 0;
						MD->BedloadFlux[i][j] = 0;
					#endif
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
					Dif_Y_Sub = (MD->GW[i] + MD->Ele[i].zmin) - Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
					Avg_Y_Sub = avgY(Dif_Y_Sub, MD->GW[i], (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin));
					//Avg_Y_Sub = (MD->GW[i] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
					/*
					 * Minimum Distance from circumcenter
					 * to the edge of the triangle on
					 * which BDD. condition is defined
					 */
					Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
					effK = effKH(MD->Ele[i].Macropore, MD->GW[i], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
					Avg_Ksat = effK;
					Grad_Y_Sub = Dif_Y_Sub / Distance;
					MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
				} 
				else if (MD->Ele[i].BC[j] == 3)
				{
						Dif_Y_Surf = (MD->SurfDepth[i]+MD->DummyY_LE[i] - MD->zmax_init[i]);
						Avg_Y_Surf = MD->SurfDepth[i];
						Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
						Grad_Y_Surf = Dif_Y_Surf/Distance;
						Avg_Rough = MD->Ele[i].Rough;
						CrossA = MD->SurfDepth[i]*MD->Ele[i].edge[j];        
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
						
							
							#ifdef LE_PIHM
								/*********************Regolith lateral flux movement************************/				
								Dif_Elev = (MD->DummyY_LE[i] - MD->zmax_init[i]);
								Grad_Elev = Dif_Elev/Distance;
								MD->SoilFluxSub[i][j] = (MD->Ele[i].DReg)*((Grad_Elev>0)?Grad_Elev:0)*MD->Ele[i].edge[j]*MD->MSF/525600;
								if (MD->SoilFluxSub[i][j]>0)
								{
									MD->Total_soil_out[i] = MD->Total_soil_out[i] + MD->SoilFluxSub[i][j]/MD->Ele[i].area;
								}
								/****************************************************************************************/
							#endif
							
							// #ifdef LE_PIHM_SED
								
							// #endif															
					}
				else {	/* Neumann BC (Note:
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
		#ifdef LE_PIHM
			if (AquiferDepth < 0)
			AquiferDepth = 0;
			MD->Soilproduction[i] = MD->Ele[i].P0*exp(-MD->Ele[i].CoefP0*(AquiferDepth))*MD->MSF/525600;
			MD->Uplift[i] = MD->Ele[i].Uplift*MD->MSF/525600;;
		#endif
		/*
		 * Note: Assumption is OVL flow depth less than EPS/100 is
		 * immobile water
		 */		
		Deficit = AquiferDepth - MD->GW[i];
		if (AquiferDepth>0.1)
		{
			if (MD->UGW[i]<Deficit&&Deficit>MD->Ele[i].infD)
			{
				MD->infil_mode[i] = 1;
				elemSatn = MD->UGW[i]/Deficit;
				elemSatn = (elemSatn < multF * EPS) ? multF * EPS : elemSatn;
				Avg_Y_Sub = (-(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha) < MINpsi) ? MINpsi : -(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha);
				
				#ifdef LE_PIHM
					TotalY_Ele = Avg_Y_Sub + MD->DummyY_LE[i] - Deficit;
					Grad_Y_Sub = (MD->SurfDepth[i] + MD->DummyY_LE[i] - TotalY_Ele) / Deficit;
				#endif
				Grad_Y_Sub = ((MD->SurfDepth[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
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
				
				effK = (MD->Ele[i].Macropore == 1) ? ((MD->GW[i] > AquiferDepth - MD->Ele[i].macD) ? effK : MD->Ele[i].KsatV * satKfunc) : MD->Ele[i].KsatV * satKfunc;
				MD->Recharge[i] = (elemSatn == 0.0) ? 0 : (Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->GW[i] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow(-1 + pow(elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow(Deficit + MD->GW[i], 2));
				MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->UGW[i] <= 0) ? 0 : MD->Recharge[i];
				MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->GW[i] <= 0) ? 0 : MD->Recharge[i];							
			}
			else
			{
				MD->infil_mode[i] = -1;				
				#ifdef LE_PIHM
					Grad_Y_Sub = (MD->SurfDepth[i] + MD->DummyY_LE[i] - (MD->GW[i] + MD->DummyY_LE[i]-AquiferDepth)) / AquiferDepth;
				#endif
				Grad_Y_Sub = ((MD->SurfDepth[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
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
		
///////////////////Ajust surface water flux ////////////////////
	for (i = 0; i < MD->NumEle; i++)  
	{
		if ((MD->EleViR[i] - MD->EleNetPrep[i])/UNIT_C>MD->SurfDepth[i])
		{
			MD->EleViR[i] = MD->EleNetPrep[i]+MD->SurfDepth[i]*UNIT_C;
		}
		if (MD->infil_mode[i]==-1)
		{
			MD->Recharge[i] = MD->EleViR[i];
		}
		
		for (j = 0; j < 3; j++)
		{	
			inabr = MD->Ele[i].nabr[j] - 1;
			// if (MD->Total_water_out[i]>0&&(MD->Total_water_out[i])/UNIT_C > (MD->SurfDepth[i]+(MD->EleNetPrep[i]-MD->EleViR[i])/UNIT_C))
			{
				if (MD->FluxSurf[i][j]>0)
				{
					// MD->FluxSurf[i][j] = ((MD->FluxSurf[i][j]/MD->Ele[i].area)/MD->Total_water_out[i])*(MD->SurfDepth[i] * UNIT_C + MD->EleNetPrep[i] - MD->EleViR[i])*MD->Ele[i].area;
					if (inabr>=0)
					{
						MD->FluxSurf[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSurf[i][j];
					}
				}
			}	
		}
		
	}
	
	#ifdef LE_PIHM
	/*********************Calculate sediment flux***********************************************/	
	for (i = 0; i < MD->NumEle; i++) 
	{
		for (j = 0; j < 3; j++)
		{
			inabr = MD->Ele[i].nabr[j] - 1;			
			/****************************************************************************************/
			/* Bedload Sediment transport
			*****************************************************************************************/
			if (MD->Ele[i].nabr[j] <= 0)
			{
				Avg_Y_Surf = avgY(MD->FluxSurf[i][j], MD->SurfDepth[i], MD->SurfDepth[i]);
			}else
			{
				Avg_Y_Surf = avgY(MD->FluxSurf[i][j], MD->SurfDepth[i], MD->SurfDepth[inabr]);
			}
			CrossA = Avg_Y_Surf * MD->Ele[i].edge[j];
			if (Avg_Y_Surf>0)
			{
				Tau0 = 1000*Cf*pow(fabs(MD->FluxSurf[i][j])/(CrossA*60*60*24),2);//shear stress, unit: kg/m*s^2
				//Tau0 = 1000*9.8*((Avg_Y_Surf>0)?Avg_Y_Surf:0)*fabs(Grad_Y_Surf);
				U_star = pow(Tau0/1000,0.5);//U star, unit: m/s
				Renolds = U_star*MD->Ele[i].Diameter/Niu; //Renolds number, dimensionless
					
				Critical_ShieldsStress = 0.5*(0.22*pow(Renolds,(-0.6))+0.06*exp(-17.77*pow(Renolds,(-0.6))));
					
				ShieldsStress = Tau0/((MD->Ele[i].RhoReg-1000)*9.8*MD->Ele[i].Diameter);
					
				if (ShieldsStress>0)
				{
					q_star = 8*pow(ShieldsStress-0,1.5);
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
	}
	 
/****************************End of Bedload Sediment transport**************************************************************/

#ifdef BEDROCK
/*********************Adjust sediment flux***********************************************/	
	for (i = 0; i < MD->NumEle; i++) 
	{
		AquiferDepth = (MD->DummyY_LE[i] - MD->DummyY_LE[i + MD->NumEle]);
		
		if (AquiferDepth < 0)
			AquiferDepth = 0;
			
		for (q=0;q<MD->NumOutlet;q++)
		{
			if (MD->Outlet_location[q][0]-1==i)
			{
				is_outlet = 1;
			}
		}
		
		for (j = 0; j < 3; j++) 
		{
			if (is_outlet==0 && MD->Total_soil_out[i]>0 && (MD->Total_soil_out[i]) > (AquiferDepth + (MD->Ele[i].RhoBed/MD->Ele[i].RhoReg-1)*MD->Soilproduction[i]))
			{
				if (MD->BedloadFlux[i][j]>0)
				{
					inabr = MD->Ele[i].nabr[j] - 1;
					MD->BedloadFlux[i][j] = ((MD->BedloadFlux[i][j]/MD->Ele[i].area)/MD->Total_soil_out[i])*(AquiferDepth + (MD->Ele[i].RhoBed/MD->Ele[i].RhoReg-1)*MD->Soilproduction[i])*MD->Ele[i].area;
					if (MD->SoilFluxSub[i][j]>0)
					{
						MD->SoilFluxSub[i][j] = ((MD->SoilFluxSub[i][j]/MD->Ele[i].area)/MD->Total_soil_out[i])*(AquiferDepth + (MD->Ele[i].RhoBed/MD->Ele[i].RhoReg-1)*MD->Soilproduction[i])*MD->Ele[i].area;
					}
					if (inabr>=0)
					{
						MD->BedloadFlux[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->BedloadFlux[i][j];
						MD->SoilFluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->SoilFluxSub[i][j];
					} 					
				}
			}
		}
		is_outlet = 0;
	}
#endif
	/******************************************************************************************/		
	
	for (q=0;q<MD->NumOutlet;q++)
	{
		out_i = MD->Outlet_location[q][0]-1;
		
		for (j = 0; j < 3; j++) 
		{
			if (MD->BedloadFlux[out_i][j]<0)
				MD->Total_soil_in[out_i] = MD->Total_soil_in[out_i] - MD->BedloadFlux[out_i][j]/MD->Ele[out_i].area;
			if (MD->SoilFluxSub[out_i][j]<0)
				MD->Total_soil_in[out_i] = MD->Total_soil_in[out_i] - MD->SoilFluxSub[out_i][j]/MD->Ele[out_i].area;
		}
		MD->Total_soil_in[out_i] = MD->Total_soil_in[out_i] + (MD->Ele[out_i].RhoBed/MD->Ele[out_i].RhoReg-1)*MD->Soilproduction[out_i];
		
		if (MD->Total_soil_out[out_i]>0 && MD->Total_soil_out[out_i]-MD->Total_soil_in[out_i] > MD->DummyY[out_i+3*MD->NumEle]-MD->zmax_init[out_i])
		{
			for (j = 0; j < 3; j++) 
			{
				if (MD->BedloadFlux[out_i][j]>0)
				{
					MD->BedloadFlux[out_i][j] = (MD->BedloadFlux[out_i][j]/MD->Ele[out_i].area)/(MD->Total_soil_out[out_i]) * (MD->DummyY[out_i + 3 * MD->NumEle]-MD->zmax_init[out_i] + MD->Total_soil_in[out_i])*MD->Ele[out_i].area;
				}
				if (MD->SoilFluxSub[out_i][j]>0)
				{
					MD->SoilFluxSub[out_i][j] = (MD->SoilFluxSub[out_i][j]/MD->Ele[out_i].area)/(MD->Total_soil_out[out_i]) * (MD->DummyY[out_i + 3 * MD->NumEle]-MD->zmax_init[out_i] + MD->Total_soil_in[out_i])*MD->Ele[out_i].area;
				}
			}
		}
	}	
	#endif		
		
	//////////////////////calculate new state variables////////////////////
	for (i = 0; i < MD->NumEle; i++) 
	{			
		for (j = 0; j < 3; j++) 
		{			
			#ifdef LE_PIHM
				DY[i] = DY[i] - MD->SoilFluxSub[i][j]/MD->Ele[i].area - MD->BedloadFlux[i][j]/MD->Ele[i].area;
			#endif
		}
		
		#ifdef LE_PIHM
			DY[i] = DY[i] + (MD->Uplift[i] + (MD->Ele[i].RhoBed/MD->Ele[i].RhoReg-1)*MD->Soilproduction[i]);
			DY[i+MD->NumEle] = DY[i+MD->NumEle] + (MD->Uplift[i] - MD->Soilproduction[i]);		
		#endif		      
	}
	return 0;
}