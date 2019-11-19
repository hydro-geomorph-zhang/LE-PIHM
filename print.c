/*******************************************************************************
 * File        : print.c	                                               *
 * Version     : Nov, 2007 (2.0)                                               *
 * Function    : print out model results output files                          *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) This file is downgraded from Version 1.0, as no ancillary results are    *
 *    being output			                                       *
 * b) Only state variables and flux to/in/accross river and its bed are being  *
 *    output							               *
 * c) Addition of Average Function to output average variables at regular time *
 *    intervals								       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundials_types.h"
#include "pihm.h"
#include "cvode.h"
#include "cvode_dense.h"
/* Temporal average of State vectors */
void
avgResults_NV(FILE * fpin, realtype * tmpVarCal, N_Vector tmpNV, int tmpIntv, int tmpNumObj, realtype tmpt, int tmpInitObj, Model_Data tmpDS)
{
	int             j;
	int             TmpIntv;
	//? ? TmpIntv = tmpIntv / 60;
	//USE THIS IF RUNNING AT 60 MIN STEP
	                TmpIntv = tmpIntv;
	///1;
	//? ?
		for (j = 0; j < tmpNumObj; j++) {
		if (tmpDS->Keyword_MEAN == 1)
		{
			tmpVarCal[j] = tmpVarCal[j] + NV_Ith_S(tmpNV, j + tmpInitObj);
		}
		if (tmpDS->Keyword_REALTIME == 1)
		{
			tmpVarCal[j] = NV_Ith_S(tmpNV, j + tmpInitObj);
		}
	}
	if (((int) tmpt % tmpIntv) == 0) {
		fprintf(fpin, "%lf\t", tmpt);
		for (j = 0; j < tmpNumObj; j++) {
		
		if (tmpDS->Keyword_MEAN == 1)
		{
			fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
		}
		if (tmpDS->Keyword_REALTIME == 1)
		{
			fprintf(fpin, "%lf\t", tmpVarCal[j]);
		}
			
		tmpVarCal[j] = 0;
		}
		fprintf(fpin, "\n");
		fflush(fpin);
	}
}

void
avgResults_NV_binary(FILE * fpin, realtype * tmpVarCal, N_Vector tmpNV, int tmpIntv, int tmpNumObj, realtype tmpt, int tmpInitObj, Model_Data tmpDS)
{
	int             j;
	int             TmpIntv;
	//? ? TmpIntv = tmpIntv / 60;
	//USE THIS IF RUNNING AT 60 MIN STEP
	                TmpIntv = tmpIntv;
	///1;
	//? ?
		for (j = 0; j < tmpNumObj; j++) {
		if (tmpDS->Keyword_MEAN == 1)
		{
			tmpVarCal[j] = tmpVarCal[j] + NV_Ith_S(tmpNV, j + tmpInitObj);
		}
		if (tmpDS->Keyword_REALTIME == 1)
		{
			tmpVarCal[j] = NV_Ith_S(tmpNV, j + tmpInitObj);
		}
	}
	if (((int) tmpt % tmpIntv) == 0) {
		int a;
		realtype temp_binary;
		a = fwrite(&tmpt,sizeof(realtype),1,fpin);
		
		for (j = 0; j < tmpNumObj; j++) 
		{		
			if (tmpDS->Keyword_MEAN == 1)
			{
				temp_binary = (realtype)  tmpVarCal[j] / TmpIntv;
				a = fwrite(&temp_binary,sizeof(realtype),1,fpin);
			}
			if (tmpDS->Keyword_REALTIME == 1)
			{
				temp_binary = (realtype) tmpVarCal[j];
				a = fwrite(&temp_binary,sizeof(realtype),1,fpin);
			}
				
			tmpVarCal[j] = 0;
		}
		fflush(fpin);
	}
}

/* Temporal average of Derived states */
void
avgResults_MD(FILE * fpin, realtype * tmpVarCal, Model_Data tmpDS, int tmpIntv, int tmpNumObj, realtype tmpt, int tmpFC)
{
	int             j;
	int             TmpIntv;
	//? ? TmpIntv = tmpIntv / 60;
	//USE THIS IF RUNNING AT 60 MIN STEP
	                TmpIntv = tmpIntv;
	//1;
	//? ?
		switch (tmpFC) {
	case 3 :
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][tmpFC - 3];
		}
		break;
	case 4 :
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][tmpFC - 3];
		}
		break;
	case 5:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][tmpFC - 3];
		}
		break;
	case 6:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleIS[j];
		}
		break;
	case 7:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + (tmpDS->EleSnowCanopy[j] + tmpDS->EleSnowGrnd[j]);
		}
		break;
	case 8:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->Recharge[j];
		}
		break;
	case 9:
	case 10:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleViR[j];
		}
		break;
	case 11:
			
	case 12:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->ElePrep[j]);
                }
            break;
	case 13:
			
	case 14:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][0];
			}
			break;
	case 15:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][1];
			}
			break;
	case 16:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][2];
			}
			break;
	case 17:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][0];
		}
		break;
	case 18:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][1];
		}
		break;
	case 19:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][2];
		}
		break;
	case 20:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->EleNetPrep[j]);
                }
            break;
	case 23:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->SoilFluxSub[j][0];
			}
			break;
			
	case 24:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->SoilFluxSub[j][1];
			}
			break;
	case 25:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->SoilFluxSub[j][2];
			}
			break;
	case 26:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->Soilproduction[j]);
                }
            break;
	case 27:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->Uplift[j]);
                }
            break;
	case 28:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->BedloadFlux[j][0];
			}
			break;
	case 29:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->BedloadFlux[j][1];
			}
			break;
	case 30:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->BedloadFlux[j][2];
			}
			break;
	default:
		break;
	}
	if (((int) tmpt % tmpIntv) == 0) {
		fprintf(fpin, "%lf\t", tmpt);
		for (j = 0; j < tmpNumObj; j++) {
			fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
			tmpVarCal[j] = 0;
		}
		fprintf(fpin, "\n");
		fflush(fpin);
	}
}

void
avgResults_MD_binary(FILE * fpin, realtype * tmpVarCal, Model_Data tmpDS, int tmpIntv, int tmpNumObj, realtype tmpt, int tmpFC)
{
	int             j;
	int             TmpIntv;
	//? ? TmpIntv = tmpIntv / 60;
	//USE THIS IF RUNNING AT 60 MIN STEP
	                TmpIntv = tmpIntv;
	//1;
	//? ?
		switch (tmpFC) {
	case 3 :
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][tmpFC - 3];
		}
		break;
	case 4 :
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][tmpFC - 3];
		}
		break;
	case 5:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][tmpFC - 3];
		}
		break;
	case 6:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleIS[j];
		}
		break;
	case 7:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + (tmpDS->EleSnowCanopy[j] + tmpDS->EleSnowGrnd[j]);
		}
		break;
	case 8:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->Recharge[j];
		}
		break;
	case 9:
	case 10:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleViR[j];
		}
		break;
	case 11:
			
	case 12:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->ElePrep[j]);
                }
            break;
	case 13:
			
	case 14:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][0];
			}
			break;
	case 15:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][1];
			}
			break;
	case 16:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][2];
			}
			break;
	case 17:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][0];
		}
		break;
	case 18:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][1];
		}
		break;
	case 19:
		for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][2];
		}
		break;
	case 20:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->EleNetPrep[j]);
                }
            break;
	case 23:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->SoilFluxSub[j][0];
			}
			break;
			
	case 24:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->SoilFluxSub[j][1];
			}
			break;
	case 25:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->SoilFluxSub[j][2];
			}
			break;
	case 26:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->Soilproduction[j]);
                }
            break;
	case 27:
			for(j=0;j<tmpNumObj;j++)
                {            
					tmpVarCal[j]=tmpVarCal[j]+(tmpDS->Uplift[j]);
                }
            break;
	case 28:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->BedloadFlux[j][0];
			}
			break;
	case 29:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->BedloadFlux[j][1];
			}
			break;
	case 30:
			for (j = 0; j < tmpNumObj; j++) {
			tmpVarCal[j] = tmpVarCal[j] + tmpDS->BedloadFlux[j][2];
			}
			break;
	default:
		break;
	}
	if (((int) tmpt % tmpIntv) == 0) 
	{
		int a;
		realtype	temp_binary;
		a = fwrite(&tmpt,sizeof(realtype),1,fpin);
		for (j = 0; j < tmpNumObj; j++) 
		{
			temp_binary = (realtype) tmpVarCal[j] / TmpIntv;
			a = fwrite(&temp_binary,sizeof(realtype),1,fpin);
			tmpVarCal[j] = 0;
		}
		fflush(fpin);
	}
}


/* print individual states */
void
PrintData(FILE ** outp, Control_Data * cD, Model_Data DS, N_Vector CV_Y, realtype t)
{
	FILE *Oinit_file;
	int i,k;
	if (DS->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		if (cD->gwDInt >0) {
			avgResults_NV(outp[0], DS->PrintVar[0], CV_Y, cD->gwDInt, DS->NumEle, t, 2 * DS->NumEle,DS);
		}
		if (cD->usDInt >0) {
			avgResults_NV(outp[7], DS->PrintVar[7], CV_Y, cD->usDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
		}
	}
	if (cD->surfDInt >0) {
		avgResults_NV(outp[1], DS->PrintVar[1], CV_Y, cD->surfDInt, DS->NumEle, t, 0 * DS->NumEle,DS);
	}
	for (k = 0; k < 3; k++) {
		if (cD->etInt >0) {
			avgResults_MD(outp[2 + k], DS->PrintVar[2 + k], DS, cD->etInt, DS->NumEle, t, 3 + k);
		}
	}
	if (cD->IsDInt >0) {
		avgResults_MD(outp[5], DS->PrintVar[5], DS, cD->IsDInt, DS->NumEle, t, 6);
	}
	if (cD->snowDInt >0) {
		avgResults_MD(outp[6], DS->PrintVar[6], DS, cD->snowDInt, DS->NumEle, t, 7);
	}
	
	if (cD->RechInt >0) {
		avgResults_MD(outp[8], DS->PrintVar[8], DS, cD->RechInt, DS->NumEle, t, 8);
		avgResults_MD(outp[9], DS->PrintVar[9], DS, cD->RechInt, DS->NumEle, t, 10);
	}
	if (cD->surfDInt >0) {
		avgResults_MD(outp[11], DS->PrintVar[11], DS, cD->surfDInt, DS->NumEle, t, 14);
		avgResults_MD(outp[12], DS->PrintVar[12], DS, cD->surfDInt, DS->NumEle, t, 15);
		avgResults_MD(outp[13], DS->PrintVar[13], DS, cD->surfDInt, DS->NumEle, t, 16);
		avgResults_MD(outp[14], DS->PrintVar[14], DS, cD->surfDInt, DS->NumEle, t, 17);
		avgResults_MD(outp[15], DS->PrintVar[15], DS, cD->surfDInt, DS->NumEle, t, 18);
		avgResults_MD(outp[16], DS->PrintVar[16], DS, cD->surfDInt, DS->NumEle, t, 19);
	}	
	if (cD->prepDInt >0)
		{
			avgResults_MD(outp[10],DS->PrintVar[10],DS,cD->prepDInt,DS->NumEle,t,12);
			avgResults_MD(outp[17],DS->PrintVar[17],DS,cD->prepDInt,DS->NumEle,t,20);
		}
}

void
PrintData_binary(FILE ** outp,FILE ** outp_binary, Control_Data * cD, Model_Data DS, N_Vector CV_Y, realtype t)
{
	FILE *Oinit_file;
	int i,k;
	if (DS->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		if (cD->gwDInt >0) {
			avgResults_NV_binary(outp_binary[0], DS->PrintVar_binary[0], CV_Y, cD->gwDInt, DS->NumEle, t, 2 * DS->NumEle,DS);
		}
		if (cD->usDInt >0) {
			avgResults_NV_binary(outp_binary[7], DS->PrintVar_binary[7], CV_Y, cD->usDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
		}
		if (cD->gwDInt >0) {
			avgResults_NV(outp[0], DS->PrintVar[0], CV_Y, cD->gwDInt, DS->NumEle, t, 2 * DS->NumEle,DS);
		}
		if (cD->usDInt >0) {
			avgResults_NV(outp[7], DS->PrintVar[7], CV_Y, cD->usDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
		}
	}
	if (cD->surfDInt >0) {
		avgResults_NV_binary(outp_binary[1], DS->PrintVar_binary[1], CV_Y, cD->surfDInt, DS->NumEle, t, 0 * DS->NumEle,DS);
		avgResults_NV(outp[1], DS->PrintVar[1], CV_Y, cD->surfDInt, DS->NumEle, t, 0 * DS->NumEle,DS);
	}
	for (k = 0; k < 3; k++) {
		if (cD->etInt >0) {
			avgResults_MD_binary(outp_binary[2 + k], DS->PrintVar_binary[2 + k], DS, cD->etInt, DS->NumEle, t, 3 + k);
		}
	}
	if (cD->IsDInt >0) {
		avgResults_MD_binary(outp_binary[5], DS->PrintVar_binary[5], DS, cD->IsDInt, DS->NumEle, t, 6);
	}
	if (cD->snowDInt >0) {
		avgResults_MD_binary(outp_binary[6], DS->PrintVar_binary[6], DS, cD->snowDInt, DS->NumEle, t, 7);
	}
	
	if (cD->RechInt >0) {
		avgResults_MD_binary(outp_binary[8], DS->PrintVar_binary[8], DS, cD->RechInt, DS->NumEle, t, 8);
		avgResults_MD_binary(outp_binary[9], DS->PrintVar_binary[9], DS, cD->RechInt, DS->NumEle, t, 10);
	}
	if (cD->surfDInt >0) {
		avgResults_MD_binary(outp_binary[11], DS->PrintVar_binary[11], DS, cD->surfDInt, DS->NumEle, t, 14);
		avgResults_MD_binary(outp_binary[12], DS->PrintVar_binary[12], DS, cD->surfDInt, DS->NumEle, t, 15);
		avgResults_MD_binary(outp_binary[13], DS->PrintVar_binary[13], DS, cD->surfDInt, DS->NumEle, t, 16);
		avgResults_MD_binary(outp_binary[14], DS->PrintVar_binary[14], DS, cD->surfDInt, DS->NumEle, t, 17);
		avgResults_MD_binary(outp_binary[15], DS->PrintVar_binary[15], DS, cD->surfDInt, DS->NumEle, t, 18);
		avgResults_MD_binary(outp_binary[16], DS->PrintVar_binary[16], DS, cD->surfDInt, DS->NumEle, t, 19);
	}	
	if (cD->prepDInt >0)
		{
			avgResults_MD_binary(outp_binary[10],DS->PrintVar_binary[10],DS,cD->prepDInt,DS->NumEle,t,12);
			avgResults_MD_binary(outp_binary[17],DS->PrintVar_binary[17],DS,cD->prepDInt,DS->NumEle,t,20);
		}
}


void
PrintData_LE(FILE ** outp_LE, Control_Data * cD, Model_Data DS, N_Vector CV_Y, realtype t)
{
	
	FILE *Oinit_file_LE;
	int i;
	if (DS->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		if (cD->groundelevDInt>0)
		{
			avgResults_NV(outp_LE[0], DS->PrintVar_LE[0], CV_Y, cD->groundelevDInt, DS->NumEle, t, 3 * DS->NumEle,DS);
		}
		if (cD->rockelevDInt>0)
		{
			avgResults_NV(outp_LE[1], DS->PrintVar_LE[1], CV_Y, cD->rockelevDInt, DS->NumEle, t, 4 * DS->NumEle,DS);
		}
	}
	else
	{
		if (cD->groundelevDInt>0)
		{
			avgResults_NV(outp_LE[0], DS->PrintVar_LE[0], CV_Y, cD->groundelevDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
		}
		if (cD->rockelevDInt>0)
		{
			avgResults_NV(outp_LE[1], DS->PrintVar_LE[1], CV_Y, cD->rockelevDInt, DS->NumEle, t, 2 * DS->NumEle,DS);
		}
	}
	if (cD->soilsubDInt>0)
		{
			avgResults_MD(outp_LE[2],DS->PrintVar_LE[2],DS,cD->soilsubDInt,DS->NumEle,t,23);
			avgResults_MD(outp_LE[3],DS->PrintVar_LE[3],DS,cD->soilsubDInt,DS->NumEle,t,24);
			avgResults_MD(outp_LE[4],DS->PrintVar_LE[4],DS,cD->soilsubDInt,DS->NumEle,t,25);
		}
	if (cD->weatheringDInt>0)
		{
			avgResults_MD(outp_LE[5],DS->PrintVar_LE[5],DS,cD->weatheringDInt,DS->NumEle,t,26);
		}
	if (cD->upliftDInt>0)
		{
			avgResults_MD(outp_LE[6],DS->PrintVar_LE[6],DS,cD->upliftDInt,DS->NumEle,t,27);
		}
	if (cD->bedloadDInt>0)
		{
			avgResults_MD(outp_LE[7],DS->PrintVar_LE[7],DS,cD->bedloadDInt,DS->NumEle,t,28);
			avgResults_MD(outp_LE[8],DS->PrintVar_LE[8],DS,cD->bedloadDInt,DS->NumEle,t,29);
			avgResults_MD(outp_LE[9],DS->PrintVar_LE[9],DS,cD->bedloadDInt,DS->NumEle,t,30);
		}
}


void
PrintData_LE_binary(FILE ** outp_LE, FILE ** outp_LE_binary, Control_Data * cD, Model_Data DS, N_Vector CV_Y, realtype t)
{
	
	FILE *Oinit_file_LE;
	int i;
	if (DS->Keyword_LE_PIHM_SUBHYDRO == 1)
	{
		if (cD->groundelevDInt>0)
		{
			avgResults_NV_binary(outp_LE_binary[0], DS->PrintVar_LE_binary[0], CV_Y, cD->groundelevDInt, DS->NumEle, t, 3 * DS->NumEle,DS);
			avgResults_NV(outp_LE[0], DS->PrintVar_LE[0], CV_Y, cD->groundelevDInt, DS->NumEle, t, 3 * DS->NumEle,DS);
		}
		if (cD->rockelevDInt>0)
		{
			avgResults_NV_binary(outp_LE_binary[1], DS->PrintVar_LE_binary[1], CV_Y, cD->rockelevDInt, DS->NumEle, t, 4 * DS->NumEle,DS);
			avgResults_NV(outp_LE[1], DS->PrintVar_LE[1], CV_Y, cD->rockelevDInt, DS->NumEle, t, 4 * DS->NumEle,DS);
		}
	}
	else
	{
		if (cD->groundelevDInt>0)
		{
			avgResults_NV_binary(outp_LE_binary[0], DS->PrintVar_LE_binary[0], CV_Y, cD->groundelevDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
		}
		if (cD->rockelevDInt>0)
		{
			avgResults_NV_binary(outp_LE_binary[1], DS->PrintVar_LE_binary[1], CV_Y, cD->rockelevDInt, DS->NumEle, t, 2 * DS->NumEle,DS);
		}
	}
	if (cD->soilsubDInt>0)
		{
			avgResults_MD_binary(outp_LE_binary[2],DS->PrintVar_LE_binary[2],DS,cD->soilsubDInt,DS->NumEle,t,23);
			avgResults_MD_binary(outp_LE_binary[3],DS->PrintVar_LE_binary[3],DS,cD->soilsubDInt,DS->NumEle,t,24);
			avgResults_MD_binary(outp_LE_binary[4],DS->PrintVar_LE_binary[4],DS,cD->soilsubDInt,DS->NumEle,t,25);
		}
	if (cD->weatheringDInt>0)
		{
			avgResults_MD_binary(outp_LE_binary[5],DS->PrintVar_LE_binary[5],DS,cD->weatheringDInt,DS->NumEle,t,26);
		}
	if (cD->upliftDInt>0)
		{
			avgResults_MD_binary(outp_LE_binary[6],DS->PrintVar_LE_binary[6],DS,cD->upliftDInt,DS->NumEle,t,27);
		}
	if (cD->bedloadDInt>0)
		{
			avgResults_MD_binary(outp_LE_binary[7],DS->PrintVar_LE_binary[7],DS,cD->bedloadDInt,DS->NumEle,t,28);
			avgResults_MD_binary(outp_LE_binary[8],DS->PrintVar_LE_binary[8],DS,cD->bedloadDInt,DS->NumEle,t,29);
			avgResults_MD_binary(outp_LE_binary[9],DS->PrintVar_LE_binary[9],DS,cD->bedloadDInt,DS->NumEle,t,30);
		}
}


void
PrintData_no_fully(FILE ** outp, Control_Data * cD, Model_Data DS, N_Vector CV_Y, realtype t)
{
	FILE *Oinit_file;
	int i,k;
	if (cD->gwDInt >0) {
		avgResults_NV(outp[0], DS->PrintVar[0], CV_Y, cD->gwDInt, DS->NumEle, t, 2 * DS->NumEle,DS);
	}
	if (cD->surfDInt >0) {
		avgResults_NV(outp[1], DS->PrintVar[1], CV_Y, cD->surfDInt, DS->NumEle, t, 0 * DS->NumEle,DS);
	}
	for (k = 0; k < 3; k++) {
		if (cD->etInt >0) {
			avgResults_MD(outp[2 + k], DS->PrintVar[2 + k], DS, cD->etInt, DS->NumEle, t, 3 + k);
		}
	}
	if (cD->IsDInt >0) {
		avgResults_MD(outp[5], DS->PrintVar[5], DS, cD->IsDInt, DS->NumEle, t, 6);
	}
	if (cD->snowDInt >0) {
		avgResults_MD(outp[6], DS->PrintVar[6], DS, cD->snowDInt, DS->NumEle, t, 7);
	}
	if (cD->usDInt >0) {
		avgResults_NV(outp[7], DS->PrintVar[7], CV_Y, cD->usDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
	}
	if (cD->RechInt >0) {
		avgResults_MD(outp[8], DS->PrintVar[8], DS, cD->RechInt, DS->NumEle, t, 8);
		avgResults_MD(outp[9], DS->PrintVar[9], DS, cD->RechInt, DS->NumEle, t, 10);
	}
	if (cD->surfDInt >0) {
		avgResults_MD(outp[11], DS->PrintVar[11], DS, cD->surfDInt, DS->NumEle, t, 14);
		avgResults_MD(outp[12], DS->PrintVar[12], DS, cD->surfDInt, DS->NumEle, t, 15);
		avgResults_MD(outp[13], DS->PrintVar[13], DS, cD->surfDInt, DS->NumEle, t, 16);
		avgResults_MD(outp[14], DS->PrintVar[14], DS, cD->surfDInt, DS->NumEle, t, 17);
		avgResults_MD(outp[15], DS->PrintVar[15], DS, cD->surfDInt, DS->NumEle, t, 18);
		avgResults_MD(outp[16], DS->PrintVar[16], DS, cD->surfDInt, DS->NumEle, t, 19);
	}	
	if (cD->prepDInt >0)
		{
			avgResults_MD(outp[10],DS->PrintVar[10],DS,cD->prepDInt,DS->NumEle,t,12);
			avgResults_MD(outp[17],DS->PrintVar[17],DS,cD->prepDInt,DS->NumEle,t,20);
		}
}

void
PrintData_no_fully_LE(FILE ** outp_LE, Control_Data * cD, Model_Data DS, N_Vector CV_Y, realtype t)
{
	
	FILE *Oinit_file_LE;
	int i;
	if (cD->groundelevDInt>0)
		{
			avgResults_NV(outp_LE[0], DS->PrintVar_LE[0], CV_Y, cD->groundelevDInt, DS->NumEle, t, 0 * DS->NumEle,DS);
		}
	if (cD->rockelevDInt>0)
		{
			avgResults_NV(outp_LE[1], DS->PrintVar_LE[1], CV_Y, cD->rockelevDInt, DS->NumEle, t, 1 * DS->NumEle,DS);
		}
		
	if (cD->soilsubDInt>0)
		{
			avgResults_MD(outp_LE[2],DS->PrintVar_LE[2],DS,cD->soilsubDInt,DS->NumEle,t,23);
			avgResults_MD(outp_LE[3],DS->PrintVar_LE[3],DS,cD->soilsubDInt,DS->NumEle,t,24);
			avgResults_MD(outp_LE[4],DS->PrintVar_LE[4],DS,cD->soilsubDInt,DS->NumEle,t,25);
		}
	if (cD->weatheringDInt>0)
		{
			avgResults_MD(outp_LE[5],DS->PrintVar_LE[5],DS,cD->weatheringDInt,DS->NumEle,t,26);
		}
	if (cD->upliftDInt>0)
		{
			avgResults_MD(outp_LE[6],DS->PrintVar_LE[6],DS,cD->upliftDInt,DS->NumEle,t,27);
		}
	if (cD->bedloadDInt>0)
		{
			avgResults_MD(outp_LE[7],DS->PrintVar_LE[7],DS,cD->bedloadDInt,DS->NumEle,t,28);
			avgResults_MD(outp_LE[8],DS->PrintVar_LE[8],DS,cD->bedloadDInt,DS->NumEle,t,29);
			avgResults_MD(outp_LE[9],DS->PrintVar_LE[9],DS,cD->bedloadDInt,DS->NumEle,t,30);
		}
}

void
copy_inputs(char *directory, char *filename)
{
	FILE *orig_file, *desti_file;
	char *temp_char,*temp_string;
	char ch;
    int flag, pos;
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	flag = mkdir(temp_string,0755);
/********* pihm.h**********/ 	
	orig_file = fopen("pihm.h","r");
	desti_file = fopen(strcat(temp_string,"pihm.h"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* pihm.c**********/ 
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("pihm.c","r");
	desti_file = fopen(strcat(temp_string,"pihm.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* f.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("f.c","r");
	desti_file = fopen(strcat(temp_string,"f.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
	/********* f_Hydro.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("f_Hydro.c","r");
	desti_file = fopen(strcat(temp_string,"f_Hydro.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
	/********* f_LE.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("f_LE.c","r");
	desti_file = fopen(strcat(temp_string,"f_LE.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* initialize.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char)); 
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("initialize.c","r");
	desti_file = fopen(strcat(temp_string,"initialize.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* is_sm_et.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("is_sm_et.c","r");
	desti_file = fopen(strcat(temp_string,"is_sm_et.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
  
/********* print.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("print.c","r");
	desti_file = fopen(strcat(temp_string,"print.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* read_alloc.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("read_alloc.c","r");
	desti_file = fopen(strcat(temp_string,"read_alloc.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* update.c**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("update.c","r");
	desti_file = fopen(strcat(temp_string,"update.c"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* Makefile**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("Makefile","r");
	desti_file = fopen(strcat(temp_string,"Makefile"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .calib**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.calib",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.calib",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .forc**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.forc",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.forc",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .geol**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.geol",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.geol",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .ibc**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.ibc",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.ibc",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .init**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.init",filename);
	orig_file = fopen(temp_string,"r");
	if (orig_file!=NULL)
	{
		free(temp_string);
		temp_string = (char *) malloc(1024 * sizeof(char));
		sprintf(temp_string,"%sinput/%s.init",directory,filename);
		desti_file = fopen(temp_string, "w");  
		fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
		pos = ftell(orig_file);
		fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
		while (pos--)
		{
			ch = fgetc(orig_file);  // copying file character by character
			fputc(ch, desti_file);
		}    
		fclose(orig_file);
	fclose(desti_file);
		free(temp_string);
	}
	
/********* .lc**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.lc",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.lc",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .mesh**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.mesh",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.mesh",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .att**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.att",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.att",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .initelev**********/
	
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.initelev",filename);
	orig_file = fopen(temp_string,"r");
	if (orig_file!=NULL)
	{
		free(temp_string);
		temp_string = (char *) malloc(1024 * sizeof(char));
		sprintf(temp_string,"%sinput/%s.initelev",directory,filename);
		desti_file = fopen(temp_string, "w");  
		fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
		pos = ftell(orig_file);
		fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
		while (pos--)
		{
			ch = fgetc(orig_file);  // copying file character by character
			fputc(ch, desti_file);
		}    
		fclose(orig_file);
		fclose(desti_file);
		free(temp_string);
	}
	
/********* .newelev**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.newelev",filename);
	orig_file = fopen(temp_string,"r");
	if (orig_file!=NULL)
	{
		free(temp_string);
		temp_string = (char *) malloc(1024 * sizeof(char));
		sprintf(temp_string,"%sinput/%s.newelev",directory,filename);
		desti_file = fopen(temp_string, "w");  
		fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
		pos = ftell(orig_file);
		fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
		while (pos--)
		{
			ch = fgetc(orig_file);  // copying file character by character
			fputc(ch, desti_file);
		}    
		fclose(orig_file);
		fclose(desti_file);
		free(temp_string);
	}
/********* .soil**********/	
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.soil",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.soil",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .LE_soil**********/	
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.LE_soil",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.LE_soil",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .LE_bedrock**********/	
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.LE_bedrock",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.LE_bedrock",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);

/********* projectName.txt**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("projectName.txt","r");
	desti_file = fopen(strcat(temp_string,"projectName.txt"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
  
/********* qsub.txt**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/",directory);
	orig_file = fopen("qsub.txt","r");
	desti_file = fopen(strcat(temp_string,"qsub.txt"), "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
  
/********* .para**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.para",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.para",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .LE_att**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.LE_att",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.LE_att",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
/********* .bid**********/
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%s.bid",filename);
	orig_file = fopen(temp_string,"r");
	free(temp_string);
	temp_string = (char *) malloc(1024 * sizeof(char));
	sprintf(temp_string,"%sinput/%s.bid",directory,filename);
	desti_file = fopen(temp_string, "w");  
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file);  // copying file character by character
        fputc(ch, desti_file);
    }    
    fclose(orig_file);
	fclose(desti_file);
	free(temp_string);
}
