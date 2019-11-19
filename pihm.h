/**********************************************************************************
 * File        : pihm.h                                                           *
 * Function    : Declaration and Definition of global variables and data structure*
 * Developer of PIHM 2.0: Mukesh Kumar (muk139@psu.edu)                           *
 * Developer of PIHM 1.0: Yizhong Qu   (quyizhong@gmail.com)                      *
 * Version     : Nov, 2007 (2.0)                                                  *
 *--------------------------------------------------------------------------------*
 *                                                                                *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0...............................*
 * a) Definition of new variables for ELEMENT data model related to 		  *
 *	i. Subsurface: KsatH, KsatV, infKsatV,infD, RzD, macD, macKsatH, macKsatV *
 *	vAreaF, hAreaF								  *
 *	ii. ET: LAImax, vegFrac, Albedo, Rs_ref, Rmin				  *
 *	iii. Snow: meltF							  *
 *	iv. Surface: dhBydx, dhBYdy, surfH					  *
 * b) Definition of new variables for RIVER data model				  *
 *	i. KsatH, KsatV, bedThick 						  *
 * c) Definition of new structures:						  *
 *	i. geology								  *
 *	ii. Land Cover								  *
 *	iii. Calibration							  *
 * d) Definition of New Control Parameters for each file			  *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                	  *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 **********************************************************************************/

#include "sundials_types.h"
#include "nvector_serial.h"

/* Definition of Global Variable Types */

// #define FULLY_COUPLE
// //#define LE_PIHM
// //#define BEDROCK
// //#define LE_PIHM_SED
// #define LE_PIHM_HYDRO
// #define LE_PIHM_SUBHYDRO
// //#define REALTIME
// #define MEAN
		


typedef struct element_type {	/* Data model for a triangular element */
	int             index;	/* Element No. */
	int             node[3];/* anti-clock-wise */
	int             nabr[3];/* neighbor i shares edge i (0: on boundary) */

	realtype        edge[3];/* edge i is from node i to node i+1 */
	realtype        area;	/* area of element */

	realtype        x;	/* x of centroid */
	realtype        y;	/* y of centroid */
	realtype        zmin;	/* z_min of centroid */
	realtype        zmax;	/* z_max of centroid */
	

	realtype        KsatH;	/* horizontal geologic saturated hydraulic
				 * conductivity */
	realtype        KsatV;	/* vertical geologic saturated hydraulic
				 * conductivity */
	realtype        infKsatV;	/* vertical surface saturated
					 * hydraulic conductivity */
	realtype        Porosity;
	realtype        infD;	/* depth from ground surface accross which
				 * head is calculated during infiltration */
	realtype        Alpha;	/* Alpha from van-genuchten eqn which is
				 * given by satn =
				 * 1/pow(1+pow(abs(Alpha*psi),Beta),1-1/Beta) */
	realtype        Beta;
	realtype        RzD;	/* Root zone depth */
	realtype        macD;	/* macropore Depth */
	realtype        macKsatH;	/* macropore horizontal saturated
					 * hydraulic conductivity */
	realtype        macKsatV;	/* macropore vertical saturated
					 * hydraulic conductivity */
	realtype        vAreaF;	/* macropore area fraction on a vertical
				 * cross-section */
	realtype        hAreaF;	/* macropore area fraction on a horizontal
				 * cross-section */
	int             Macropore;	/* 1: macropore; 0: regular soil */

	realtype        LAImax;	/* maxm. LAI accross all seasons for a
				 * vegetation type */
	realtype        VegFrac;/* areal vegetation fraction in a triangular
				 * element */
	realtype        Albedo;	/* albedo of a triangular element */
	realtype        Rs_ref;	/* reference incoming solar flux for
				 * photosynthetically active canopy */
	realtype        Rmin;	/* minimum canopy resistance */
	realtype        Rough;	/* surface roughness of an element */

	realtype        windH;	/* wind measurement height */


	int             soil;	/* soil type */
	int             geol;	/* geology type */
	int             LC;	/* Land Cover type  */
	int             IC;	/* initial condition type */
	int             BC[3];	/* boundary type. 0:natural bc (no flow);
				 * 1:Dirichlet BC; 2:Neumann BC */
	int             prep;	/* precipitation (forcing) type */
	int             temp;	/* temperature (forcing) type   */
	int             humidity;	/* humidity type */
	int             WindVel;/* wind velocity type  */
	int             Rn;	/* net radiation input */
	int             G;	/* radiation into ground */
	int             pressure;	/* pressure type */
	int             source;	/* source (well) type */
	int             meltF;	/* meltFactor */
	/* for calculation of dh/ds */
	realtype        surfH[3];	/* Total head in neighboring cells */
	realtype        surfX[3];	/* Center X location of neighboring
					 * cells */
	realtype        surfY[3];	/* Center Y location of neighboring
					 * cells */
	realtype        dhBYdx;	/* Head gradient in x dirn. */
	realtype        dhBYdy;	/* Head gradient in y dirn. */
	int	 			*Nabr_index; /* the index of the Nabour element, 0,1 or 2*/
	
	//#ifdef LE_PIHM
		realtype        DReg;   /*diffusivity of regolith*/
		realtype        P0; /*Regolith production rate in the absence of soil above bedrock*/
		realtype        Uplift; /*Bedrock uplift rate*/
		realtype        CoefP0; /*Fitting constant in weather equation*/
		realtype		RhoReg;
		realtype		RhoBed;
		realtype		Diameter; //The diameter of sediment particle
		realtype		Tan_slope;
		int 			LE_soil;
		int				LE_geol;
	//#endif
	
	// #ifdef LE_PIHM_SED
		// realtype		RhoReg;
		// realtype		Diameter; //The diameter of sediment particle
	// #endif
	
}               element;

typedef struct nodes_type {	/* Data model for a node */
	int             index;	/* Node no. */

	realtype        x;	/* x coordinate */
	realtype        y;	/* y coordinate */
	realtype        zmin;	/* z bed rock elevation */
	realtype        zmax;	/* z surface elevation  */

}               nodes;

typedef struct element_IC_type {/* Initial state variable conditions on each
				 * element */
	int             index;

	realtype        interception;	/* Interception storage (Note all
					 * these variables have dimension of
					 * L */
	realtype        snow;	/* Snow depth */
	realtype        surf;	/* Overland flow depth */
	realtype        unsat;	/* unsaturated zone depth */
	realtype        sat;	/* saturated zone depth */
	
	//#ifdef LE_PIHM
		realtype        groundelev;	/* elevation of ground surface*/
		realtype        bedrockelev;	/* elevation of bedrock */
	//#endif
	
	// #ifdef LE_PIHM_SED
		// realtype        groundelev;	/* elevation of ground surface*/
	// #endif

}               element_IC;

typedef struct soils_type {
	int             index;	/* index */

	realtype        KsatV;	/* vertical saturated soil conductivity */
	realtype        ThetaS;	/* soil porosity */
	realtype        ThetaR;	/* soil moisture residual */
	realtype        Alpha;	/* soil curve parameter 1 */
	realtype        Beta;	/* soil curve parameter 2 */

	realtype        hAreaF;	/* macroporous area fraction on horizontal
				 * section */
	realtype        macKsatV;	/* macroporous saturated vertical
					 * conductivity */

	realtype        infD;	/* depth from ground surface accross which
				 * head is calculated during infiltration */
}               soils;

typedef struct LE_soils {
	int             index;	/* index */

	//#ifdef LE_PIHM
		realtype        RhoReg; /*Density of regolith*/
		realtype        DReg; /*Diffusivity*/
		realtype		Diameter; //The diameter of sediment particle
	//#endif
	
	// #ifdef LE_PIHM_SED
		// realtype        RhoReg; /*Density of regolith*/
		// realtype		Diameter; //The diameter of sediment particle
	// #endif
}               LE_soils;

typedef struct geol_type {
	int             index;	/* index */

	realtype        KsatH;	/* horizontal saturated geology conductivity */
	realtype        KsatV;	/* vertical saturated geology conductivity */
	realtype        ThetaS;	/* geology porosity */
	realtype        ThetaR;	/* residual porosity */
	realtype        Alpha;	/* van genuchten parameter */
	realtype        Beta;	/* van genuchten parameter */

	realtype        vAreaF;	/* macroporous area fraction on vertical
				 * section */
	realtype        macKsatH;	/* macroporous saturated horizontal
					 * conductivity */
	realtype        macD;
	
	// #ifdef LE_PIHM_SED
		
	// #endif
	
	

}               geol;

typedef struct LE_bedrock

{
	int             index;	/* index */
	//#ifdef LE_PIHM
		realtype        RhoBed; /*Density of bedrock*/
		realtype        P0; /*Regolith production rate in the absence of soil above bedrock*/
		realtype        Uplift; /*Bedrock uplift rate*/
		realtype        CoefP0; /*Fitting constant in weather equation*/
	//#endif	
	// #ifdef LE_PIHM_SED		
	// #endif
}               LE_bedrock;


typedef struct lc_type {
	int             index;	/* index */

	realtype        LAImax;	/* max LAI */
	realtype        VegFrac;/* Canopy Fracn */
	realtype        Albedo;	/* Albedo */
	realtype        Rs_ref;
	realtype        Rmin;	/* Minimum stomatal resistance */
	realtype        Rough;	/* surface roughness factor  */
	realtype        RzD;	/* rootZone Depth */
}               LC;

typedef struct TSD_type {
	char            name[5];
	int             index;
	int             length;	/* length of time series */
	int             iCounter;	/* interpolation counter */
	realtype      **TS;	/* 2D time series data */

}               TSD;

typedef struct global_calib {
	realtype        KsatH;	/* For explanation of each calibration
				 * variable, look for corresponding variables
				 * above */
	realtype        KsatV;
	realtype        infKsatV;
	realtype        macKsatH;
	realtype        macKsatV;
	realtype        infD;
	realtype        RzD;
	realtype        macD;
	realtype        Porosity;
	realtype        Alpha;
	realtype        Beta;
	realtype        vAreaF;
	realtype        hAreaF;
	realtype        Temp;
	realtype        Prep;
	realtype        VegFrac;
	realtype        Albedo;
	realtype        Rough;
	
	//#ifdef LE_PIHM
		realtype        DReg;   /*diffusivity of regolith*/
		realtype        P0; /*Regolith production rate in the absence of soil above bedrock*/
		realtype        CoefP0; /*Fitting constant in weather equation*/
	//#endif
	
	// #ifdef LE_PIHM_SED
		
	// #endif
	
}               globalCal;

typedef struct process_control {
	realtype        Et0;
	realtype        Et1;
	realtype        Et2;
}               processCal;

typedef struct model_data_structure {	/* Model_data definition */
	int             UnsatMode;	/* Unsat Mode */
	int             SurfMode;	/* Surface Overland Flow Mode */

	int             NumEle;	/* Number of Elements */
	int             NumNode;/* Number of Nodes    */

	int             NumPrep;/* Number of Precipatation time series types  */
	int             NumTemp;/* Number of Temperature time series types      */
	int             NumHumidity;	/* Number of Humidity time series
					 * types         */
	int             NumWindVel;	/* Number of Wind Velocity time
					 * series types    */
	int             NumRn;	/* Number of Net Radiation time series types    */
	int             NumG;	/* Number of Ground Heat time series types      */
	int             NumP;	/* Number of Pressure time series types         */
	int             NumSource;	/* Number of Source time series types           */

	int             NumSoil;/* Number of Soils           */
	int             NumLESoil;/* Number of Soils           */
	int             NumLEBedrock;/* Number of Soils           */
	int             NumGeol;/* Number of Geologies           */
	int             NumRes;	/* Number of Reservoir       */
	int             NumLC;	/* Number of Land Cover Index Data */

	int             NumMeltF;	/* Number of Melt Factor Time series */

	int             Num1BC;	/* Number of Dirichlet BC    */
	int             Num2BC;	/* Number of Numann BC       */
	int             NumEleIC;	/* Number of Element Initial Condtion */

	int	   			NumOutlet; /*Number of element at outlet */
	double          MSF;	/* Morphological Scaling Factor */
	int 			MSF_multiplyer; /*multiplyer for MSF */
	int 			Dynamic_interval; /*The interval of using  MSF_multiplyer*/
	int				computtime; 

	element        *Ele;	/* Store Element Information  */
	nodes          *Node;	/* Store Node Information     */
	element_IC     *Ele_IC;	/* Store Element Initial Condtion */
	soils          *Soil;	/* Store Soil Information     */
	LE_soils	   *LE_soil;
	geol           *Geol;	/* Store Geology Information     */
	LE_bedrock	   *LE_bedrock;
	LC             *LandC;	/* Store Land Cover Information */

	TSD            *TSD_Inc;/* Infiltration Capacity Time Series Data */
	TSD            *TSD_LAI;/* Leaves Area Index Time Series Data     */
	              //TSD * TSD_DH;	/* Zero plane Displacement Height */
	TSD            *TSD_RL;	/* Roughness Length */
	realtype       *ISFactor;	/* ISFactor is used to calculate
					 * ISMax from LAI */
	realtype       *windH;	/* Height at which wind velocity is measured */
	TSD            *TSD_MeltF;	/* Monthly Varying Melt Factor for
					 * Temperature Index model */

	TSD            *TSD_EleBC;	/* Element Boundary Condition Time
					 * Series Data  */
	TSD            *TSD_Prep;	/* RainFall Time Series Data       */
	TSD            *TSD_Temp;	/* Temperature Time Series Data    */
	TSD            *TSD_Humidity;	/* Humidity Time Series Data       */
	TSD            *TSD_WindVel;	/* Wind Velocity Time Series Data  */
	TSD            *TSD_Rn;	/* Net Radiation Time Series Data  */
	TSD            *TSD_G;	/* Radiation into Ground Time Series Data */
	TSD            *TSD_Pressure;	/* Vapor Pressure Time Series data       */
	TSD            *TSD_Source;	/* Source (well) Time Series data  */

	realtype      **FluxSurf;	/* Overland Flux   */
	realtype      **FluxSub;	/* GW Flux   */
	realtype      **Seepage;	/* GW Flux   */
	
	realtype       *ElePrep;/* Precep. on each element */
	realtype       *EleETloss;
	realtype       *EleNetPrep;	/* Net precep. on each elment */
	realtype       *EleViR;	/* Variable infiltration rate */
	realtype       *Recharge;	/* Recharge rate to GW */
	realtype       *EleSnow;/* Snow depth on each element */
	realtype       *EleSnowGrnd;	/* Snow depth on ground element */
	realtype       *EleSnowCanopy;	/* Snow depth on canopy element */
	realtype       *EleIS;	/* Interception storage */
	realtype       *EleISmax;	/* Maximum interception storage
					 * (liquid precep) */
	realtype       *EleISsnowmax;	/* Maximum interception storage
					 * (snow) */
	realtype       *EleTF;	/* Through Fall */
	realtype      **EleET;	/* Evapo-transpiration (from canopy, ground,
				 * subsurface, transpiration) */
    realtype		**DY; /*store the DY for surface water, surface elevation and bed elevation*/
	realtype       *DummyY;
	realtype       *DummyY_Hydro;
	realtype       *DummyY_LE;
	realtype       *PrintVar[18];
	realtype       *PrintVar_binary[18];
	realtype	   *Total_Surfwater_out;
	realtype	   *Total_Unsatwater_out;
	realtype	   *Total_Satwater_out;
	realtype	   *infil_mode;
	int	   **Outlet_location;
	processCal      pcCal;
	realtype *SurfDepth;     /*Surface water*/
	realtype *GW;      /*Ground water*/
	realtype *UGW;     /*Unsaturated water*/
	realtype *grdelev;     /*Unsaturated water*/
	realtype *bedelev;     /*Unsaturated water*/
	realtype		t; /*current simulation time*/
	realtype		Constant_precip;
	int 			SedMode;
	//#ifdef LE_PIHM
		realtype	  *Total_soil_out;
		realtype	  *Total_soil_in;
		realtype      *PrintVar_LE[10];
		realtype      *PrintVar_LE_binary[10];
		realtype      **SoilFluxSurf;	/* Overland Soil Flux   */
		realtype      **SoilFluxSub;/* Subsurface Soil Flux */
		realtype	  **BedloadFlux; /* Bedload sediment transport flux */
		realtype	  *Soilproduction; /*Soil production*/
		realtype	  *Uplift;//Bedrock uplift
		realtype      *zmax_init;	/* init elevation of z_max of centroid, used to calculate boundary problem*/
		realtype      DReg;   /*diffusivity of regolith*/
		realtype      P0; /*Regolith production rate in the absence of soil above bedrock*/
		realtype      CoefP0; /*Fitting constant in weather equation*/
	//#endif
	int				Keyword_FULLY_COUPLE;
	int				Keyword_LE_PIHM;
	int				Keyword_BEDROCK;
	int				Keyword_LE_PIHM_SED;
	int				Keyword_LE_PIHM_HYDRO;
	int				Keyword_LE_PIHM_SUBHYDRO;
	int				Keyword_REALTIME;
	int				Keyword_MEAN;
	int				Keyword_BINARY;
}              *Model_Data;

typedef struct control_data_structure {
	int             Verbose;
	int             Debug;

	int             Solver;	/* Solver type */
	int             NumSteps;	/* Number of external time steps
					 * (when results can be printed) for
					 * the whole simulation */

	int             gwD;	/* File boolean, Choose 1 if your want to
				 * print ground water */
	int prepD;  /* File boolean, Choose 1 if your want to print precipitation */
	int             surfD;	/* File boolean, Choose 1 if your want to
				 * print overland flow */
	int             snowD;	/* File boolean, Choose 1 if your want to
				 * print snow Depth */
	int             Rech;	/* File boolean, Choose 1 if your want to
				 * print recharge to ground water */
	int             IsD;	/* File boolean, Choose 1 if your want to
				 * print interception depth */
	int             usD;	/* File boolean, Choose 1 if your want to
				 * print unsaturated depth */
	int             et[3];	/* File boolean, Choose 1 if your want to
				 * print individual evapo-transpiration
				 * components */
	
	int             rivFlx;	/* File boolean, Choose 1 if your
					 * want to print river/river bed
					 * fluxes */

	int             gwDInt;	/* Time interval to output average val of
				 * variables */
	int 			prepDInt;
	int             surfDInt;
	int             snowDInt;
	int             RechInt;
	int             IsDInt;
	int             usDInt;
	int             etInt;
	int				computtimeDInt;
	
	int             init_type;	/* initialization mode */
	
	

	realtype        abstol;	/* absolute tolerance */
	realtype        reltol;	/* relative tolerance */
	realtype        InitStep;	/* initial step size */
	realtype        MaxStep;/* Maximum step size */
	realtype        ETStep;	/* Step for et from interception */

	int             GSType, MaxK;	/* Maximum Krylov order */
	realtype        delt;

	realtype        StartTime;	/* Start time of simulation */
	realtype        EndTime;/* End time of simulation */
	realtype		FinishTime; /*Time of the total run, if FinishTim is greater than EndTime, time t will start from StartTime again. */


	int             outtype;
	realtype        a;	/* External time stepping controls */
	realtype        b;

	realtype       *Tout;

	globalCal       Cal;	/* Convert this to pointer for localized
				 * calibration */
	//#ifdef LE_PIHM

		int				groundelevDInt;
		int				rockelevDInt;
		int				bedloadDInt;
		int				soilsubDInt;
		int				upliftDInt;
		int				weatheringDInt;
		int				continue_LEM; /* Use previous results as initial condition running Landscape evolution model again. 0 don't use, 1 use.*/
	//#endif

	
	// #ifdef LE_PIHM_SED
		
	// #endif
}               Control_Data;

/* Function Declarations */
void            initialize(char *, Model_Data, Control_Data *, N_Vector);
void			initialize_no_fully(char *, Model_Data, Control_Data *, N_Vector, N_Vector);
void            is_sm_et(realtype, realtype, Model_Data, N_Vector);
/* Function to calculate right hand side of ODE systems */
int             f(realtype, N_Vector, N_Vector, void *);
realtype* initiation(Model_Data,realtype*,realtype*, realtype);
void 			HydroSubLateral(Model_Data, int,int,int);
void 			HydroSurfLateral(Model_Data, int,int,int);
void 			SedDiffu(Model_Data, int,int,int);
void 			ET(Model_Data, int,realtype);
void 			Infiltration(Model_Data, int);
void 			Adjust_water(Model_Data);
void 			SedAdvec(Model_Data, int,int,int);
void 			Adjust_sed(Model_Data);
realtype* 		Calculate_DY(Model_Data,realtype*);
void 			SoilProduction(Model_Data,  int);
void 			BedUplift(Model_Data, int);

int             f_LE(realtype, N_Vector, N_Vector, void *);
int             f_Hydro(realtype, N_Vector, N_Vector, void *);
void            read_alloc(char *, Model_Data, Control_Data *);	/* Variable definition */
void            update(realtype, void *);
void            PrintData(FILE **, Control_Data *, Model_Data, N_Vector, realtype);
void            PrintData_LE(FILE **, Control_Data *, Model_Data, N_Vector, realtype);
void            PrintData_binary(FILE **,FILE **, Control_Data *, Model_Data, N_Vector, realtype);
void            PrintData_LE_binary(FILE **,FILE **, Control_Data *, Model_Data, N_Vector, realtype);
void			PrintData_no_fully(FILE **, Control_Data *, Model_Data, N_Vector, realtype);
void			PrintData_no_fully_LE(FILE **, Control_Data *, Model_Data, N_Vector, realtype);
void            FreeData(Model_Data, Control_Data *);
void 			copy_inputs(char *, char *);
realtype        Interpolation(TSD * Data, realtype t);
realtype 		avgY(realtype diff, realtype yi, realtype yinabr);
realtype 		effKV(realtype ksatFunc, realtype gradY, realtype macKV, realtype KV, realtype areaF);
realtype 		effKH(int mp, realtype tmpY, realtype aqDepth, realtype MacD, realtype MacKsatH, realtype areaF, realtype ksatH);