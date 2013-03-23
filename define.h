#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include<time.h>

//physical constants:
#define Pi             (3.1415926)
#define yr             (3.1536e7)        /****s***/
#define myr            ((1e6)*yr)
#define kyr            (1e3*yr)
#define Msolar         (1.989e30)        /****kg***/
#define pc             (3.1e16)          /****m****/
#define Grav           (6.67e-11)        /****the constant of gravity****/

//observables:
#define mBE            (2*Msolar)  //remember to change function frag_rate_gen and sf_rate_gen when change mBE!!!!
#define mBEnew         (0.84*Msolar)
#define v_m			   (-0.2)        //initial core v~m^0.1
//#define Rb             (5*pc)   //for concentrated, 2.75pc; for realistic, 5pc; now use total mass instead.
#define Backmassinsolar  (1e3)
#define Tmass_core0solar (350.0)
#define Totalmass_core0 (Tmass_core0solar*Msolar)
#define vc             (1e3)  //(2.0e3)   //before 0209-05(3e3)            /****velocity dispersion of cores.m/s,¾ßÌåÊý×Ö´ý¶¨****/
#define v0             (0.8e3)      //before 0209-05(4e3)         /****thermal velocity of molecules in background******/
#define rhoc           (3e-17)  //(3.187e-17)      /****kg/m^3****/
#define Drho            15 //7              /******density contrast, 100 or 1000****/
//#define fifa            0.0006             /***initial filling factor***/
/**for realistic, 0.01/Drho=0.0006; for concentrated, just keep the total core mass~150Msolar;**/
#define Rhob           (rhoc/Drho)

#define alpha_init      (0.7)       //in initialize(); dn/dlogm=A*m^(-alpha_init);
#define alpha           0.3//0.5        //in star_formation(); the efficiency of sf;
#define beta            (-2.5)		   //in coagulation(); the modified coefficient of cross section;
#define beta_frag       (5.0)      //in fragmentation(); fragment spectra powerlaw index;
#define kappa           0.8 //0.5        //in conden_nucle(); modified coefficient;         
#define gamma           0.5        //in sf_rate_gen(); modified coefficient;
#define co_frag         4.0        //in frag_rate_gen; modified coefficient;
#define co_evap			2.0  //2.0  //8.0		   //in evapration; powerlaw index;
#define evap_mod        3 //0.5       //in evapration; coefficient;
#define q_beta		    (-1.2) 
#define lambda_comp     (-1.5)      // in companion IMF, the IMF slope for the secondary stars;
#define idens 40				   // mean rhoc at m = 1.05Msun
#define X_fac				10			//coefficient for cooling function;
#define NUCLEAT		   (3.36e-2/kyr)   //coefficient for nucleation;

//#define R_H				(1e-10)	   //hydrogen radius
//#define mu_H			(1.6736e-27)		   //hydrogen molecular mass

//initial values:
//#define Backmassinsolar (4*Pi*rhoc*pow(Rb,3)/(3*Msolar*Drho))
#define Back_mass       (Backmassinsolar*Msolar)
#define Rb              (pow(Back_mass*3/(4*Pi*Rhob),1.0/3.0))
#define fifa            (Totalmass_core0/Backmassinsolar)
//#define Totalmass_core0 (Back_mass*fifa*Drho)

//derived constants:
#define Vcircular0      (pow(Grav*Back_mass/Rb, 0.5)) /*****circular speed******/

//runtime parameters:

#define unit           (20000)	//timestep interval of output record; 
#define Timestep       (10*yr)
#define totaltime      (2.0*myr)
#define bin             40          //parameter for initial mass distribution: truncation bin number;
#define SF_rate         (0.03/myr)           //in sf_rate_gen();
#define Frag_rate       1.0			//in frag_rate_gen();
#define Evap_rate       ((7e-3)/kyr)     //calculated according to energy conservation; try /kyr and /yr both;
//outdated !! #define cn_ratio        0.8              /****condensation to nucleation ratio***/
#define Gauss1          1                /****sigma for sf_rate_gen()*****/
#define Gauss2          20                /****sigma for frag_rate_gen()****/




/************function declaration***************/

	int coagulation(int i, int j, double ni, double nj, double ri, double rj);
	int fragmentation(int i);
	int star_formation_binary(int i);
	int star_formation(int i);
	int star_formation_companion(int i);
	int star_formation_hierarchy(int i);
	int evaporation();
	int conden_nucle();
	int update();    
	int initialize(double *time, double *TimeContinue, int *ff, int *ff2);
	int filename_gen(char f_txt[25], char *f_binary, char*f_distr, char*f_star, char* f_companion, char*cp);
	int examine(double *n, double *star_mass, double *back_mass, double *total_mass, double *total_mass_core);
    int finalize(double time, FILE *fp, char *cp, int ff, int ff2);
	int dataoutput2(double result[15][100], double result_star[15][100], double Ek_evo[15][100], int countcolumn, char *f_distr, char *f_star, char *f_Ek, int countcolumncf, char *filefrag, char*filecoag, double frag[15][100], double coag[15][100], double TimeContinue);
	int readdata(double *time, double *TimeContinue);
	int print_paras(FILE *fp, double back_mass, double TimeContinue, double total_mass_core);
	void sf_rate_gen(double *sf_rate, double *m);     /*to generate the star formation rate for every mass bin;*/
    double *r_gen(double *m, double *dens);                 /*to generate the radius of every mass bin;*/
	double *dens_gen(double *m);						 /*to generate the density of every mass bin*/
	void taoKH_gen(double *taoKH, double *r);  
	double *frag_rate_gen(double *m);
	double *evap_rate_gen();
	int output_star(FILE* fstar);
	int output_binary(FILE* fbinary);
	int output_companion(FILE* fcompanion);
	int output_distr(FILE* fdistr);

/**********************************/

/***************external variables*******************/

//main solution array:
double n[100], delta_n[100], rk[4][100], Ek[100], delta_Ek[100]; //rk[4][100] is for RK4 algorithm;
double n_star[100], n_binary[100], n_star_sec[100], n_star_prim[100];

//constant array:
double m[101], m_insolar[100], sum[100]/*** auxiliary array for m[]***/, v_d[100]/*** core velocity ***/, sf_rate[100];
double *evap_rate, *dens, *frag_rate/**smoothing function for frag*/; /***initially used, then don't need to update***/
int Evap_lim;


//need to update every timestep!!!!!!!!:
double *r, taoKH[100];  /**** auxiliary for calculation ****/
double delta_n_cond[100], delta_n_coag[100], delta_n_frag[100], delta_n_evap[100];
double rhob;   /*** rhob is variable with evolution ***/
double star_mass, back_mass, total_mass, total_mass_core, ghost_mass, totEk;

//runtime:
clock_t start, finish, start2, end2;
char *name;







