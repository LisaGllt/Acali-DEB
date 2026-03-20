/* /Users/lisagollot/Library/CloudStorage/OneDrive-Personnel/Documents/0_These/0_RepoGit/ew-deb-tktd/mod/DEB_deSolve/DEBcali.c for R deSolve package
   ___________________________________________________

   Model File:  /Users/lisagollot/Library/CloudStorage/OneDrive-Personnel/Documents/0_These/0_RepoGit/ew-deb-tktd/mod/DEB_deSolve/DEBcali.model

   Date:  Mon Feb  2 14:47:24 2026

   Created by:  "/Users/lisagollot/mcsim-6.2.0/mod/mod v6.2.0"
    -- a model preprocessor by Don Maszle
   ___________________________________________________

   Copyright (c) 1993-2020 Free Software Foundation, Inc.

   Model calculations for compartmental model:

   8 States:
     E = 0.0,
     L = 0.0,
     Eh = 0.0,
     R = 0.0,
     OM_soil = 0.0,
     OM_horse = 0.0,
     W = 0.0,
     OM = 0.0,

   5 Outputs:
    "Energy",
    "Maturity",
    "Reproduction",
    "Weight",
    "Organic_matter",

   1 Input:
     Dens (forcing function)

   25 Parameters:
     muOM = 0,
     rOM_ClxHorse = 0,
     Fm = 0,
     kapX = 0,
     pAm = 0,
     v = 0,
     kap = 0,
     pM = 0,
     pT = 0,
     Eg = 0,
     Shape = 0,
     w = 0,
     kJ = 0,
     Ehb = 0,
     Ehp = 0,
     L_coc = 0,
     E_coc = 0,
     kapR = 0,
     TAH = 0,
     TH = 0,
     TA = 0,
     Tref = 0,
     Winit = 0,
     Texp = 0,
     WeightSoilCosm = 0,
*/
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Model variables: States */
#define ID_E 0x00000
#define ID_L 0x00001
#define ID_Eh 0x00002
#define ID_R 0x00003
#define ID_OM_soil 0x00004
#define ID_OM_horse 0x00005
#define ID_W 0x00006
#define ID_OM 0x00007

/* Model variables: Outputs */
#define ID_Energy 0x00000
#define ID_Maturity 0x00001
#define ID_Reproduction 0x00002
#define ID_Weight 0x00003
#define ID_Organic_matter 0x00004

/* Parameters */
static double parms[25];

#define muOM parms[0]
#define rOM_ClxHorse parms[1]
#define Fm parms[2]
#define kapX parms[3]
#define pAm parms[4]
#define v parms[5]
#define kap parms[6]
#define pM parms[7]
#define pT parms[8]
#define Eg parms[9]
#define Shape parms[10]
#define w parms[11]
#define kJ parms[12]
#define Ehb parms[13]
#define Ehp parms[14]
#define L_coc parms[15]
#define E_coc parms[16]
#define kapR parms[17]
#define TAH parms[18]
#define TH parms[19]
#define TA parms[20]
#define Tref parms[21]
#define Winit parms[22]
#define Texp parms[23]
#define WeightSoilCosm parms[24]

/* Forcing (Input) functions */
static double forc[1];

#define Dens forc[0]

/* Function definitions for delay differential equations */

int Nout=1;
int nr[1]={0};
double ytau[1] = {0.0};

static double yini[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /*Array of initial state variables*/

void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

double CalcDelay(int hvar, double dTime, double delay) {
  double T = dTime-delay;
  if (dTime > delay){
    nr[0] = hvar;
    lagvalue( T, nr, Nout, ytau );
}
  else{
    ytau[0] = yini[hvar];
}
  return(ytau[0]);
}

/*----- Initializers */
void initmod (void (* odeparms)(int *, double *))
{
  int N=25;
  odeparms(&N, parms);
}

void initforc (void (* odeforcs)(int *, double *))
{
  int N=1;
  odeforcs(&N, forc);
}


/* Calling R code will ensure that input y has same
   dimension as yini */
void initState (double *y)
{
  int i;

  for (i = 0; i < sizeof(yini) / sizeof(yini[0]); i++)
  {
    yini[i] = y[i];
  }
}

void getParms (double *inParms, double *out, int *nout) {
/*----- Model scaling */

  /* local */ double kapreal;
  /* local */ double sA;
  /* local */ double srH;
  /* local */ double Tc;
  /* local */ double pAm_t;
  /* local */ double v_t;
  /* local */ double kJ_t;
  /* local */ double Km;
  /* local */ double Km_t;
  /* local */ double Em;
  /* local */ double g;
  /* local */ double Lm;
  /* local */ double freal;
  int i;

  for (i = 0; i < *nout; i++) {
    parms[i] = inParms[i];
  }


  for (i = 0; i < *nout; i++) {
    out[i] = parms[i];
  }
  }
/*----- Dynamics section */

void derivs (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{
  /* local */ double sA;
  /* local */ double srH;
  /* local */ double Tc;
  /* local */ double pAm_t;
  /* local */ double v_t;
  /* local */ double kJ_t;
  /* local */ double kapreal;
  /* local */ double Km;
  /* local */ double Km_t;
  /* local */ double Em;
  /* local */ double g;
  /* local */ double Lm;
  /* local */ double Qfmax;
  /* local */ double Qf_soil;
  /* local */ double Qf_horse;
  /* local */ double Qf;
  /* local */ double OM_soil_out;
  /* local */ double OM_horse_out;
  /* local */ double X;
  /* local */ double freal;
  /* local */ double tmp_E;
  /* local */ double pC;
  /* local */ double tmp_Eh;
  /* local */ double tmp_R;

  sA = exp ( TA / Tref - TA / ( Texp +273.15 ) ) ;

  srH = ( 1 + exp ( TAH / TH - TAH / Tref ) ) / ( 1 + exp ( TAH / TH - TAH / ( Texp +273.15 ) ) ) ;

  Tc = sA * ( ( Texp +273.15 >= Tref ) * srH + ( Texp +273.15 < Tref ) ) ;

  pAm_t = pAm * Tc ;

  v_t = v * Tc ;

  kJ_t = kJ * Tc ;

  kapreal = ( y[ID_Eh] < Ehp ? kap : ( Dens > 1 ? ( 1 - kap ) : kap ) ) ;

  Km = pM / Eg ;

  Km_t = Km * Tc ;

  Em = pAm / v ;

  g = Eg / ( kapreal * Em ) ;

  Lm = ( v / ( Km * g ) ) ;

  Qfmax = pAm_t * pow ( y[ID_L] , 2 ) / kapX ;

  Qf_soil = Fm * yout[ID_Weight] * ( muOM * rOM_ClxHorse * y[ID_OM_soil] / Dens ) ;

  Qf_horse = Fm * yout[ID_Weight] * ( muOM * y[ID_OM_horse] / Dens ) ;

  Qf = Qf_soil + Qf_horse ;

  OM_soil_out = ( Qf_soil / ( muOM * rOM_ClxHorse ) ) * Dens ;

  OM_horse_out = ( Qf_horse / muOM ) * Dens ;

  ydot[ID_OM_horse] = ( ( y[ID_OM_horse] - OM_horse_out ) < 1e-12 ? ( - y[ID_OM_horse] ) : ( - OM_horse_out ) ) ;

  ydot[ID_OM_soil] = ( ( y[ID_OM_soil] - OM_soil_out ) < 1e-12 ? ( - y[ID_OM_soil] ) : ( - OM_soil_out ) ) ;

  y[ID_OM] = y[ID_OM_soil] + y[ID_OM_horse] ;

  X = Qf / ( 0.5 * Qfmax ) ;

  freal = X / ( 1 + X ) ;

  tmp_E = ( pAm_t / y[ID_L] ) * ( freal - y[ID_E] / Em ) ;

  ydot[ID_E] = ( y[ID_E] < 1e-12 ? ( - y[ID_E] ) : tmp_E ) ;

  ydot[ID_L] = ( v_t / ( 3 * ( ( ( y[ID_E] / Em ) + g ) ) ) ) * ( ( y[ID_E] / Em ) - ( y[ID_L] / Lm ) ) ;

  y[ID_W] = pow ( y[ID_L] , 3 ) * ( 1 + y[ID_E] / Em * w ) ;

  pC = ( ( g * y[ID_E] ) / ( g + ( y[ID_E] / Em ) ) ) * ( ( v_t * pow ( y[ID_L] , 2 ) ) + ( Km_t * pow ( y[ID_L] , 3 ) ) ) ;

  tmp_Eh = ( ( 1 - kapreal ) * pC ) - ( kJ_t * y[ID_Eh] ) ;

  ydot[ID_Eh] = ( y[ID_Eh] < Ehp ? tmp_Eh : 0 ) ;

  tmp_R = kapR * ( ( ( 1 - kapreal ) * pC ) - ( kJ_t * Ehp ) ) ;

  ydot[ID_R] = ( y[ID_Eh] >= Ehp ? tmp_R / ( E_coc ) : 0 ) ;

  yout[ID_Energy] = y[ID_E] ;
  yout[ID_Maturity] = ( y[ID_Eh] <= 0 ? 1E-9 : y[ID_Eh] ) ;
  yout[ID_Reproduction] = ( y[ID_R] <= 0 ? 1E-9 : y[ID_R] ) ;
  yout[ID_Weight] = y[ID_W] ;
  yout[ID_Organic_matter] = y[ID_OM] ;

} /* derivs */


/*----- Jacobian calculations: */
void jac (int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{

} /* jac */


/*----- Events calculations: */
void event (int *n, double *t, double *y)
{

} /* event */

/*----- Roots calculations: */
void root (int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip)
{

} /* root */

