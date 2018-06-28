//
//  GeoMag.c
//  a
//
//  Modified by Nathan Zimmerberg on 6/24/18.
//  Modified from WMM2015 code from https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml
//Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, B. Hamilton, A. Woods, V. Ridley, S. Maus and A. Thomson, 2014. The US/UK //World Magnetic Model for 2015-2020, NOAA National Geophysical Data Center, Boulder, CO, doi: 10.7289/V5TH8JNW [6/24/18].
//Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, B. Hamilton, A. Woods, V. Ridley, S. Maus and A. Thomson, 2015. The US/UK World Magnetic Model for 2015-2020: Technical Report, NOAA National Geophysical Data Center, Boulder, CO, doi: 10.7289/V5TB14V7.
//
//Almost all of the code, and data is from WMM2015, this is not subject to copyright protection.
//Modifications are using ECEF coordinates, conversion from nT to gauss, use of radians instead of degrees, constant model struct instead of a WMM.COF file reader, and removal of some error checking, removal of special summation at poles, removal of high(>20) degree model support, and removal of geodetic stuff.


/* LICENSE of the WMM2015 source code
*
*  The WMM source code is in the public domain and not licensed or under copyright.
*    The information and software may be used freely by the public. As required by 17 U.S.C. 403,
*    third parties producing copyrighted works consisting predominantly of the material produced by
*    U.S. government agencies must provide notice with such work(s) identifying the U.S. Government material
*    incorporated and stating that such material is not subject to copyright protection.
*/

#include "GeoMag.h"
#include "GeoMagDATA.h"
#include <math.h>
//#include <stdio.h>
#include <stdlib.h>
//#include <ctype.h>




int MAG_geomag(MAGtype_CoordECEF coordECEF, double dyear, MAGtype_MagneticResults *results)
/*Return if success. Get the magnetic feild in ECEF coords units Gauss and puts the answer in *results
 The main subroutine that calls a sequence of WMM sub-functions to calculate the magnetic field elements for a single point.
 The function expects the point coordinates and dyear as input and returns the magnetic field elements
 The function uses constants from GeoMag.h and GeoMagDATA.h
 
 INPUT:
    coordECEF(Above the surface of earth, Not at the geographic poles): The location where the field is predicted, units m
    dyear(should be around 2015.0): The decimal year, years far from 2015.0 are more innaccurate, units years CE
    results(): Pointer to where to store the results.
 
 
 */
{
    MAGtype_LegendreFunction *LegendreFunction;
    MAGtype_SphericalHarmonicVariables *SphVariables;
    MAGtype_MagneticModel  *TimedMagneticModel;
    MAGtype_MagneticResults MagneticResultsSph;
    MAGtype_CoordSpherical CoordSpherical;


    
    //allocate memory
    LegendreFunction = MAG_AllocateLegendreFunctionMemory(MAG_NUMCOF); /* For storing the ALF functions */
    SphVariables = MAG_AllocateSphVarMemory(NMAX);
    TimedMagneticModel = MAG_AllocateModelMemory(MAG_NUMCOF); /* For storing the time modified WMM Model parameters */
    //check if mem allocation failed here
    if (!(LegendreFunction && SphVariables && TimedMagneticModel)){
        MAG_FreeLegendreMemory(LegendreFunction);
        MAG_FreeSphVarMemory(SphVariables);
        MAG_FreeMagneticModelMemory(TimedMagneticModel);
        return FALSE;
    }
    //check to see if it is at the earth center or poles here.
    if (fabs(coordECEF.x)<=MAG_GEO_POLE_TOLERANCE && fabs(coordECEF.y)<=MAG_GEO_POLE_TOLERANCE){
        MAG_FreeLegendreMemory(LegendreFunction);
        MAG_FreeSphVarMemory(SphVariables);
        MAG_FreeMagneticModelMemory(TimedMagneticModel);
        return FALSE;
    }
    CoordSpherical= MAG_ConvertECEFtoSph(coordECEF);
    MAG_TimelyModifyMagneticModel(dyear, &WMM2015, TimedMagneticModel); /* Time adjust the coefficients, Equation 9, WMM Technical report */
    
    MAG_ComputeSphericalHarmonicVariables(EARTH_R, CoordSpherical, NMAX, SphVariables); /* Compute Spherical Harmonic variables  */
    MAG_AssociatedLegendreFunction(CoordSpherical, NMAX, LegendreFunction); /* Compute ALF  */
    MAG_Summation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical, &MagneticResultsSph); /* Accumulate the spherical harmonic coefficients*/
    *results= MAG_ConvertNEDtoECEF(MagneticResultsSph, coordECEF);
    //converion to gauss
    results->Bx*=1.0e-5;
    results->By*=1.0e-5;
    results->Bz*=1.0e-5;
    

    MAG_FreeLegendreMemory(LegendreFunction);
    MAG_FreeSphVarMemory(SphVariables);
    MAG_FreeMagneticModelMemory(TimedMagneticModel);
    
    return TRUE;
} /*MAG_Geomag*/





int MAG_AssociatedLegendreFunction(MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_LegendreFunction *LegendreFunction)

/* Computes  all of the Schmidt-semi normalized associated Legendre
 functions up to degree nMax. If nMax <= 16, function MAG_PcupLow is used.
 Otherwise MAG_PcupHigh is called.
 INPUT  CoordSpherical     A data structure with the following elements
        double lambda; ( longitude)
        double phig; ( geocentric latitude )
        double r;        ( distance from the center of the ellipsoid)
    nMax            (positive integer must be <= 16)      Maxumum degree of spherical harmonic secular model.
    LegendreFunction Pointer to data structure with the following elements
        double *Pcup;  (  pointer to store Legendre Function  )
        double *dPcup; ( pointer to store  Derivative of Lagendre function )
 
 OUTPUT  LegendreFunction  Calculated Legendre variables in the data structure
 
 */
{
    double sin_phi;
    int FLAG = 1;
    
    sin_phi = sin(CoordSpherical.phig); /* sin  (geocentric latitude) */
    FLAG = MAG_PcupLow(LegendreFunction->Pcup, LegendreFunction->dPcup, sin_phi, nMax);
    if(FLAG == 0) /* Error while computing  Legendre variables*/
        return FALSE;
    
    
    return TRUE;
} /*MAG_AssociatedLegendreFunction */

int MAG_ComputeSphericalHarmonicVariables(double re, MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_SphericalHarmonicVariables *SphVariables)

/* Computes Spherical variables
 Variables computed are (a/r)^(n+2), cos_m(lamda) and sin_m(lambda) for spherical harmonic
 summations. (Equations 10-12 in the WMM Technical Report)
 INPUT
    re( greater than zero); mean radius of  ellipsoid
    CoordSpherical     A data structure with the following elements
        double lambda; ( longitude) units rad
        double phig; ( geocentric latitude ) units rad
        double r;        ( distance from the center of the ellipsoid greater than 0.0) units m
    nMax   integer      ( Maxumum degree of spherical harmonic secular model)\
 
 OUTPUT  SphVariables  Pointer to the   data structure with the following elements
        double RelativeRadiusPower[MAG_MAX_MODEL_DEGREES+1];   [earth_reference_radius_km  sph. radius ]^n
        double cos_mlambda[MAG_MAX_MODEL_DEGREES+1]; cp(m)  - cosine of (mspherical coord. longitude)
        double sin_mlambda[MAG_MAX_MODEL_DEGREES+1];  sp(m)  - sine of (mspherical coord. longitude)
 CALLS : none
 */
{
    double cos_lambda, sin_lambda;
    int m, n;
    cos_lambda = cos(CoordSpherical.lambda);
    sin_lambda = sin(CoordSpherical.lambda);
    /* for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2)
     for n  1..nMax-1 (this is much faster than calling pow MAX_N+1 times).      */
    SphVariables->RelativeRadiusPower[0] = (re / CoordSpherical.r) * (re / CoordSpherical.r);
    for(n = 1; n <= nMax; n++)
    {
        SphVariables->RelativeRadiusPower[n] = SphVariables->RelativeRadiusPower[n - 1] * (re / CoordSpherical.r);
    }
    
    /*
     Compute cos(m*lambda), sin(m*lambda) for m = 0 ... nMax
     cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
     sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
     */
    SphVariables->cos_mlambda[0] = 1.0;
    SphVariables->sin_mlambda[0] = 0.0;
    
    SphVariables->cos_mlambda[1] = cos_lambda;
    SphVariables->sin_mlambda[1] = sin_lambda;
    for(m = 2; m <= nMax; m++)
    {
        SphVariables->cos_mlambda[m] = SphVariables->cos_mlambda[m - 1] * cos_lambda - SphVariables->sin_mlambda[m - 1] * sin_lambda;
        SphVariables->sin_mlambda[m] = SphVariables->cos_mlambda[m - 1] * sin_lambda + SphVariables->sin_mlambda[m - 1] * cos_lambda;
    }
    return TRUE;
} /*MAG_ComputeSphericalHarmonicVariables*/

int MAG_Summation(MAGtype_LegendreFunction *LegendreFunction, MAGtype_MagneticModel *MagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults *MagneticResults)
{
    /* Computes Geomagnetic Field Elements X, Y and Z in Spherical coordinate system using
     spherical harmonic summation.
     
     
     The vector Magnetic field is given by -grad V, where V is Geomagnetic scalar potential
     The gradient in spherical coordinates is given by:
     
     dV ^     1 dV ^        1     dV ^
     grad V = -- r  +  - -- t  +  -------- -- p
     dr       r dt       r sin(t) dp
     
     
     INPUT :  LegendreFunction
        MagneticModel
        SphVariables
        CoordSpherical(not be at center of earth or at the poles)
     OUTPUT : MagneticResults
     
     CALLS : MAG_SummationSpecial
     
     
     
     Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
     */
    int m, n, index;
    double cos_phi;
    MagneticResults->Bz = 0.0;
    MagneticResults->By = 0.0;
    MagneticResults->Bx = 0.0;
    for(n = 1; n <= MagneticModel->nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            
            /*            nMax      (n+2)       n     m            m           m
             Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
             n=1                m=0   n            n           n  */
            /* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
            MagneticResults->Bz -= SphVariables.RelativeRadiusPower[n] *
            (MagneticModel->Main_Field_Coeff_G[index] * SphVariables.cos_mlambda[m] +
             MagneticModel->Main_Field_Coeff_H[index] * SphVariables.sin_mlambda[m])
            * (double) (n + 1) * LegendreFunction-> Pcup[index];
            
            /*          1 nMax  (n+2)    n     m            m           m
             By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
             n=1             m=0   n            n           n  */
            /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
            MagneticResults->By += SphVariables.RelativeRadiusPower[n] *
            (MagneticModel->Main_Field_Coeff_G[index] * SphVariables.sin_mlambda[m] -
             MagneticModel->Main_Field_Coeff_H[index] * SphVariables.cos_mlambda[m])
            * (double) (m) * LegendreFunction-> Pcup[index];
            /*           nMax  (n+2) n     m            m           m
             Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
             n=1         m=0   n            n           n  */
            /* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */
            
            MagneticResults->Bx -= SphVariables.RelativeRadiusPower[n] *
            (MagneticModel->Main_Field_Coeff_G[index] * SphVariables.cos_mlambda[m] +
             MagneticModel->Main_Field_Coeff_H[index] * SphVariables.sin_mlambda[m])
            * LegendreFunction-> dPcup[index];
            
            
            
        }
    }
    
    cos_phi = cos(CoordSpherical.phig);
    MagneticResults->By = MagneticResults->By / cos_phi;
    /* Special calculation for component - By - at Geographic poles.
     * If the user wants to avoid using this function,  please make sure that
     * the latitude is not exactly +/-pi/2. An option is to make use the function
     * MAG_CheckGeographicPoles.
     */
    return TRUE;
}/*MAG_Summation */

int MAG_PcupLow(double *Pcup, double *dPcup, double x, int nMax)

/*   This function evaluates all of the Schmidt-semi normalized associated Legendre
 functions up to degree nMax.
 
 Calling Parameters:
 INPUT
 nMax:     Maximum spherical harmonic degree to compute.
 x:        cos(colatitude) or sin(latitude).
 
 OUTPUT
 Pcup:    A vector of all associated Legendgre polynomials evaluated at
 x up to nMax.
 dPcup: Derivative of Pcup(x) with respect to latitude
 
 Notes: Overflow may occur if nMax > 20 , especially for high-latitudes.
 Use MAG_PcupHigh for large nMax.
 
 Written by Manoj Nair, June, 2009 . Manoj.C.Nair@Noaa.Gov.
 
 Note: In geomagnetism, the derivatives of ALF are usually found with
 respect to the colatitudes. Here the derivatives are found with respect
 to the latitude. The difference is a sign reversal for the derivative of
 the Associated Legendre Functions.
 */
{
    int n, m, index, index1, index2, NumTerms;
    double k, z, *schmidtQuasiNorm;
    Pcup[0] = 1.0;
    dPcup[0] = 0.0;
    /*sin (geocentric latitude) - sin_phi */
    z = sqrt((1.0 - x) * (1.0 + x));
    
    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    schmidtQuasiNorm = (double *) malloc((NumTerms + 1) * sizeof ( double));
    
    if(schmidtQuasiNorm == NULL)
    {
        //MAG_Error(19);
        return FALSE;
    }
    
    /*     First,    Compute the Gauss-normalized associated Legendre  functions*/
    for(n = 1; n <= nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if(n == m)
            {
                index1 = (n - 1) * n / 2 + m - 1;
                Pcup [index] = z * Pcup[index1];
                dPcup[index] = z * dPcup[index1] + x * Pcup[index1];
            } else if(n == 1 && m == 0)
            {
                index1 = (n - 1) * n / 2 + m;
                Pcup[index] = x * Pcup[index1];
                dPcup[index] = x * dPcup[index1] - z * Pcup[index1];
            } else if(n > 1 && n != m)
            {
                index1 = (n - 2) * (n - 1) / 2 + m;
                index2 = (n - 1) * n / 2 + m;
                if(m > n - 2)
                {
                    Pcup[index] = x * Pcup[index2];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2];
                } else
                {
                    k = (double) (((n - 1) * (n - 1)) - (m * m)) / (double) ((2 * n - 1) * (2 * n - 3));
                    Pcup[index] = x * Pcup[index2] - k * Pcup[index1];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];
                }
            }
        }
    }
    /* Compute the ration between the the Schmidt quasi-normalized associated Legendre
     * functions and the Gauss-normalized version. */
    
    schmidtQuasiNorm[0] = 1.0;
    for(n = 1; n <= nMax; n++)
    {
        index = (n * (n + 1) / 2);
        index1 = (n - 1) * n / 2;
        /* for m = 0 */
        schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (double) (2 * n - 1) / (double) n;
        
        for(m = 1; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            index1 = (n * (n + 1) / 2 + m - 1);
            schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * sqrt((double) ((n - m + 1) * (m == 1 ? 2 : 1)) / (double) (n + m));
        }
        
    }
    
    /* Converts the  Gauss-normalized associated Legendre
     functions to the Schmidt quasi-normalized version using pre-computed
     relation stored in the variable schmidtQuasiNorm */
    
    for(n = 1; n <= nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            Pcup[index] = Pcup[index] * schmidtQuasiNorm[index];
            dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index];
            /* The sign is changed since the new WMM routines use derivative with respect to latitude
             insted of co-latitude */
        }
    }
    
    if(schmidtQuasiNorm)
        free(schmidtQuasiNorm);
    return TRUE;
} /*MAG_PcupLow */

int MAG_TimelyModifyMagneticModel(double dyear, MAGtype_ConstModel *MagneticModel, MAGtype_MagneticModel *TimedMagneticModel)

/* Time change the Model coefficients from the base year of the model using secular variation coefficients.
 Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
 Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.
 If the address of "TimedMagneticModel" is the same as the address of "MagneticModel", then this procedure overwrites
 the given item "MagneticModel".
 
 INPUT: dyear
 MagneticModel
 OUTPUT:TimedMagneticModel
 CALLS : none
 */
{
    int n, m, index, a, b;
    //TimedMagneticModel->EditionDate = MagneticModel->EditionDate;
    TimedMagneticModel->epoch = (double) MagneticModel->epoch;
    TimedMagneticModel->nMax = NMAX;
    TimedMagneticModel->nMaxSecVar = NMAXSECVAR;
    a = TimedMagneticModel->nMaxSecVar;
    b = (a * (a + 1) / 2 + a);
    //strcpy(TimedMagneticModel->ModelName, MagneticModel->ModelName);
    for(n = 1; n <= NMAX; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if(index <= b)
            {
                TimedMagneticModel->Main_Field_Coeff_H[index] = (double)MagneticModel->Main_Field_Coeff_H[index] + (dyear -  MagneticModel->epoch) * (double)MagneticModel->Secular_Var_Coeff_H[index];
                TimedMagneticModel->Main_Field_Coeff_G[index] = (double)MagneticModel->Main_Field_Coeff_G[index] + (dyear - MagneticModel->epoch) * (double)MagneticModel->Secular_Var_Coeff_G[index];
                TimedMagneticModel->Secular_Var_Coeff_H[index] = (double)MagneticModel->Secular_Var_Coeff_H[index]; /* We need a copy of the secular var coef to calculate secular change */
                TimedMagneticModel->Secular_Var_Coeff_G[index] = (double)MagneticModel->Secular_Var_Coeff_G[index];
            } else
            {
                TimedMagneticModel->Main_Field_Coeff_H[index] = (double)MagneticModel->Main_Field_Coeff_H[index];
                TimedMagneticModel->Main_Field_Coeff_G[index] = (double)MagneticModel->Main_Field_Coeff_G[index];
            }
        }
    }
    return TRUE;
} /* MAG_TimelyModifyMagneticModel */


/*Memory management*/

MAGtype_LegendreFunction *MAG_AllocateLegendreFunctionMemory(int NumTerms)

/* Allocate memory for Associated Legendre Function data types.
 Should be called before computing Associated Legendre Functions.
 
 INPUT: NumTerms : int : Total number of spherical harmonic coefficients in the model
 
 
 OUTPUT:    Pointer to data structure MAGtype_LegendreFunction with the following elements
 double *Pcup;  (  pointer to store Legendre Function  )
 double *dPcup; ( pointer to store  Derivative of Legendre function )
 
 FALSE: Failed to allocate memory
 
 CALLS : none
 
 */
{
    MAGtype_LegendreFunction *LegendreFunction;
    
    LegendreFunction = (MAGtype_LegendreFunction *) calloc(1, sizeof (MAGtype_LegendreFunction));
    
    if(!LegendreFunction)
    {
        //MAG_Error(1);
        return NULL;
    }
    LegendreFunction->Pcup = (double *) malloc((NumTerms + 1) * sizeof ( double));
    if(LegendreFunction->Pcup == 0)
    {
        //MAG_Error(1);
        return NULL;
    }
    LegendreFunction->dPcup = (double *) malloc((NumTerms + 1) * sizeof ( double));
    if(LegendreFunction->dPcup == 0)
    {
        //MAG_Error(1);
        return NULL;
    }
    return LegendreFunction;
} /*MAGtype_LegendreFunction*/

MAGtype_MagneticModel *MAG_AllocateModelMemory(int NumTerms)

/* Allocate memory for WMM Coefficients
 * Should be called before reading the model file *
 
 INPUT: NumTerms : int : Total number of spherical harmonic coefficients in the model
 
 
 OUTPUT:    Pointer to data structure MAGtype_MagneticModel with the following elements
 double EditionDate;
 double epoch;       Base time of Geomagnetic model epoch (yrs)
 char  ModelName[20];
 double *Main_Field_Coeff_G;          C - Gauss coefficients of main geomagnetic model (nT)
 double *Main_Field_Coeff_H;          C - Gauss coefficients of main geomagnetic model (nT)
 double *Secular_Var_Coeff_G;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
 double *Secular_Var_Coeff_H;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
 int nMax;  Maximum degree of spherical harmonic model
 int nMaxSecVar; Maxumum degree of spherical harmonic secular model
 int SecularVariationUsed; Whether or not the magnetic secular variation vector will be needed by program
 
 FALSE: Failed to allocate memory
 CALLS : none
 */
{
    MAGtype_MagneticModel *MagneticModel;
    int i;
    
    
    MagneticModel = (MAGtype_MagneticModel *) calloc(1, sizeof (MAGtype_MagneticModel));
    
    if(MagneticModel == NULL)
    {
        //MAG_Error(2);
        return NULL;
    }
    
    MagneticModel->Main_Field_Coeff_G = (double *) malloc((NumTerms + 1) * sizeof ( double));
    
    if(MagneticModel->Main_Field_Coeff_G == NULL)
    {
        //MAG_Error(2);
        return NULL;
    }
    
    MagneticModel->Main_Field_Coeff_H = (double *) malloc((NumTerms + 1) * sizeof ( double));
    
    if(MagneticModel->Main_Field_Coeff_H == NULL)
    {
        //MAG_Error(2);
        return NULL;
    }
    MagneticModel->Secular_Var_Coeff_G = (double *) malloc((NumTerms + 1) * sizeof ( double));
    if(MagneticModel->Secular_Var_Coeff_G == NULL)
    {
        //MAG_Error(2);
        return NULL;
    }
    MagneticModel->Secular_Var_Coeff_H = (double *) malloc((NumTerms + 1) * sizeof ( double));
    if(MagneticModel->Secular_Var_Coeff_H == NULL)
    {
        //MAG_Error(2);
        return NULL;
    }
    MagneticModel->CoefficientFileEndDate = 0;
    MagneticModel->EditionDate = 0;
    //strcpy(MagneticModel->ModelName, "");
    MagneticModel->SecularVariationUsed = 0;
    MagneticModel->epoch = 0;
    MagneticModel->nMax = 0;
    MagneticModel->nMaxSecVar = 0;
    
    for(i=0; i<NumTerms; i++) {
        MagneticModel->Main_Field_Coeff_G[i] = 0;
        MagneticModel->Main_Field_Coeff_H[i] = 0;
        MagneticModel->Secular_Var_Coeff_G[i] = 0;
        MagneticModel->Secular_Var_Coeff_H[i] = 0;
    }
    
    return MagneticModel;
    
} /*MAG_AllocateModelMemory*/

MAGtype_SphericalHarmonicVariables* MAG_AllocateSphVarMemory(int nMax)
{
    MAGtype_SphericalHarmonicVariables* SphVariables;
    SphVariables  = (MAGtype_SphericalHarmonicVariables*) calloc(1, sizeof(MAGtype_SphericalHarmonicVariables));
    SphVariables->RelativeRadiusPower = (double *) malloc((nMax + 1) * sizeof ( double));
    SphVariables->cos_mlambda = (double *) malloc((nMax + 1) * sizeof (double));
    SphVariables->sin_mlambda = (double *) malloc((nMax + 1) * sizeof (double));
    return SphVariables;
} /*MAG_AllocateSphVarMemory*/

int MAG_FreeMagneticModelMemory(MAGtype_MagneticModel *MagneticModel)

/* Free the magnetic model memory used by WMM functions.
 INPUT :  MagneticModel    pointer to data structure with the following elements
 
 double EditionDate;
 double epoch;       Base time of Geomagnetic model epoch (yrs)
 char  ModelName[20];
 double *Main_Field_Coeff_G;          C - Gauss coefficients of main geomagnetic model (nT)
 double *Main_Field_Coeff_H;          C - Gauss coefficients of main geomagnetic model (nT)
 double *Secular_Var_Coeff_G;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
 double *Secular_Var_Coeff_H;  CD - Gauss coefficients of secular geomagnetic model (nT/yr)
 int nMax;  Maximum degree of spherical harmonic model
 int nMaxSecVar; Maxumum degree of spherical harmonic secular model
 int SecularVariationUsed; Whether or not the magnetic secular variation vector will be needed by program
 
 OUTPUT  none
 CALLS : none
 
 */
{
    if(MagneticModel->Main_Field_Coeff_G)
    {
        free(MagneticModel->Main_Field_Coeff_G);
        MagneticModel->Main_Field_Coeff_G = NULL;
    }
    if(MagneticModel->Main_Field_Coeff_H)
    {
        free(MagneticModel->Main_Field_Coeff_H);
        MagneticModel->Main_Field_Coeff_H = NULL;
    }
    if(MagneticModel->Secular_Var_Coeff_G)
    {
        free(MagneticModel->Secular_Var_Coeff_G);
        MagneticModel->Secular_Var_Coeff_G = NULL;
    }
    if(MagneticModel->Secular_Var_Coeff_H)
    {
        free(MagneticModel->Secular_Var_Coeff_H);
        MagneticModel->Secular_Var_Coeff_H = NULL;
    }
    if(MagneticModel)
    {
        free(MagneticModel);
        MagneticModel = NULL;
    }
    
    return TRUE;
} /*MAG_FreeMagneticModelMemory */

int MAG_FreeLegendreMemory(MAGtype_LegendreFunction *LegendreFunction)

/* Free the Legendre Coefficients memory used by the WMM functions.
 INPUT : LegendreFunction Pointer to data structure with the following elements
 double *Pcup;  (  pointer to store Legendre Function  )
 double *dPcup; ( pointer to store  Derivative of Lagendre function )
 
 OUTPUT: none
 CALLS : none
 
 */
{
    if(LegendreFunction->Pcup)
    {
        free(LegendreFunction->Pcup);
        LegendreFunction->Pcup = NULL;
    }
    if(LegendreFunction->dPcup)
    {
        free(LegendreFunction->dPcup);
        LegendreFunction->dPcup = NULL;
    }
    if(LegendreFunction)
    {
        free(LegendreFunction);
        LegendreFunction = NULL;
    }
    
    return TRUE;
} /*MAG_FreeLegendreMemory */

int MAG_FreeSphVarMemory(MAGtype_SphericalHarmonicVariables *SphVar)

/* Free the Spherical Harmonic Variable memory used by the WMM functions.
 INPUT : LegendreFunction Pointer to data structure with the following elements
 double *RelativeRadiusPower
 double *cos_mlambda
 double *sin_mlambda
 OUTPUT: none
 CALLS : none
 */
{
    if(SphVar->RelativeRadiusPower)
    {
        free(SphVar->RelativeRadiusPower);
        SphVar->RelativeRadiusPower = NULL;
    }
    if(SphVar->cos_mlambda)
    {
        free(SphVar->cos_mlambda);
        SphVar->cos_mlambda = NULL;
    }
    if(SphVar->sin_mlambda)
    {
        free(SphVar->sin_mlambda);
        SphVar->sin_mlambda = NULL;
    }
    if(SphVar)
    {
        free(SphVar);
        SphVar = NULL;
    }
    
    return TRUE;
} /*MAG_FreeSphVarMemory*/






/*coordinate transforms*/

MAGtype_MagneticResults MAG_ConvertNEDtoECEF(MAGtype_MagneticResults resultNED, MAGtype_CoordECEF coordECEF){
    /* Return the resultNED vector in ECEF coordinates units nT.
     
     
     INPUT:
        resultNED(): The Magnetic feild in NED(north east down) coords units nT,
        coordECEF(not at the poles or center of earth): The location of the results in ECEF coords units m
    */
    
    double r;
    double rxy;
    MAGtype_MagneticResults resultECEF;
    r= sqrt(coordECEF.x*coordECEF.x+coordECEF.y*coordECEF.y+coordECEF.z*coordECEF.z);
    rxy= sqrt(coordECEF.x*coordECEF.x+coordECEF.y*coordECEF.y);
    resultECEF.Bx= -resultNED.Bx*coordECEF.x*coordECEF.z/r/rxy - resultNED.By*coordECEF.y/rxy - resultNED.Bz*coordECEF.x/r ;
    resultECEF.By= -resultNED.Bx*coordECEF.y*coordECEF.z/r/rxy + resultNED.By*coordECEF.x/rxy - resultNED.Bz*coordECEF.y/r ;
    resultECEF.Bz= resultNED.Bx*(coordECEF.x*coordECEF.x+coordECEF.y*coordECEF.y)/r/rxy - resultNED.Bz*coordECEF.z/r ;
    return resultECEF;
}


MAGtype_CoordSpherical MAG_ConvertECEFtoSph(MAGtype_CoordECEF coordECEF){
    /* Return the Spherical coordinants from ECEF coordinates units m and rad.
        If at the geographic poles, lambda is 0.0
     
     INPUT:
        coordECEF(radius cannot be 0): The location of the results in ECEF coords units m
     */
    
    MAGtype_CoordSpherical sph;
    sph.r= sqrt(coordECEF.x*coordECEF.x+coordECEF.y*coordECEF.y+coordECEF.z*coordECEF.z);
    sph.phig= asin(coordECEF.z/sph.r);
    sph.lambda= atan2(coordECEF.y, coordECEF.x);
    return sph;
}






