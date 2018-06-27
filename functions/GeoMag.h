//
//  GeoMag.h
//  a
//
//  Created by Nathan Zimmerberg on 6/24/18.
//  Copyright © 2018 Nathan Zimmerberg. All rights reserved.

//  Using  *      Revision Number: $Revision: 1288 $
//*      Last changed by: $Author: awoods $
//*      Last changed on: $Date: 2014-12-09 16:43:07 -0700 (Tue, 09 Dec 2014) $
/*
 Data types and prototype declaration for
 World Magnetic Model (WMM) subroutines.
 
 July 28, 2009
 
 manoj.c.nair@noaa.gov*/
//

#ifndef GeoMag_h
#define GeoMag_h

//mean radius of  ellipsoid in meters from section 1.2 of the WMM2015 Technical report
#define EARTH_R (6371200.0)



#define NMAX (12)
#define NMAXSECVAR (12)
#define MAG_NUMCOF ((NMAX+1)*(NMAX+2)/2)
#define MAG_NUMSECCOF ((NMAXSECVAR+1)*(NMAXSECVAR+2)/2)

/*These error values come from the ISCWSA error model:
 *http://www.copsegrove.com/Pages/MWDGeomagneticModels.aspx
 */
#define INCL_ERROR_BASE (0.20)
#define DECL_ERROR_OFFSET_BASE (0.36)
#define F_ERROR_BASE (130)
#define DECL_ERROR_SLOPE_BASE (5000)
#define WMM_ERROR_MULTIPLIER 1.21
#define IGRF_ERROR_MULTIPLIER 1.21

/*These error values are the NGDC error model
 *
 */
#define WMM_UNCERTAINTY_X 138
#define WMM_UNCERTAINTY_Y 89
#define WMM_UNCERTAINTY_Z 165




#define TRUE            ((int)1)
#define FALSE           ((int)0)

#define MAG_GEO_POLE_TOLERANCE  1e-6

/*
 Data types and prototype declaration for
 World Magnetic Model (WMM) subroutines.
 
 July 28, 2009
 
 manoj.c.nair@noaa.gov*/

typedef struct {
    const double epoch;
    const double Main_Field_Coeff_G[MAG_NUMCOF];
    const double Main_Field_Coeff_H[MAG_NUMCOF];
    const double Secular_Var_Coeff_G[MAG_NUMSECCOF];
    const double Secular_Var_Coeff_H[MAG_NUMSECCOF];
} MAGtype_ConstModel;



typedef struct {
    double EditionDate;
    double epoch; /*Base time of Geomagnetic model epoch (yrs)*/
    char ModelName[32];
    double *Main_Field_Coeff_G; /* C - Gauss coefficients of main geomagnetic model (nT) Index is (n * (n + 1) / 2 + m) */
    double *Main_Field_Coeff_H; /* C - Gauss coefficients of main geomagnetic model (nT) */
    double *Secular_Var_Coeff_G; /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */
    double *Secular_Var_Coeff_H; /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */
    int nMax; /* Maximum degree of spherical harmonic model */
    int nMaxSecVar; /* Maximum degree of spherical harmonic secular model */
    int SecularVariationUsed; /* Whether or not the magnetic secular variation vector will be needed by program*/
    double CoefficientFileEndDate;
    
} MAGtype_MagneticModel;



typedef struct {
    double lambda; /* longitude*/
    double phig; /* geocentric latitude*/
    double r; /* distance from the center of the ellipsoid*/
} MAGtype_CoordSpherical;

typedef struct {
    double x; /* x*/
    double y; /* y*/
    double z; /* z*/
} MAGtype_CoordECEF;

typedef struct {
    double *Pcup; /* Legendre Function */
    double *dPcup; /* Derivative of Legendre fcn */
} MAGtype_LegendreFunction;

typedef struct {
    double Bx; /* North */
    double By; /* East */
    double Bz; /* Down */
} MAGtype_MagneticResults;

typedef struct {
    double *RelativeRadiusPower; /* [earth_reference_radius_km / sph. radius ]^n  */
    double *cos_mlambda; /*cp(m)  - cosine of (m*spherical coord. longitude)*/
    double *sin_mlambda; /* sp(m)  - sine of (m*spherical coord. longitude) */
} MAGtype_SphericalHarmonicVariables;





int MAG_geomag(MAGtype_CoordECEF coordECEF, double dyear, MAGtype_MagneticResults *results);

/*Spherical Harmonics*/

int MAG_AssociatedLegendreFunction(MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_LegendreFunction *LegendreFunction);

int MAG_ComputeSphericalHarmonicVariables(double re, MAGtype_CoordSpherical CoordSpherical,
                                          int nMax,
                                          MAGtype_SphericalHarmonicVariables * SphVariables);

int MAG_Summation(MAGtype_LegendreFunction *LegendreFunction, MAGtype_MagneticModel *MagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults *MagneticResults);

int MAG_PcupLow(double *Pcup, double *dPcup, double x, int nMax);

int MAG_TimelyModifyMagneticModel(double dyear, MAGtype_ConstModel *MagneticModel, MAGtype_MagneticModel *TimedMagneticModel);

/*Memory and File Processing*/

MAGtype_LegendreFunction *MAG_AllocateLegendreFunctionMemory(int NumTerms);

MAGtype_MagneticModel *MAG_AllocateModelMemory(int NumTerms);

MAGtype_SphericalHarmonicVariables *MAG_AllocateSphVarMemory(int nMax);

int MAG_FreeMagneticModelMemory(MAGtype_MagneticModel *MagneticModel);

int MAG_FreeLegendreMemory(MAGtype_LegendreFunction *LegendreFunction);

int MAG_FreeSphVarMemory(MAGtype_SphericalHarmonicVariables *SphVar);

/*Coordinate transformations*/

MAGtype_MagneticResults MAG_ConvertNEDtoECEF(MAGtype_MagneticResults resultNED, MAGtype_CoordECEF coordECEF);

MAGtype_CoordSpherical MAG_ConvertECEFtoSph(MAGtype_CoordECEF coordECEF);










#endif /* GeoMag_h */
