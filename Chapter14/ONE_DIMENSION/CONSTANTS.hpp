/*
 * CONSTANTS.hpp
 *
 *  Created on: 12 May 2022
 *      Author: julio
 */

#ifndef CONSTANTS
#define CONSTANTS

// Entering concentration
double cae = 0.01;// gmol/cm3
// Entering temperature
double Tke = 305.0;//  K
// Initial concentration
double ca0 = 0.0;//   gmol/cm3
// Initial temperature
double Tk0 = 305.0;//   K
// Wall temperature
double Tw = 355.0;//   K
// reactor radius
double r0 = 1.0;//   cm
// Linear fluid velocity
double v = 1.0 ;//  cm/s
// Mass diffusivity
double Dc = 0.1 ;//  cm2/s
// Thermal diffusivity
double Dt = 0.1;//  = k/(ρC p) cm2/s
// Fluid density
double rho = 1.0 ;//  g/cm3
// Fluid specific heat
double Cp =  0.5 ;//  cal/(g·K)
// Heat of reaction
double dH = -10000.0; //  cal/gmol
// Specific rate constant
double rk0 = 1.5e9 ;//  cm3/(gmol·s)
// Activation energy
double E = 15000.0 ;//  cal/(gmol·K)
// Gas constant
double R = 1.987;//   cal/gmol·K
// Thermal conductivity
double k = 0.01;//   (cal·cm)/(s·cm2·K)
// Heat transfer coefficient
double h = 0.01;//   cal/(s·cm2·K)

double NTH=11;

double NR=11;

double ZL=100;

double D = 0.1;
double STD = 1.0;
double tau = 1.0;

#endif /* CONSTANTS_HPP_ */
