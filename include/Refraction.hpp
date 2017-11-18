#ifndef REFRACTION_H
#define REFRACTION_H
#include <math.h>
#include "Constants.hpp"

#define SELLMEIER_SAPPHIRE(lamda) ( sqrt( 1.0 + ( 1.023798 * gsl_pow_2(lambda/1000.0)/( gsl_pow_2(lambda/1000.0) -  0.00377588) ) + (1.058264 * gsl_pow_2(lambda/1000.0)/( gsl_pow_2(lambda/1000.0) - 0.0122544  ) ) + ( 5.2807922 *gsl_pow_2(lambda/1000.0)/( gsl_pow_2(lambda/1000.0) - 321.3616  ) )  ) )

#define SAPPHIRE_fs(length_mm,lamda,lamda_0) ( ( SELLMEIER_SAPPHIRE(lambda) - SELLMEIER_SAPPHIRE(lambda_0) ) * length_mm / C_mmPfs )
// now add a good way to phase the complex vectors properly.
#endif
