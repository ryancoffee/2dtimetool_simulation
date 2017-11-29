#ifndef REFRACTION_H
#define REFRACTION_H
#include <cmath>
#include <algorithm>
#include "Constants.hpp"

#define SELLMEIER_SAPPHIRE(lamda) ( sqrt( 1.0 + ( 1.023798 * gsl_pow_2(lambda/1000.0)/( gsl_pow_2(lambda/1000.0) -  0.00377588) ) + (1.058264 * gsl_pow_2(lambda/1000.0)/( gsl_pow_2(lambda/1000.0) - 0.0122544  ) ) + ( 5.2807922 *gsl_pow_2(lambda/1000.0)/( gsl_pow_2(lambda/1000.0) - 321.3616  ) )  ) )

#define SAPPHIRE_fs(length_mm,lamda,lamda_0) ( ( SELLMEIER_SAPPHIRE(lambda) - SELLMEIER_SAPPHIRE(lambda_0) ) * length_mm / C_mmPfs )

namespace Refraction{
	template <typename T>
		inline T& sellmeier_sapphire(T& lamda){
			std::vector<T> coeffs = {1.023798,0.00377588,1.058264,0.0122544,5.2807922,321.3616};
			T val = std::sqrt( 
					T(1.0) 
					+ ( coeffs[0] * std::pow(lambda/T(1000),int(2))/( std::pow(lambda/T(1000),int(2)) -  coeffs[1]) ) 
					+ ( coeffs[2] * std::pow(lambda/1000.0,int(2))/( std::pow(lambda/1000.0,int(2)) - coeffs[3]) ) 
					+ ( coeffs[4] * std::pow(lambda/1000.0,int(2))/( std::pow(lambda/1000.0,int(2)) - coeffs[5]) )   
					);
		}
	template <typename T>
		inline T& sapphire_fs(T& length_mm,T& lamda,T& lamda_0) { 
			T val = (sellmeier_sapphire(lambda) - sellmeier_sapphire(lambda_0) ) * length_mm / C_mmPfs<T>();
			return val;
		}

}

// now add a good way to phase the complex vectors properly.
#endif
