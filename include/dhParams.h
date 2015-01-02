//**************************************************
//This file contains the dh-parameters of the WAM
#ifndef DHPARAMS_H_
#define DHPARAMS_H_

#include <cmath>

#define PI 3.14
const float a[7]={0, 0, 0.045, -0.045, 0, 0, 0};
const float alpha[7]={-PI/2, PI/2, -PI/2, PI/2, -PI/2, PI/2, 0};
const float d[7]={0, 0, 0.55, 0, 0.3, 0, 0.06};
  // in radians

//extern float theta[7];
// Mass parameterd
/*
extern float d11= 3.01* ((cos(theta[1])* (0.0416* sin(theta[3])
				   + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
				   (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))*(cos(theta[1])* (0.0416* sin(theta[3])
					+ 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
					(0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
		     + (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
			   + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))*(cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
					   + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)))

    + 5.677* (((0.119* sin(theta[1]) + 0.0005* cos(theta[1])))*(cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
			   + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))

	      + ((0.119* cos(theta[1]) - 0.0005* sin(theta[1])))*((0.119* cos(theta[1]) - 0.0005* sin(theta[1]))) )
    + 0.174;

extern  float d12=3.01*((cos(theta[1])* (0.041* sin(theta[3])
			      + 0.171* cos(theta[3])) + sin(theta[1])* (0.041* cos(theta[3])
							    - 0.171* sin(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
											   + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
										 (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
		   + (sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
		      - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
		   * (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
		      + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))) + 0.093254778214;


   extern float d21=3.010* ((cos(theta[1])* (0.041* sin(theta[3])
				+ 0.171* cos(theta[3])) + sin(theta[1])* (0.041* cos(theta[3])
							      - 0.171* sin(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
											     + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
										   (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
		     + (sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
			- cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*
		     (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
		      + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))) + 0.093254778214;

extern  float d22=3.01* ((sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
		      - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*(sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
				      - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
		  + (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
			+ sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*(cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
					+ sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))))+0.0932547782144;
*/

#endif /*DHPARAMS_H*/
