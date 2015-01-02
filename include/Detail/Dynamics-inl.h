/*
 * Dynamics.cpp
 *
 *  Created on: 25-Aug-2012
 *      Author: mobman
 */

//#include <Dynamics.h>
#include <dhParams.h>
#include <typeinfo>
#include <math.h>
#include <sstream>

namespace Sam {

template<typename InputType1, typename InputType2>
Dynamics<InputType1, InputType2>::Dynamics():massMatrix(Eigen::Matrix2d()), Cmatrix(Eigen::Matrix2d()), phi(Eigen::Vector2d()), massMatrix4(Eigen::Matrix4d()), fx(Eigen::Matrix4d())
					,Cmatrix4(Eigen::Matrix4d()), phi4(Eigen::Vector4d()),gx(Eigen::Matrix<double, 4,2>()) {
	phi.Zero();
	massMatrix4.fill(0);
	Cmatrix4.fill(0);
	phi4.fill(0);
	gx.fill(0);

}

// Mass matrix definition/////////////////////////////////////////
template<typename InputType1, typename InputType2>
void Dynamics<InputType1, InputType2>::makemassMat( const InputType1& theta){

	double  d11=3.01* ((cos(theta[1])* (0.0416* sin(theta[3])
	   	   + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
	   	   (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))*(cos(theta[1])* (0.0416* sin(theta[3])
			+ 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
			(0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
			+ (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
					+ sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))*(cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
							+ sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)))

							+ 5.677* (((0.119* sin(theta[1]) + 0.0005* cos(theta[1])))*((0.119* sin(theta[1]) + 0.0005* cos(theta[1])))

									+ ((0.119* cos(theta[1]) - 0.0005* sin(theta[1])))*((0.119* cos(theta[1]) - 0.0005* sin(theta[1]))) )
									+ 0.174;

	  double d12=3.01*((cos(theta[1])* (0.041* sin(theta[3])
				      + 0.171* cos(theta[3])) + sin(theta[1])* (0.041* cos(theta[3])
								    - 0.171* sin(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
												   + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
											 (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
			   + (sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
			      - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
			   * (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
			      + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))) + 0.093254778214;


	  double d21=3.010* ((cos(theta[1])* (0.041* sin(theta[3])
					+ 0.171* cos(theta[3])) + sin(theta[1])* (0.041* cos(theta[3])
								      - 0.171* sin(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
												     + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
											   (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
			     + (sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
				- cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*
			     (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
			      + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))) + 0.093254778214;

	  double d22=3.01* ((sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
			      - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*(sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
					      - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
			  + (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
				+ sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*(cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
						+ sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))))+0.0932547782144;


	massMatrix<<d11, d12, d21, d22;
}

// End of mass matrix////////////////////////////////////////////////////

// Definition of C matrix//////////////////////////////////////////////////////////////
template<typename InputType1, typename InputType2>
void Dynamics<InputType1, InputType2>::makeCmat(const InputType1& theta, const InputType2& thetaDot){

	double  c111=0.5* (3.011* (2* (cos(theta[1])* (0.171* sin(theta[3])
	 - 0.041* cos(theta[3]) + 0.045) + sin(theta[1])*
	 (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))*
	 (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)
	 - sin(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
	 + 2* (- cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
	 - sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))
				 * (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)
	 - sin(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)))
	 + 5.677* (2* (0.119* cos(theta[1]) - 0.0005* sin(theta[1]))*
	 (0.119* sin(theta[1]) + 0.0005* cos(theta[1]))
	 + 2* (- 0.119* sin(theta[1]) - 0.0005* cos(theta[1]))*
		   (0.119* cos(theta[1]) - 0.0005* sin(theta[1]))));

	 double c211= 1.505* (2* (cos(theta[1])* (0.041* cos(theta[3])
	 - 0.171* sin(theta[3])) - sin(theta[1])* (0.041* sin(theta[3])
	 + 0.171* cos(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
	 + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
	 (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
	 + 2* (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*
	 (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
	  + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)));

	 double c112=3.011* ((sin(theta[1])* (0.041* sin(theta[3])
	 + 0.171* cos(theta[3])) - cos(theta[1])* (0.041* cos(theta[3])
	 - 0.171* sin(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
	 + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
	 (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
	 + (cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))
	 - sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3])))
	* (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)
	 - sin(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
	 + (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	* (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
	 + sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))
	 + (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	* (- cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
	 - sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)))
	 - 1.505* (2* (cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))
	 - sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3])))
	 *(cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55)
	 - sin(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
	 + 2 *(cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	* (cos(theta[1])* (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
	   + sin(theta[1])* (0.041* sin(theta[3]) + 0.171*cos(theta[3]) + 0.55)));

	 double c212= 1.505*( 2* (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	 *(sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	 + 2* (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	* (cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))
	   - sin(theta[1]) *(0.0410*sin(theta[3]) + 0.171* cos(theta[3]))));

	 double c121=c211;
	 double c221=3.011* ((sin(theta[1])* (- 0.041* sin(theta[3])
	 - 0.171* cos(theta[3])) + cos(theta[1]) *(0.041* cos(theta[3])
	 - 0.171* sin(theta[3])))* (cos(theta[1])* (0.041* sin(theta[3])
	 + 0.171* cos(theta[3]) + 0.55) - sin(theta[1])*
	 (0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045))
	 + (sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))
	 - cos(theta[1])* (- 0.041* sin(theta[3]) - 0.171* cos(theta[3])))
	* (cos(theta[1]) *(0.171* sin(theta[3]) - 0.041* cos(theta[3]) + 0.045)
	 + sin(theta[1])* (0.041* sin(theta[3]) + 0.171*cos(theta[3]) + 0.55))
	 + (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*
	 (sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 - cos(theta[1]) *(0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	 + (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*
	 (cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))
	 - sin(theta[1]) *(0.041* sin(theta[3]) + 0.171* cos(theta[3]))))
	 - 1.505* (2 *(cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))*
	 (sin(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 - cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	 + 2* (cos(theta[1])* (0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	 + sin(theta[1]) *(0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	* (cos(theta[1]) *(0.041* cos(theta[3]) - 0.171* sin(theta[3]))
	   - sin(theta[1]) *(0.041* sin(theta[3]) + 0.171* cos(theta[3]))));

	 double c122=c212;
	 double c222=1.505* (2* (sin(theta[1])* (0.041* cos(theta[3])
	 - 0.171* sin(theta[3])) - cos(theta[1])* (- 0.041* sin(theta[3])
	 - 0.171* cos(theta[3])))* (sin(theta[1]) *(0.041* sin(theta[3])
	 + 0.171* cos(theta[3])) - cos(theta[1])* (0.041* cos(theta[3])
	 - 0.171* sin(theta[3]))) + 2* (sin(theta[1])*
	 (- 0.041* sin(theta[3]) - 0.171* cos(theta[3]))
	 + cos(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3])))
	* (cos(theta[1]) *(0.041* sin(theta[3]) + 0.171* cos(theta[3]))
	   + sin(theta[1])* (0.041* cos(theta[3]) - 0.171* sin(theta[3]))));



	 Eigen::Matrix<double, 4, 2> Ctmp;
	 Ctmp<<c111, c211, c112, c212, c121, c221, c122, c222;
	 Eigen::Vector2d thetaDot2jnt;
	 thetaDot2jnt<<thetaDot[1], thetaDot[3];




	 Eigen::Vector4d tmp=Ctmp*thetaDot2jnt;

	 int cnt=0;
	 for (int i=0; i<=1; i++){
		 for (int j=0; j<=1; j++){
		 Cmatrix(i,j)=tmp[cnt];
		 cnt=cnt+1;
		 }

	 }

	}

// Cmatrix Definition ends///////////////////////////////////////////////////////


// Phi Vector definition ///////////////////////////////////////////////////////////
template<typename InputType1, typename InputType2>
void Dynamics<InputType1, InputType2>::makePhi(const InputType1& theta){
	double phi11=29.5034 *(- cos(theta[1])* (0.171* sin(theta[3])
	 - 0.041* cos(theta[3]) + 0.045) - sin(theta[1])*
	 (0.041* sin(theta[3]) + 0.171* cos(theta[3]) + 0.55))
	 + 55.637* (- 0.119* sin(theta[1])
		    - 0.0005* cos(theta[1]));
	 double phi21=29.5034* (cos(theta[1])* (0.041* cos(theta[3])
	 - 0.171* sin(theta[3])) - sin(theta[1])* (0.041* sin(theta[3])
				       + 0.171* cos(theta[3])));

phi<<phi11, phi21;
}

template<typename InputType1, typename InputType2>
void Dynamics<InputType1, InputType2>::computeModel(const InputType1& theta, const InputType2& thetaDot/*, Eigen::Vector2d& torque, Eigen::Vector4d& X*/){

			Eigen::Matrix2d tmp_gx;
			//Eigen::Vector4d sta, u;
			//sta<<theta(1), theta(3), thetaDot(1), thetaDot(3);
			//u<<0,0,torque(0), torque(1);
			makemassMat(theta);
			makeCmat(theta, thetaDot);
//			makePhi(theta); // to make it unaffected from gravity
			//phi.Zero();

			//massMatrix4.block<2,2>(0,0)=tmp;
			massMatrix4.block<2,2>(0,0).setIdentity();

			massMatrix4.block<2,2>(2,2)=massMatrix;
			Cmatrix4.block<2,2>(0,2).setIdentity();   //when this comes right side, it's  POSITIVE
			Cmatrix4.block<2,2>(2,2) = -Cmatrix;
			phi4.block<2,1>(2,0) = -phi;
//			massMatrix.inverse();
			gx.block<2,2>(2,0)=massMatrix.inverse();


			//svd representation of previous line
//			massMatrix.svd().solve(Eigen::Matrix2d::Identity(), &tmp_gx);
//			gx.block<2,2>(2,0)=tmp_gx;
			//fx = massMatrix4.inverse()*Cmatrix4;
//			massMatrix4.svd().solve(Cmatrix4, &fx);
			//X= massMatrix4.inverse()*(u +Cmatrix4*sta + phi4);
			//massMatrix4.svd().solve(u + Cmatrix4*sta + phi4, &Xdot);

			return;
}



} /* namespace Sam */
