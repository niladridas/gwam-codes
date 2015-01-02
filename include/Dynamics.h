/*
 * Dynamics.h
 *
 *  Created on: 25-Aug-2012
 *      Author: mobman
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_
#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/Array>
//#include <eigen2/Eigen/Dense>
#include <eigen2/Eigen/LU>
#include <math.h>
namespace Sam {

template<typename InputType1, typename InputType2>
class Dynamics {
public:
	Dynamics();
	virtual ~Dynamics(){}


void makemassMat( const InputType1& theta);
void makeCmat(const InputType1& theta, const InputType2& thetaDot);
void makePhi(const InputType1& theta);

// Twolink model

void computeModel(const InputType1& theta, const InputType2& thetaDot/*, Eigen::Vector2d& torque, Eigen::Vector4f& X*/);
const Eigen::Matrix2d& getmassMat() const
{return massMatrix;}
const Eigen::Matrix2d& getCmat() const {return Cmatrix;}
const Eigen::Vector2d& getPhi() const {return phi;}
const Eigen::Matrix<double, 4, 2>& getGx() const {return gx;}
const Eigen::Matrix4d& getFx() const {return fx;}



protected:
Eigen::Matrix2d massMatrix;
Eigen::Matrix2d Cmatrix;
Eigen::Vector2d phi;

Eigen::Matrix4d massMatrix4, fx;   //for state space model;
Eigen::Matrix4d Cmatrix4;
Eigen::Vector4d phi4;
Eigen::Matrix<double, 4, 2>  gx;
};

} /* namespace Sam */
#include<Detail/Dynamics-inl.h>
#endif /* DYNAMICS_H_ */
