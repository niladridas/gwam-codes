/*
 * errDotGen-inl.h
 *
 *  Created on: 15-Sep-2013
 *      Author: mobman
 */

#ifndef ERRDOTGEN_INL_H_
#define ERRDOTGEN_INL_H_

#include <cassert>
#include <barrett/math/utils.h>
#include <barrett/thread/abstract/mutex.h>
#include <cmath>
#include <vector>
namespace barrett {
namespace systems {

template<typename InputType1, typename MathTraits>
errDotGen<InputType1, MathTraits>::errDotGen(const InputType1& StartPos, const InputType1& Des, const double maxVel, const double& mFact, const double vel_start, const double vel_final,
		const double acc_start, const double acc_final,	 const std::string& sysName) :
		System(sysName), feedbackInput(this),timeInput(this), errorDotOut(this, &errorDotOutVal), errorOut(this, &errorOutVal), refOut(this, &refOutVal), currentJp(0.0), jp_s(StartPos),
		jp_f(Des), refJp(0.0), refJp_b(0.0), err(0.0),   prevErr(0.0), eDot(0.0), velocity_s(vel_start), velocity_f(vel_final),acc_s(acc_start), acc_f(acc_final),
		Ts(0.0), max_vel(maxVel), t(0.0), totalTime(0.0), designMat(Eigen::Matrix<double, 6, 6>()), timePolynom(Eigen::Matrix<double, 1, 6>()), tmp1(Eigen::Matrix<double, 1,1>()), tmp2(0.0)
{

	a = new Eigen::Matrix<double, 6, 1>[DOF];
	dataVec = new Eigen::Matrix<double, 6, 1>[DOF];
	getSamplePeriodFromEM();
	calcReachingTime(mFact);

//	calcTbend();
//	calAcc();
//	setCoeffs();

		designMat<< 1, 0, 0, 0, 0, 0,
					0, 1, 0, 0, 0, 0,
					0, 0, 2, 0, 0, 0,
					1, totalTime, pow(totalTime, 2), pow(totalTime,3), pow(totalTime,4), pow(totalTime, 5),
					0, 1, 2*totalTime, 3*pow(totalTime,2), 4*pow(totalTime,3), 5*pow(totalTime,4),
					0, 0, 2, 6*totalTime, 12*pow(totalTime,2), 20*pow(totalTime,3);

		for(size_t i=0; i<DOF; i++){
			dataVec[i] << jp_s[i], velocity_s[i], acc_s[i], jp_f[i], velocity_f[i], acc_f[i];

			a[i] = designMat.inverse()*dataVec[i];
		}


}

//template<typename InputType1, typename MathTraits>
//void errDotGen<InputType1, MathTraits>::setCurrentPos(const InputType1& J){}
//
//template<typename InputType1, typename MathTraits>
//void errDotGen<InputType1, MathTraits>::setDesPos(const InputType1& D){}


template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::calcReachingTime(const double& mFact){
	size_t jointId =0;
	double jVal = 0;
	for(size_t i=0; i<DOF; i++){

		if(std::abs(jp_f[i]-jp_s[i])>jVal){
			jVal = std::abs(jp_f[i]-jp_s[i]);
			jointId = i;
		}
	}
	totalTime = mFact*std::abs(mFact*(jp_f[jointId] - jp_s[jointId])/max_vel);
}

//template<typename InputType1, typename MathTraits>
//void errDotGen<InputType1, MathTraits>::calAcc(){
//	for(size_t i=0; i<DOF; i++ ){
//	acc_actual[i] = max_vel/tb[i];
//	}
//}
//
//template<typename InputType1, typename MathTraits>
//void errDotGen<InputType1, MathTraits>::calcTbend(){
//	for(size_t i=0; i<DOF; i++ ){
//		tb[i] = std::abs((jp_s[i]- jp_f[i] + max_vel*totalTime)/max_vel);
//	}
//}

template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::setSamplePeriod(double timeStep)
{
	Ts = timeStep;
}

template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::getCurrentPos(){
	currentJp = feedbackInput.getValue();
}

template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::getSamplePeriodFromEM()
{
	if (this->hasExecutionManager()) {
		assert(this->getExecutionManager()->getPeriod() > 0.0);
		setSamplePeriod(this->getExecutionManager()->getPeriod());
	} else {
		setSamplePeriod(0.0);
	}
}


template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::setCoeffs(){
	for(size_t i=0; i<DOF;i++){

	}
}


template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::operate(){
	typedef MathTraits MT;
	t=Ts*(timeInput.getValue());
	currentJp = feedbackInput.getValue();
	timePolynom << 1, t, pow(t,2), pow(t,3), pow(t,4), pow(t,5);
	for(size_t i=0; i<DOF; i++){
		 tmp1 = timePolynom*a[i];
		 tmp2 = tmp1(0,0);
		 refJp[i]=tmp2;
	}

	if (t>totalTime){
		refJp = unitless_type(jp_f);
	}

	err = MT::sub(refJp, currentJp);
	eDot=MT::div(MT::sub(err, prevErr), Ts);
errorOutVal->setData(&err);
errorDotOutVal->setData(&eDot);
refOutVal->setData(&refJp);
prevErr = err;


}

}
}
#endif /* ERRDOTGEN_INL_H_ */
