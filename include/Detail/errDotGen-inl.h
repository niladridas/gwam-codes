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
errDotGen<InputType1, MathTraits>::errDotGen():System("sysName"), feedbackInput(this),timeInput(this), errorDotOut(this, &errorDotOutVal), errorOut(this, &errorOutVal), refAcc(this, &refAccVal),
refOut(this, &refOutVal), timeOutput(this, &timeOutputVal), currentJp(0.0), jp_s(0.0), jp_f(0.0), refJp(0.0), refJp_b(0.0), velocity_s(0.0), velocity_f(0.0),acc_s(0.0),
acc_f(0.0), err(0.0),   prevErr(0.0), eDot(0.0), refacc(0.0), Ts(0.0), max_vel(0.0), m_fact(0.0), t(0.0), totalTime(0.0), designMat(Eigen::Matrix<double, 6, 6>()),
timePolynom(Eigen::Matrix<double, 1, 6>()), timePolNome2(Eigen::Matrix<double, 1, 6>()),tmp1(Eigen::Matrix<double, 1,1>()), tmp2(0.0)
{
	a = new Eigen::Matrix<double, 6, 1>[DOF];
	dataVec = new Eigen::Matrix<double, 6, 1>[DOF];
	getSamplePeriodFromEM();


}



template<typename InputType1, typename MathTraits>
errDotGen<InputType1, MathTraits>::errDotGen(const InputType1& StartPos, const InputType1& Des, const double mFact, const double& maxVel, const double vel_start, const double vel_final,
		const double acc_start, const double acc_final,	 const std::string& sysName) :
		System("sysName"), feedbackInput(this),timeInput(this), errorDotOut(this, &errorDotOutVal), errorOut(this, &errorOutVal),refAcc(this, &refAccVal),
		refOut(this, &refOutVal), timeOutput(this, &timeOutputVal), currentJp(0.0), jp_s(0.0), jp_f(0.0), refJp(0.0), refJp_b(0.0), velocity_s(0.0), velocity_f(0.0),acc_s(0.0),
		acc_f(0.0), err(0.0),   prevErr(0.0), eDot(0.0),refacc(0.0), Ts(0.0), max_vel(0.0), m_fact(0.0), t(0.0), totalTime(0.0), designMat(Eigen::Matrix<double, 6, 6>()),
		timePolynom(Eigen::Matrix<double, 1, 6>()),timePolNome2(Eigen::Matrix<double, 1, 6>()), tmp1(Eigen::Matrix<double, 1,1>()), tmp2(0.0)
{

	a = new Eigen::Matrix<double, 6, 1>[DOF];
	dataVec = new Eigen::Matrix<double, 6, 1>[DOF];
	getSamplePeriodFromEM();

	setAll(StartPos, Des, mFact, maxVel);
//	calcReachingTime();
//
//
//		designMat<< 1, 0, 0, 0, 0, 0,
//					0, 1, 0, 0, 0, 0,
//					0, 0, 2, 0, 0, 0,
//					1, totalTime, pow(totalTime, 2), pow(totalTime,3), pow(totalTime,4), pow(totalTime, 5),
//					0, 1, 2*totalTime, 3*pow(totalTime,2), 4*pow(totalTime,3), 5*pow(totalTime,4),
//					0, 0, 2, 6*totalTime, 12*pow(totalTime,2), 20*pow(totalTime,3);
//
//		for(size_t i=0; i<DOF; i++){
//			dataVec[i] << jp_s[i], velocity_s[i], acc_s[i], jp_f[i], velocity_f[i], acc_f[i];
//
//			a[i] = designMat.inverse()*dataVec[i];
//		}


}

//template<typename InputType1, typename MathTraits>
//void errDotGen<InputType1, MathTraits>::setCurrentPos(const InputType1& J){}
//
//template<typename InputType1, typename MathTraits>
//void errDotGen<InputType1, MathTraits>::setDesPos(const InputType1& D){}


template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::calcReachingTime(){
	size_t jointId =0;
	double jVal = 0;
	for(size_t i=0; i<DOF; i++){

		if(std::abs(jp_f[i]-jp_s[i])>jVal){
			jVal = std::abs(jp_f[i]-jp_s[i]);
			jointId = i;
		}
	}
	totalTime = m_fact*std::abs((jp_f[jointId] - jp_s[jointId])/max_vel);
	std::cout<<"Total Time Calculated: "<<totalTime<<std::endl;
}



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
void errDotGen<InputType1, MathTraits>::setAll(const InputType1& startingPos, const InputType1& newDesp, const double mFact, const double& maxVel){
	jp_s= startingPos;
	jp_f = newDesp;

	if(maxVel<=MAX_VEL_LIMIT)
		max_vel = maxVel;
	else
		max_vel = MAX_VEL_LIMIT;

	if(mFact>=MIN_MFACT_LIMIT)
		m_fact=mFact;
	else
		m_fact =  MIN_MFACT_LIMIT;

	for(size_t i=0; i<DOF; i++){

		refJp[i]=0;
		err[i]=0;
		prevErr[i]=0;
		eDot[i]=0;
	}
//std::cout<<"=== prev==="<<prevErr<<std::endl;
	totalTime =0.0;
	t= 0.0;
	calcReachingTime();

	if(std::abs(startingPos[1]-newDesp[1])< 0.3 && std::abs((startingPos[1]-newDesp[1])/totalTime)>0.13)
	{totalTime = 2*totalTime;}

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

//	print();
}



template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::print() const
{
		std::cout<<"Current jP: "<<currentJp<<std::endl<<" jP_s: "<<jp_s<<std::endl<<" jP_f: "<<jp_f<<std::endl<<" refJp: "<<refJp<<std::endl<<" err: "<<err<<std::endl
				<<" prevErr: "<<prevErr<<std::endl<<" eDot: "<<eDot<<std::endl<<"Ts, max vel, totalTime, mfact"<<Ts<<", "<<max_vel<<", "<<totalTime<<", "<<m_fact<<std::endl;

		std::cout<<"==========================Coeff============================="<<std::endl;
		for (size_t i=0; i<DOF; i++){
			std::cout<<"a"<<i<<" : "<<a[i].transpose()<<std::endl;
		}
		std::cout<<"==========================Desmat============================="<<std::endl;
		std::cout<<designMat<<std::endl;
		std::cout<<"==========================Desmat End============================="<<std::endl;
	}

template<typename InputType1, typename MathTraits>
void errDotGen<InputType1, MathTraits>::operate(){
	typedef MathTraits MT;
//	t=Ts*(timeInput.getValue());
	t=timeInput.getValue();
	if (t>totalTime){
		t=totalTime;
	}

	currentJp = feedbackInput.getValue();
	timePolynom << 1, t, pow(t,2), pow(t,3), pow(t,4), pow(t,5);
	timePolNome2 << 0, 0, 2, 6*t, 12*pow(t,2), 20*pow(t,3);
	for(size_t i=0; i<DOF; i++){
		 tmp1 = timePolynom*a[i];
		 tmp2 = tmp1(0,0);
		 refJp[i]=tmp2;
		 tmp1 = timePolNome2*a[i];
		 tmp2 = tmp1(0,0);
		 refacc[i] = tmp2;
	}




	err = MT::sub(refJp, currentJp);
	eDot=MT::div(MT::sub(err, prevErr), Ts);
errorOutVal->setData(&err);
errorDotOutVal->setData(&eDot);
refAccVal->setData(&refacc);
refOutVal->setData(&refJp);
timeOutputVal->setData(&t);
prevErr = err;


}

}
}
#endif /* ERRDOTGEN_INL_H_ */
