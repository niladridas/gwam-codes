/*
 * errDotGen.h
 *
 *  Created on: 14-Sep-2013
 *      Author: mobman
 */

#ifndef ERRDOTGEN_H_
#define ERRDOTGEN_H_

#include <eigen3/Eigen/Core>

#include <barrett/detail/ca_macro.h>
#include <barrett/math/traits.h>
#include <barrett/systems/abstract/execution_manager.h>
//#include <barrett/systems/callback.h>

namespace barrett {
namespace systems {
template<typename InputType1,
		 typename MathTraits = math::Traits<InputType1> >
class errDotGen : public System{//Controller<InputType1> {

public:
	typedef typename MathTraits::unitless_type unitless_type;


public:
	System::Input<InputType1> feedbackInput;
	System::Input<double> timeInput;
	System::Output<unitless_type> errorDotOut;
	System::Output<unitless_type> errorOut;
	System::Output<unitless_type> refAcc;
	System::Output<InputType1> refOut;
	System::Output<double> timeOutput;
protected:
	typename System::Output<unitless_type>::Value* errorDotOutVal;
	typename System::Output<unitless_type>::Value* errorOutVal;
	typename System::Output<unitless_type>::Value* refAccVal;
	typename System::Output<InputType1>::Value* refOutVal;
	typename System::Output<double>::Value* timeOutputVal;
private:
	InputType1 currentJp, jp_s,  jp_f, refJp, refJp_b;
	unitless_type velocity_s, velocity_f, acc_s, acc_f,err, prevErr, eDot, refacc;
//	unitless_type a0, a1, a2, a3, a4, a5, tb;
	double Ts, max_vel, m_fact;
	double t, totalTime;
	static const size_t DOF = 7;

	Eigen::Matrix<double, 6, 6> designMat;
	Eigen::Matrix<double, 6, 1>* a, *dataVec;
	Eigen::Matrix<double, 1, 6> timePolynom, timePolNome2;
	Eigen::Matrix<double, 1,1> tmp1;
	double tmp2;

	static const double MAX_VEL_LIMIT = 2.0;
	static const double MIN_MFACT_LIMIT = 1.0;
public:
	explicit errDotGen();
	explicit errDotGen(const InputType1& StartPos, const InputType1& Des, const double mFact = 1.5, const double& maxVel=2,  const double vel_start=0,
			const double vel_final=0, const double acc_start=0, const double acc_final= 0, const std::string& sysName = "errDotGen");
	virtual ~errDotGen(){
		this->mandatoryCleanUp();

		delete [] a;
		delete [] dataVec;
	}

	void setAll(const InputType1& startingPos, const InputType1& newDesp, const double mFact, const double& maxVel);

	void print() const;
	void setSamplePeriod(double timeStep);

protected:

//	void setCurrentPos(const InputType1& J);
//	void setDesPos(const InputType1& D);
	void calcReachingTime();
//	void calAcc();
//	void calcTbend();
	void getCurrentPos();
	void setCoeffs();
private:

	virtual void operate();
	virtual void onExecutionManagerChanged() {
		System::onExecutionManagerChanged();  // First, call super
		getSamplePeriodFromEM();
	}

	void getSamplePeriodFromEM();
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(MathTraits::RequiresAlignment)

};
}
}

#include <Detail/errDotGen-inl.h>
#endif /* ERRDOTGEN_H_ */
