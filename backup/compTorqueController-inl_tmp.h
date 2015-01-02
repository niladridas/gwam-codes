/*
 * compTorqueController-inl.h
 *
 *  Created on: 09-Sep-2013
 *      Author: mobman
 */

#ifndef COMPTORQUECONTROLLER_INL_H_
#define COMPTORQUECONTROLLER_INL_H_

#include <cassert>

#include <libconfig.h++>

#include <barrett/math/utils.h>
#include <barrett/thread/abstract/mutex.h>


namespace barrett {
namespace systems {
template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::compTorqueController(const std::string& sysName) :
	Controller<InputType1, OutputType>(sysName),
	wamJVIn(this), compConOp(this, &compConOpValue),pEout(this, &pEoutVal), vEout(this, &vEoutVal), refOut(this, &refOutVal), ctorque(0.0), rCount(0),  x(Eigen::Matrix<double, NSTATES, 1>()),
	err_states_wam(Eigen::Matrix<double, NSTATES, 1>()), u(Eigen::Matrix<double, NSTATES/2, 1>()),qDesDD(Eigen::Matrix<double, NSTATES/2, 1>()),
	eVel(Eigen::Matrix<double, NSTATES/2, 1>()), ePos(Eigen::Matrix<double, NSTATES/2, 1>()), qDot(Eigen::Matrix<double, NSTATES/2, 1>()),
	kp_comp(Eigen::Matrix<double, NSTATES/2, 1>()), kd_comp(Eigen::Matrix<double, NSTATES/2, 1>()),
	fuz("./FuzzParams/cons.txt","./FuzzParams/weight.txt", "./FuzzParams/model_gaussParams.txt"), model_A(Eigen::Matrix<double,NSTATES,NSTATES>()),
	model_B(Eigen::Matrix<double,NSTATES,NSTATES/2>()), a(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), b(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()),
	massMat(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), Cmat(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), refJp(0.0),
	T_s(0.0), error_jp(0.0), error_1_jp(0.0),currentJp(0.0), prevJp(0.0), error_jv(0.0), error_1_jv(0.0),
	currentJv(0.0),intError(0.0), intErrorLimit(0.0),	kp(0.0), ki(0.0), kd(0.0),controlSignal(0.0), controlSignalLimit(0.0)
	{
		getSamplePeriodFromEM();
		qDesDD.fill(0.0);
		u.fill(0.0);
		x.fill(0.0);
		err_states_wam.fill(0.0);



	}



template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::compTorqueController(const libconfig::Setting& setting, const std::string& sysName) :
wamJVIn(this), compConOp(this, &compConOpValue), pEout(this, &pEoutVal), vEout(this, &vEoutVal), refOut(this, &refOutVal), ctorque(0.0), rCount(0),  x(Eigen::Matrix<double, NSTATES, 1>()),
	err_states_wam(Eigen::Matrix<double, NSTATES, 1>()), u(Eigen::Matrix<double, NSTATES/2, 1>()),qDesDD(Eigen::Matrix<double, NSTATES/2, 1>()),
	eVel(Eigen::Matrix<double, NSTATES/2, 1>()), ePos(Eigen::Matrix<double, NSTATES/2, 1>()), qDot(Eigen::Matrix<double, NSTATES/2, 1>()),
	kp_comp(Eigen::Matrix<double, NSTATES/2, 1>()), kd_comp(Eigen::Matrix<double, NSTATES/2, 1>()),
	fuz("./FuzzParams/cons.txt","./FuzzParams/weight.txt", "./FuzzParams/model_gaussParams.txt"), model_A(Eigen::Matrix<double,NSTATES,NSTATES>()),
	model_B(Eigen::Matrix<double,NSTATES,NSTATES/2>()), a(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), b(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()),
	massMat(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), Cmat(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), refJp(0.0),
	T_s(0.0), error_jp(0.0), error_1_jp(0.0),currentJp(0.0), prevJp(0.0), error_jv(0.0), error_1_jv(0.0),
	currentJv(0.0),intError(0.0), intErrorLimit(0.0),	kp(0.0), ki(0.0), kd(0.0),controlSignal(0.0), controlSignalLimit(0.0)
	{
	getSamplePeriodFromEM();
	setFromConfig(setting);
	setCompTorqeGains("./compTorqueGain.txt");
	qDesDD.fill(0.0);
	u.fill(0.0);
	x.fill(0.0);
	err_states_wam.fill(0.0);
}





template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setFromConfig(const libconfig::Setting& setting)
{


  ///////////////////////////////////////////////////////////////
  // Test for compTorqueController
  //===========================================================
  if (setting.exists("kp")) {
		setKp(unitless_type_jp(setting["kp"]));
	}



  if (setting.exists("ki")) {
	  setKi(unitless_type_jp(setting["ki"]));
  }


  if (setting.exists("kd")) {
  	  setKd(unitless_type_jp(setting["kd"]));
    }


  if (setting.exists("integrator_limit")) {
	  setIntegratorLimit(unitless_type_jp(setting["integrator_limit"]));
  }

	//////////////////////////////////////////////////////////////////


	if (setting.exists("control_signal_limit")) {
		setControlSignalLimit(unitless_type_jp(setting["control_signal_limit"]));
	}




}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setSamplePeriod(double timeStep)
{
	T_s = timeStep;
}



  ///////////////////////////////////////////////////////////////////////////////
  //Test for compTorqueController
  //============================================================================
template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setKp(const unitless_type_jp& proportionalGains)
{
//	std::cout << "Our previous proportional control gain: kp = " << kp << std::endl;
	kp = proportionalGains;
//	std::cout << "Our new proportional control gain: kp = " << kp << std::endl;
}





template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setKi(const unitless_type_jp& integralGains)
{
	ki = integralGains;
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setKd(const unitless_type_jp& derivitiveGains)
{
	kd = derivitiveGains;
}



template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setIntegratorState(
		const unitless_type_jp& integratorState)
{
	// intError is written and read in operate(), so it needs to be locked.
	BARRETT_SCOPED_LOCK(this->getEmMutex());
	intError = integratorState;
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setIntegratorLimit(
		const unitless_type_jp& intSaturations)
{
	intErrorLimit = intSaturations;
}


 template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
inline void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::resetIntegrator()
{
	setIntegratorState(unitless_type_jp(0.0));
}


 template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
 void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setCompTorqeGains(const char* RmatFile){

 	std::vector<std::vector<double> > r;
 	size_t tr, tc;

 	Sam::readFile(RmatFile, r, tr, tc);
 	if(!(tr==NSTATES/2 || tc== NSTATES/2)){

 		std::cout<<"File "<<RmatFile<<" is not proper for Gain matrix"<<std::endl;
 		exit(EXIT_FAILURE);
 	}


 	for (size_t j=0; j<NSTATES/2; j++){
 		kp_comp(j)=r[0][j];
 	}
 	for (size_t j=0; j<NSTATES/2; j++){
 		kd_comp(j)=r[1][j];
 	}

 	std::cout << "Our new proportional control gain: kp_comp = " << kp_comp << std::endl;
 	std::cout << "Our new proportional control gain: kd_comp = " << kd_comp << std::endl;
 }





 template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setControlSignalLimit(
		const unitless_type_jp& csSaturations)
{
	controlSignalLimit = csSaturations;
}



template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::operate()
{
	typedef MathTraits MT;
		currentJp = this->feedbackInput.getValue();
		refJp = this->referenceInput.getValue();
		error_jp = MT::sub(refJp, currentJp);
		currentJv = wamJVIn.getValue();

		pE=unitless_type_jp(error_jp);
		vE=MT::div(MT::sub(error_jp, error_1_jp), T_s);

		// Original model

		x<<currentJp[1],currentJv[1],currentJp[3],currentJv[3];
		u<<controlSignal[1],controlSignal[3];

				fuz.computeFuzModel(x,u);




		model_A = fuz.model_A;
		model_B=fuz.model_B;

		a(0,0) = model_A(1,1);
		a(0,1) = model_A(1,3);
		a(1,0) = model_A(3,1);
		a(1,1) = model_A(3,3);

		b(0,0)=model_B(1,0);
		b(0,1)=model_B(1,1);
		b(1,0)=model_B(3,0);
		b(1,1)=model_B(3,1);
		massMat = b.inverse();
		Cmat = a.inverse()*a;

	//////////////////////////////////////////////////////////////////
  // Critic Controller
	err_states_wam(0)=  /*currentJp[1];//*/-error_jp[1];
	err_states_wam(1)=  /*currentJv[1];//*/-vE[1]; // As e=curr-des
	err_states_wam(2)=  /*currentJp[3];//*/-error_jp[3];
	err_states_wam(3)= 	/*currentJv[3];//*/-vE[3];

	qDot<< currentJv[1], currentJv[3];
	ePos << error_jp[1], error_jp[3];
//	ePos << -currentJp[1], -currentJp[3];
//	eVel << vE[1], vE[3];
	eVel = -qDot;

	u = massMat*(qDesDD + kd_comp.cwise()*eVel + kp_comp.cwise()*ePos) + Cmat*qDot;





  //=================================================================

	intError = MT::add(intError, MT::mult(ki, MT::mult(T_s, error_1_jp)));
	if (intErrorLimit != MT::zero()) {
		intError = math::saturate(intError, intErrorLimit);
	}

//	controlSignal = MT::add(MT::mult(kp, error_jp),
//							MT::add(intError,
//								MT::mult(kd, MT::div(MT::sub(error_jp, error_1_jp), T_s))));

	controlSignal = MT::add(MT::mult(kp, error_jp),
							MT::add(intError,
								MT::mult(kd, vE)));


	//===========================================================================================



		    ctorque[0] = u(0);
		    ctorque[1] =u(1);

	//=============================================================================================

	prevJp=currentJp;
	error_1_jp = error_jp;

//	if(rCount>50){
//
//	controlSignal[1]=ctorque[0];
//	controlSignal[3]=ctorque[1];
//	}
	if (controlSignalLimit != MT::zero()) {
			controlSignal = math::saturate(controlSignal, controlSignalLimit);
		}
	this->controlOutputValue->setData(&controlSignal);
	compConOpValue->setData(&ctorque);
	pEoutVal->setData(&pE);
	vEoutVal->setData(&vE);
	refOutVal->setData(&refJp);
	rCount++;

  /////////////////////////////////////////////////////////////////////////


}



  template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
  void compTorqueController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::getSamplePeriodFromEM()
{
	if (this->hasExecutionManager()) {
		assert(this->getExecutionManager()->getPeriod() > 0.0);
		setSamplePeriod(this->getExecutionManager()->getPeriod());
	} else {
		setSamplePeriod(0.0);
	}
}
}
}



#endif /* COMPTORQUECONTROLLER_INL_H_ */
