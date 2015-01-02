#include <cassert>

#include <libconfig.h++>

#include <barrett/math/utils.h>
#include <barrett/thread/abstract/mutex.h>


namespace barrett {
namespace systems {
template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::myController(const std::string& sysName) :
	Controller<InputType1, OutputType>(sysName),
	wamJVIn(this), /*gx(this),*/criticconOp(this, &criticconOpValue),//gxrow(this, &gxrowVal),//wamJVref(this),
	critic("param.txt", "Wmatrices.txt", "center.txt"),NU(NSTATES/2),gx_wam(Eigen::Matrix<double,NSTATES,NSTATES/2>()),gx_fuz(Eigen::Matrix<double,NSTATES,NSTATES/2>()),
	/*gx_matrow(0.0),*/T_s(0.0), error_jp(0.0), error_1_jp(0.0),currentJp(0.0), error_jv(0.0), error_1_jv(0.0),currentJv(0.0),intError(0.0),intErrorLimit(0.0),
	kp(0.0), ki(0.0), kd(0.0),controlSignal(0.0), controlSignalLimit(0.0)
	{
		getSamplePeriodFromEM();
		Rinv = allocate<double>((NU)*(NU));
		eye(Rinv,NU);
		gTld = allocate<double>(NU);
		gth = allocate<double>(NSTATES*(NU));
		delx = allocate<double>(NSTATES);
		err_states_wam = allocate<double>(NSTATES);
		states_wam = allocate<double>(NSTATES);
		controlSignal_critic = allocate<double>(NU);


	}



template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::myController(
	const libconfig::Setting& setting, const std::string& sysName) :
	Controller<InputType1, OutputType>(sysName),
	wamJVIn(this), /*gx(this),*/criticconOp(this, &criticconOpValue),//gxrow(this, &gxrowVal),//wamJVref(this),
	critic("param.txt", "Wmatrices.txt", "center.txt"),NU(NSTATES/2), gx_wam(Eigen::Matrix<double,NSTATES,NSTATES/2>()),
	gx_fuz(Eigen::Matrix<double,NSTATES,NSTATES/2>()),/*gx_matrow(0.0),*/ T_s(0.0), error_jp(0.0), error_1_jp(0.0), currentJp(0.0), error_jv(0.0), error_1_jv(0.0),currentJv(0.0), intError(0.0),intErrorLimit(0.0),
	kp(0.0), ki(0.0), kd(0.0),controlSignal(0.0), controlSignalLimit(0.0)
{
	getSamplePeriodFromEM();
	setFromConfig(setting);
	Rinv = allocate<double>((NU)*(NU));
	eye(Rinv,NU);
	gTld = allocate<double>(NU);
	gth = allocate<double>(NSTATES*(NU));
	delx = allocate<double>(NSTATES);
	err_states_wam = allocate<double>(NSTATES);
	states_wam = allocate<double>(NSTATES);
	controlSignal_critic = allocate<double>(NU);

}

template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setFromConfig(const libconfig::Setting& setting)
{


  ///////////////////////////////////////////////////////////////
  // Test for myController
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

template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setSamplePeriod(double timeStep)
{
	T_s = timeStep;
}


  
  ///////////////////////////////////////////////////////////////////////////////
  //Test for myController
  //============================================================================
  template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setKp(const unitless_type_jp& proportionalGains)
{
	std::cout << "Our previous proportional control gain: kp = " << kp << std::endl;
	kp = proportionalGains;
	std::cout << "Our new proportional control gain: kp = " << kp << std::endl;
}

template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setKi(const unitless_type_jp& integralGains)
{
	ki = integralGains;
}

template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setKd(const unitless_type_jp& derivitiveGains)
{
	kd = derivitiveGains;
}

template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setIntegratorState(
		const unitless_type_jp& integratorState)
{
	// intError is written and read in operate(), so it needs to be locked.
	BARRETT_SCOPED_LOCK(this->getEmMutex());
	intError = integratorState;
}

template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setIntegratorLimit(
		const unitless_type_jp& intSaturations)
{
	intErrorLimit = intSaturations;
}


 template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
inline void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::resetIntegrator()
{
	setIntegratorState(unitless_type_jp(0.0));
}
  /////////////////////////////////////////////////////////////////////////////////////////



  

 template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::setControlSignalLimit(
		const unitless_type_jp& csSaturations)
{
	controlSignalLimit = csSaturations;
}


  
template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::operate()
{
	typedef MathTraits MT;
		currentJp = this->feedbackInput.getValue();
		error_jp = MT::sub(this->referenceInput.getValue(), currentJp);
		currentJv = wamJVIn.getValue();
		//collect gx
//		gx_fuz = gx.getValue();
//		///O/P gx row
//		gx_matrow(0,0)=gx_fuz(1,0);
//		gx_matrow(1,0)=gx_fuz(1,1);
//		gx_matrow(2,0)=gx_fuz(3,0);
//		gx_matrow(3,0)=gx_fuz(3,1);
//		gxrowVal->setData(&gx_matrow);

		// Original model
		dyn.computeModel(this->feedbackInput.getValue(), currentJv);
		gx_wam=dyn.getGx();
		gth[0] = 0;
		gth[1] = 0;

		gth[2] = gx_wam(2,0);            //2nd row of gx
		gth[3] = gx_wam(2,1);


		gth[4] = 0;
		gth[5] = 0;


		gth[6] = gx_wam(3,0);
		gth[7] = gx_wam(3,1);



	//////////////////////////////////////////////////////////////////
  // Critic Controller



	states_wam[0] = currentJp[1];
	states_wam[1] = currentJv[1];
	states_wam[2] = currentJp[3];
	states_wam[3] = currentJv[3];

	err_states_wam[0]= -error_jp[1];
	err_states_wam[1]= states_wam[1]; // As ref velocity is 0
	err_states_wam[2]= -error_jp[3];
	err_states_wam[3]= states_wam[3];

	critic.computeDelX(err_states_wam, delx, states_wam);

	  multiply(gth, delx, gTld, NU,1,NSTATES, TRANS);
	  multiply(Rinv, gTld, controlSignal_critic, NU,1,NU);

	    multiply(controlSignal_critic, double(-0.5), NU);

	    ctorque[0] = controlSignal_critic[0];
	    ctorque[1] =controlSignal_critic[1];

  //=================================================================

	intError = MT::add(intError, MT::mult(ki, MT::mult(T_s, error_1_jp)));
	if (intErrorLimit != MT::zero()) {
		intError = math::saturate(intError, intErrorLimit);
	}

	controlSignal = MT::add(MT::mult(kp, error_jp),
							MT::add(intError,
								MT::mult(kd, MT::div(MT::sub(error_jp, error_1_jp), T_s))));
//	if (controlSignalLimit != MT::zero()) {
//		controlSignal = math::saturate(controlSignal, controlSignalLimit);
//	}

	error_1_jp = error_jp;

//	controlSignal[1]=ctorque[0];
//	controlSignal[3]=ctorque[1];
	if (controlSignalLimit != MT::zero()) {
			controlSignal = math::saturate(controlSignal, controlSignalLimit);
		}
	this->controlOutputValue->setData(&controlSignal);
	criticconOpValue->setData(&ctorque);


  /////////////////////////////////////////////////////////////////////////

  
}


  
  template<int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
  void myController<NSTATES, InputType1, InputType2, OutputType, MathTraits>::getSamplePeriodFromEM()
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
