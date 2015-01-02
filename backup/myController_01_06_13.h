#ifndef MYCONTROLLER_H_
#define MYCONTROLLER_H_


#include <eigen3/Eigen/Core>
#include <libconfig.h++>

#include <barrett/detail/ca_macro.h>
#include <barrett/math/traits.h>
#include <barrett/systems/abstract/execution_manager.h>
#include <barrett/systems/abstract/controller.h>
#include <QuadTS.h>
# include <Dynamics.h>
#include <system.h>

namespace barrett {
namespace systems {
template<int NSTATES, typename InputType1, typename InputType2,
		 typename OutputType = typename InputType1::actuator_type,
		 typename MathTraits = math::Traits<InputType1> >
class myController : public Controller<InputType1, OutputType> {
public:
	typedef typename MathTraits::unitless_type unitless_type_jp;

	typedef  typename units::JointTorques<NSTATES/2>::type criticTorque_type;
	typedef typename math::Vector<NSTATES, double>::type gx_type;
public:

	   System::Input<InputType2> wamJVIn;

//	   System::Input<Eigen::Matrix<double, NSTATES, NSTATES/2> > gx;
	   System::Output<criticTorque_type> criticconOp;
//	   System::Output<gx_type> gxrow;
	   //To make available critic torque value
protected:
	   typename System::Output<criticTorque_type>::Value* criticconOpValue;
//	   typename System::Output<gx_type>::Value* gxrowVal;
private:
	   criticTorque_type ctorque;


private:
    QuadTS critic;
    unsigned int NU;
    double* Rinv;
    double* gTld;
    double* gth;
    double* delx;
    double* err_states_wam;
    double* states_wam;
    double* controlSignal_critic;

    Eigen::Matrix<double,NSTATES,NSTATES/2> gx_wam, gx_fuz;
    Sam::Dynamics<InputType1, InputType2> dyn;
//    gx_type  gx_matrow;
public:

	explicit myController(const std::string& sysName = "myController");
	explicit myController(const libconfig::Setting& setting, const std::string& sysName = "myController" );

	virtual ~myController(){
		this->mandatoryCleanUp();
		del(Rinv);
		del(gTld);
		del(gth);
		del(delx);
		del(err_states_wam);
		del(states_wam);

		del(controlSignal_critic);

	}

  void setFromConfig(const libconfig::Setting& setting);

  ////////////////////////////////////////////////////////////////////////////////////////
  //  test for myController
  //================================================================================
  void setKp(const unitless_type_jp& proportionalGains);
	void setKi(const unitless_type_jp& integralGains);
	void setKd(const unitless_type_jp& derivitiveGains);
	void setIntegratorState(const unitless_type_jp& integratorState);
	void setIntegratorLimit(const unitless_type_jp& intSaturations);

  void resetIntegrator();

	unitless_type_jp& getKp() {  return kp;  }
	const unitless_type_jp& getKp() const {  return kp;  }
	unitless_type_jp& getKi() {  return ki;  }
	const unitless_type_jp& getKi() const {  return ki;  }
	unitless_type_jp& getKd() {  return kd;  }
	const unitless_type_jp& getKd() const {  return kd;  }
	const unitless_type_jp& getIntegratorState() const {  return intError;  }
	unitless_type_jp& getIntegratorLimit() {  return intErrorLimit;  }
	const unitless_type_jp& getIntegratorLimit() const {  return intErrorLimit;  }



  //////////////////////////////////////////////////////////////////////////////////////
  
	void setControlSignalLimit(const unitless_type_jp& csSaturations);
	OutputType& getControlSignalLimit() {  return controlSignalLimit;  }
  const OutputType& getControlSignalLimit() const {  return controlSignalLimit;  }
  

protected:
	void setSamplePeriod(double timeStep);

	virtual void operate();
	virtual void onExecutionManagerChanged() {
		Controller<InputType1, OutputType>::onExecutionManagerChanged();  // First, call super
		getSamplePeriodFromEM();
	}

	double T_s;
	InputType1 error_jp, error_1_jp, currentJp;
	InputType2 error_jv, error_1_jv, currentJv;
  ///////////////////////////////////////////////////////////////////////////////
  //Test for myController
  //=========================================================================
	unitless_type_jp intError, intErrorLimit;
	unitless_type_jp kp, ki, kd;
  ///////////////////////////////////////////////////////////////////////////////
  
	OutputType controlSignal, controlSignalLimit;
//	unitless_type_jp controlSignal, controlSignalLimit;

	void getSamplePeriodFromEM();

private:
	DISALLOW_COPY_AND_ASSIGN(myController);

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(MathTraits::RequiresAlignment);






};

}
}
#include <Detail/myController-inl.h>
#endif
