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
#include <tsk.h>
namespace barrett {
namespace systems {
template<int NSTATES, int NR, typename InputType1, typename InputType2,
		 typename OutputType = typename InputType1::actuator_type,
		 typename MathTraits = math::Traits<InputType1> >
class myController : public Controller<InputType1, OutputType> {
public:
	typedef typename MathTraits::unitless_type unitless_type_jp;

	typedef  typename units::JointTorques<NSTATES/2>::type criticTorque_type;
	typedef typename math::Vector<2*NSTATES, double>::type gx_type;
	typedef typename math::Vector<NR, double>::type vec_type;
public:

	   System::Input<InputType2> wamJVIn;


	   System::Output<criticTorque_type> criticconOp;
//	   System::Output<gx_type> gxrow;
//	   System::Output<unitless_type_jp> pEout, vEout, refout;
//	   System::Output<vec_type> Rout1, Rout2;
	   //To make available critic torque value
protected:
	   typename System::Output<criticTorque_type>::Value* criticconOpValue;
//	   typename System::Output<gx_type>::Value* gxrowVal;

//	   typename System::Output<unitless_type_jp>::Value* pEoutVal, *vEoutVal, *refoutVal;
//	   typename System::Output<vec_type>::Value* RoutVal1, *RoutVal2;
private:
	   criticTorque_type ctorque;


private:
//    QuadTS critic;
	size_t rCount;
    unsigned int NU;
    Eigen::Matrix<double, NSTATES/2, NSTATES/2> R, Rinv;

    Eigen::Matrix<double, NSTATES, 1> x, err_states_wam;
    Eigen::Matrix<double, NSTATES/2, 1> u;
    Sam::tsk<NSTATES, NSTATES/2, NR> fuz;

    Eigen::Matrix<double,NSTATES,NSTATES/2> gx_fuz;

//    vec_type rout1, rout2;

//    gx_type  gx_matrow;
public:

	explicit myController(const std::string& sysName = "myController");
	explicit myController(const libconfig::Setting& setting, const std::string& sysName = "myController" );

	virtual ~myController(){
		this->mandatoryCleanUp();
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


	Eigen::Matrix<double, NSTATES/2, NSTATES/2>& getRinv(){
		return Rinv;
	}



  //////////////////////////////////////////////////////////////////////////////////////
  
	void setControlSignalLimit(const unitless_type_jp& csSaturations);
	OutputType& getControlSignalLimit() {  return controlSignalLimit;  }
  const OutputType& getControlSignalLimit() const {  return controlSignalLimit;  }
  
  void resetRcount(){rCount=0;}



protected:
	void setSamplePeriod(double timeStep);

	virtual void operate();
	virtual void onExecutionManagerChanged() {
		Controller<InputType1, OutputType>::onExecutionManagerChanged();  // First, call super
		getSamplePeriodFromEM();
	}
	void setR(const char* RmatFile);
	double T_s;
	InputType1 error_jp, error_1_jp, currentJp, prevJp;
	InputType2 error_jv, error_1_jv, currentJv;
  ///////////////////////////////////////////////////////////////////////////////
  //Test for myController
  //=========================================================================
	unitless_type_jp intError, intErrorLimit;
	unitless_type_jp kp, ki, kd, pE, vE, ref;
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
