/*
 * compTorqueController.h
 *
 *  Created on: 09-Sep-2013
 *      Author: mobman
 */

#ifndef COMPTORQUECONTROLLER_H_
#define COMPTORQUECONTROLLER_H_



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
class compTorqueController : public Controller<InputType1, OutputType> {
public:
	typedef typename MathTraits::unitless_type unitless_type_jp;

	typedef  typename units::JointTorques<NSTATES/2>::type compTorque_type;

	typedef typename Eigen::Matrix<double, NSTATES/2, NSTATES/2> cmptMat_type;
	typedef typename math::Traits<compTorque_type>::unitless_type unitless_type_cmpTorque;

public:
	   System::Input<InputType2> wamJVIn;
	   System::Output<compTorque_type> compConOp;
	   System::Output<unitless_type_jp> pEout, vEout, refOut;

protected:
	   typename System::Output<compTorque_type>::Value* compConOpValue;


	   typename System::Output<unitless_type_jp>::Value* pEoutVal, *vEoutVal, *refOutVal;
//	   typename System::Output<vec_type>::Value* RoutVal1, *RoutVal2;
private:
	   compTorque_type ctorque;


private:

	size_t rCount;
	Eigen::Matrix<double, NSTATES, 1> x, err_states_wam;
    Eigen::Matrix<double, NSTATES/2, 1> u, qDesDD, eVel, ePos, qDot, kp_comp, kd_comp;

    Sam::tsk<NSTATES, NSTATES/2, NR> fuz;

    Eigen::Matrix<double,NSTATES,NSTATES> model_A;
    Eigen::Matrix<double,NSTATES,NSTATES/2> model_B;


    cmptMat_type a, b, massMat, Cmat;
    unitless_type_jp refJp;

public:

	explicit compTorqueController(const std::string& sysName = "compTorqueController");
	explicit compTorqueController(const libconfig::Setting& setting, const std::string& sysName = "compTorqueController" );

	virtual ~compTorqueController(){
		this->mandatoryCleanUp();
	}

  void setFromConfig(const libconfig::Setting& setting);

  ////////////////////////////////////////////////////////////////////////////////////////
  //  test for compTorqueController
  //================================================================================
  void setKp(const unitless_type_jp& proportionalGains);
  void setKpComp(const Eigen::Matrix<double, NSTATES/2,1>& kpGains);
  void setKi(const unitless_type_jp& integralGains);
  void setKd(const unitless_type_jp& derivitiveGains);
  void setKdComp(const Eigen::Matrix<double, NSTATES/2,1>& kdGains);
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

  void resetRcount(){rCount=0;}

private:
  void setCompTorqeGains(const char* RmatFile);


protected:
	void setSamplePeriod(double timeStep);

	virtual void operate();
	virtual void onExecutionManagerChanged() {
		Controller<InputType1, OutputType>::onExecutionManagerChanged();  // First, call super
		getSamplePeriodFromEM();
	}

	double T_s;
	InputType1 error_jp, error_1_jp, currentJp, prevJp;
	InputType2 error_jv, error_1_jv, currentJv;
  ///////////////////////////////////////////////////////////////////////////////
  //Test for compTorqueController
  //=========================================================================
	unitless_type_jp intError, intErrorLimit;
	unitless_type_jp kp, ki, kd, pE, vE;

  ///////////////////////////////////////////////////////////////////////////////

	OutputType controlSignal, controlSignalLimit;
//	unitless_type_jp controlSignal, controlSignalLimit;

	void getSamplePeriodFromEM();

private:
	DISALLOW_COPY_AND_ASSIGN(compTorqueController);

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(MathTraits::RequiresAlignment);






};

}
}
#include <Detail/compTorqueController-inl.h>




#endif /* COMPTORQUECONTROLLER_H_ */
