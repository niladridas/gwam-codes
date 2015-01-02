#ifndef CRITICCONTROLLER_H_
#define CRITICCONTROLLER_H_


#include <eigen3/Eigen/Core>
#include <libconfig.h++>

#include <barrett/detail/ca_macro.h>
#include <barrett/math/traits.h>
#include <barrett/systems/abstract/execution_manager.h>
#include <barrett/systems/abstract/controller.h>
#include<barrett/systems.h>
#include <barrett/systems/abstract/system.h>
#include <iostream>
#include <QuadTS.h>
#include <system.h>
#include <list>  // For samae type Class members. Same time the constructor of the member should be called in the constructor of the class
namespace barrett{
namespace systems{

template <size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType= typename InputType1::actuator_type,
		typename MathTraits = math::Traits<InputType1> >
class CriticController: public Controller< InputType1, OutputType > {
//	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);


public:
	typedef typename MathTraits::unitless_type unitless_type_jp;
//
public:

	   System::Input<InputType2> wamJVIn;
//	   System::Input<InputType2> wamJVref;
// public:
//
////	   Output<jt_type> wamJTOutput;
//
//   // ========================Private variables for Fuzzy u===================================================
 private:
    QuadTS critic;
    double* R;
    double* Rinv;
    double* Q;
    float* gTld;
    float* delx;
    Eigen::Vector4f states_wam;
    Eigen::Matrix<double,NSTATES,NSTATES/2> gx_wam;
    Sam::Dynamics dyn;
//==================================================================================



//     typename Output<jt_type>::Value* jtOutputValue;
   public:
      explicit CriticController( const std::string& sysName =  "CriticController");
      explicit CriticController(const libconfig::Setting& setting, const std::string& sysName = "CriticController");


     virtual ~CriticController()
     {
       this->mandatoryCleanUp();
       del(Q);
       del(R);
       del(Rinv);
       del(gTld);
       del(delx);
     }

     //========================= PID members ================================================================
     	 void setFromConfig(const libconfig::Setting& setting);
     	void setKp(const unitless_type_jp& proportionalGains);
     	void setKi(const unitless_type_jp& integralGains);
     	void setKd(const unitless_type_jp& derivitiveGains);
     	void setIntegratorState(const unitless_type_jp& integratorState);
     	void setIntegratorLimit(const unitless_type_jp& intSaturations);
     	void setControlSignalLimit(const unitless_type_jp& csSaturations);

     	void resetIntegrator();

     	unitless_type_jp& getKp() {  return kp;  }
     	const unitless_type_jp& getKp() const {  return kp;  }
     	unitless_type_jp& getKi() {  return ki;  }
     	const unitless_type_jp& getKi() const {  return ki;  }
     	unitless_type_jp& getKd() {  return kd;  }
     	const unitless_type_jp& getKd() const {  return kd;  }
     	const unitless_type_jp& getIntegratorState() const {  return intError_jp;  }
     	unitless_type_jp& getIntegratorLimit() {  return intErrorLimit_jp;  }
     	const unitless_type_jp& getIntegratorLimit() const {  return intErrorLimit_jp;  }
     	OutputType& getControlSignalLimit() {  return controlSignalLimit;  }
     	const OutputType& getControlSignalLimit() const {  return controlSignalLimit;  }


     //=================================================================================



   protected:
      	void setSamplePeriod(double timeStep);



      	double T_s;
      	InputType1 error_jp, error_1_jp;
      	InputType2 error_jv, error_1_jv;
      	unitless_type_jp intError_jp, intErrorLimit_jp;
      	unitless_type_jp kp, ki, kd;
      	OutputType controlSignal_pid, controlSignalLimit;

      	void getSamplePeriodFromEM();


   protected:

      virtual void operate();

//    public:
//      jp_type wamJP, wamjpref;
//      jv_type wamJV, wamjvref;
//      jt_type jtOutput;
    private:
      DISALLOW_COPY_AND_ASSIGN(CriticController);
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(MathTraits::RequiresAlignment)
//
};

}
}
#include <Detail/CriticController-inl.h>
#endif // CRITICCONTROLLER_H_
