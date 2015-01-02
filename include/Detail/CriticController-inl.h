#include <cassert>

#include <libconfig.h++>

#include <barrett/math/utils.h>
#include <barrett/thread/abstract/mutex.h>

namespace barrett {
namespace systems {


template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
CriticController<DOF, NSTATES, InputType1, InputType2, OutputType, MathTraits>::CriticController(const std::string& sysName ):
Controller<InputType1,OutputType>(sysName),
 wamJVIn(this), /*wamJVref(this),*/
 critic("param.txt", "Wmatrices.txt", "center.txt"), T_s(0.0), error_jp(0.0), error_1_jp(0.0),error_jv(0.0), error_1_jv(0.0),
 intError_jp(0.0),intErrorLimit_jp(0.0),kp(0.0), ki(0.0), kd(0.0), controlSignal_pid(0.0), controlSignalLimit(0.0)
{
	Q = allocate<double>(NSTATES*NSTATES);
	eye(Q, NSTATES);
	R = allocate<double>((NSTATES/2)*(NSTATES/2));
	Rinv = allocate<double>((NSTATES/2)*(NSTATES/2));
	eye(R,NSTATES/2);
	gTld = allocate<float>(NSTATES/2);
	delx = allocate<float>(NSTATES);
	getSamplePeriodFromEM();

}
//
template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
CriticController<DOF, NSTATES, InputType1, InputType2, OutputType, MathTraits>::CriticController(const libconfig::Setting& setting, const std::string& sysName):
Controller<InputType1, OutputType>(sysName),
wamJVIn(this), /*wamJVref(this),*/
critic("param.txt", "Wmatrices.txt", "center.txt"), T_s(0.0), error_jp(0.0), error_1_jp(0.0),error_jv(0.0), error_1_jv(0.0),
 intError_jp(0.0),intErrorLimit_jp(0.0),kp(0.0), ki(0.0), kd(0.0), controlSignal_pid(0.0), controlSignalLimit(0.0)
{
	Q = allocate<double>(NSTATES*NSTATES);
		eye(Q, NSTATES);
		R = allocate<double>((NSTATES/2)*(NSTATES/2));
		Rinv = allocate<double>((NSTATES/2)*(NSTATES/2));
		eye(R,NSTATES/2);
		gTld = allocate<float>(NSTATES/2);
		delx = allocate<float>(NSTATES);

		getSamplePeriodFromEM();
		setFromConfig(setting);
}
//
template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void CriticController<DOF, NSTATES, InputType1, InputType2, OutputType, MathTraits>::setFromConfig(const libconfig::Setting& setting){
	if (setting.exists("kp")) {
			setKp(unitless_type_jp(setting["kp"]));
		} else {
			setKp(unitless_type_jp(0.0));
		}
		if (setting.exists("ki")) {
			setKi(unitless_type_jp(setting["ki"]));
		} else {
			setKi(unitless_type_jp(0.0));
		}
		if (setting.exists("kd")) {
			setKd(unitless_type_jp(setting["kd"]));
		} else {
			setKd(unitless_type_jp(0.0));
		}
		if (setting.exists("integrator_limit")) {
			setIntegratorLimit(unitless_type_jp(setting["integrator_limit"]));
		} else {
			setIntegratorLimit(unitless_type_jp(0.0));
		}
		if (setting.exists("control_signal_limit")) {
			setControlSignalLimit(unitless_type_jp(setting["control_signal_limit"]));
		} else {
			setControlSignalLimit(unitless_type_jp(0.0));
		}
	}
//
//
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF, NSTATES, InputType1, InputType2, OutputType, MathTraits>::setSamplePeriod(double timeStep)
		{
			T_s = timeStep;
		}
//
//
//
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::setKp(const unitless_type_jp& proportionalGains){
			kp = proportionalGains;
		}
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::setKi(const unitless_type_jp& integralGains){
			ki = integralGains;
		}
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::setKd(const unitless_type_jp& derivitiveGains){
			kd=derivitiveGains;
		}
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::setIntegratorState(const unitless_type_jp& integratorState){
			// intError is written and read in operate(), so it needs to be locked.
				BARRETT_SCOPED_LOCK(this->getEmMutex());
				intError_jp = integratorState;
		}

		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::setIntegratorLimit(const unitless_type_jp& intSaturations){
			intErrorLimit_jp = intSaturations;
		}

		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::setControlSignalLimit(const unitless_type_jp& csSaturations){
			controlSignalLimit = csSaturations;
		}
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		inline void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::resetIntegrator(){
			setIntegratorState(unitless_type_jp(0.0));
		}
//
		template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF,NSTATES, InputType1, InputType2, OutputType, MathTraits>::operate()
		{
			typedef MathTraits MT;

			error_jp = MT::sub(this->referenceInput.getValue(), this->feedbackInput.getValue());

			intError_jp = MT::add(intError_jp, MT::mult(ki, MT::mult(T_s, error_1_jp)));
			if (intErrorLimit_jp != MT::zero()) {
				intError_jp = math::saturate(intError_jp, intErrorLimit_jp);
			}

			controlSignal_pid = MT::add(MT::mult(kp, error_jp),
									MT::add(intError_jp,
										MT::mult(kd, MT::div(MT::sub(error_jp, error_1_jp), T_s))));
			if (controlSignalLimit != MT::zero()) {
				controlSignal_pid = math::saturate(controlSignal_pid, controlSignalLimit);
			}

			error_1_jp = error_jp;

			this->controlOutputValue->setData(&controlSignal_pid);

		}

//
//
template<size_t DOF, int NSTATES, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
		void CriticController<DOF, NSTATES, InputType1, InputType2, OutputType, MathTraits>::getSamplePeriodFromEM()
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
