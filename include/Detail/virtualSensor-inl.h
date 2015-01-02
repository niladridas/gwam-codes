#include <barrett/math/utils.h>
namespace barrett {
namespace systems {


template<size_t DOF, size_t NJOINTS, typename InputType, typename MathTraits>
virtualSensor<DOF, NJOINTS, InputType, MathTraits>::virtualSensor(const std::string& sysName) :
System(sysName), jp7DOF(this), jv7DOF(this, &jv7DOFval), jpErrDot7DOF(this, &jpErrDot7DOFval), ja7DOF(this, &ja7DOFval), jp2DOF(this, &jp2DOFval), jv2DOF(this, &jv2DOFval), ja2DOF(this, &ja2DOFval),
jp7(0.0), jp7_pre(0.0), jv7(0.0), jv7_pre(0.0), ja7(0.0), jp2(0.0), jv2(0.0), ja2(0.0), TS(0.0)
{
	getSamplePeriodFromEM();
}

template<size_t DOF, size_t NJOINTS, typename InputType, typename MathTraits>
void virtualSensor<DOF, NJOINTS, InputType, MathTraits>::setSamplePeriod(double timeStep)
{
	TS = timeStep;
}

template<size_t DOF, size_t NJOINTS, typename InputType, typename MathTraits>
void virtualSensor<DOF, NJOINTS, InputType, MathTraits>::getSamplePeriodFromEM()
{
	if (this->hasExecutionManager()) {
		assert(this->getExecutionManager()->getPeriod() > 0.0);
		setSamplePeriod(this->getExecutionManager()->getPeriod());
	} else {
		setSamplePeriod(0.0);
	}
}

template<size_t DOF, size_t NJOINTS, typename InputType, typename MathTraits >
void virtualSensor<DOF, NJOINTS, InputType, MathTraits >::operate(){

	typedef MathTraits MT;
	jp7 = jp7DOF.getValue();
	jv7=MT::div(MT::sub(jp7, jp7_pre), TS);
	ja7 = MT::div(MT::sub(jv7, jv7_pre), TS);

	jp2[0]=jp7[1];
	jp2[1]=jp7[3];

	jv2[0]=jv7[1];
	jv2[1]=jv7[3];

	ja2[0]=ja7[1];
	ja2[1]=ja7[3];

	jv7DOFval->setData(&jv7);
	ja7DOFval->setData(&ja7);

	jp2DOFval->setData(&jp2);
	jv2DOFval->setData(&jv2);
	ja2DOFval->setData(&ja2);

	jp7_pre = jp7;
	jv7_pre = jv7;
}




}
}
