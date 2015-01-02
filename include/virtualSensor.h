#ifndef VIRTUALSENSOR_H_
#define VIRTUALSENSOR_H_

#include <eigen3/Eigen/Core>
#include <barrett/detail/ca_macro.h>
#include <barrett/math/traits.h>
#include <barrett/systems/abstract/execution_manager.h>
#include <system.h>

namespace barrett {
namespace systems {

template<size_t DOF, size_t NJOINTS, typename InputType = typename units::JointPositions<DOF>::type,
typename MathTraits = math::Traits<InputType> >
class virtualSensor : public System {

	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
public:

	typedef  typename units::JointPositions<NJOINTS>::type jp2DOF_type;
	typedef  typename units::JointVelocities<NJOINTS>::type jv2DOF_type;
	typedef  typename units::JointAccelerations<NJOINTS>::type ja2DOF_type;

	Input<InputType> jp7DOF;

	Output<jv_type> jv7DOF, jpErrDot7DOF;
	Output<ja_type> ja7DOF;

	Output<jp2DOF_type> jp2DOF;
	Output<jv2DOF_type> jv2DOF;
	Output<ja2DOF_type> ja2DOF;

protected:
	typename Output<jv_type>::Value* jv7DOFval, *jpErrDot7DOFval;
	typename Output<ja_type>::Value* ja7DOFval;
	typename Output<jp2DOF_type>::Value* jp2DOFval;
	typename Output<jv2DOF_type>::Value* jv2DOFval;
	typename Output<ja2DOF_type>::Value* ja2DOFval;

private:
	InputType jp7, jp7_pre;
	jv_type jv7, jv7_pre, jpEdot7;
	ja_type ja7;
	jp2DOF_type jp2;
	jv2DOF_type jv2;
	ja2DOF_type ja2;
	double TS;

public:
	explicit virtualSensor(const std::string& sysName = "virtualSensor");
	virtual ~virtualSensor(){
		this->mandatoryCleanUp();
	}

	double& getTS(){return TS;}
protected:
	virtual void operate();
	virtual void onExecutionManagerChanged() {
			System::onExecutionManagerChanged();  // First, call super
			getSamplePeriodFromEM();
		}
	void getSamplePeriodFromEM();
	void setSamplePeriod(double timeStep);





private:
	DISALLOW_COPY_AND_ASSIGN(virtualSensor);
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(MathTraits::RequiresAlignment);


};

}
}
#include <Detail/virtualSensor-inl.h>
#endif
