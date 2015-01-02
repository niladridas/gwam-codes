/*
 * GMMTrajectory.h
 *
 *  Created on: 18-Dec-2014
 *      Author: mobman
 */

#ifndef GMMTRAJECTORY_H_
#define GMMTRAJECTORY_H_
//#include <eigen3/Eigen/Core>
#include <barrett/units.h>
#include <barrett/systems.h>

#include <barrett/detail/ca_macro.h>
#include <barrett/math/traits.h>
#include <barrett/systems/abstract/system.h>
//#include <barrett/systems/abstract/single_io.h>
#include <barrett/systems/abstract/execution_manager.h>
#include <GMM.h>
using namespace barrett;


namespace isl{

class GMMTrajectory: public systems::System {

	BARRETT_UNITS_FIXED_SIZE_TYPEDEFS;



public:
	explicit GMMTrajectory(systems::ExecutionManager* em,double Rho_0, double K,const std::string& sysName="GMMTrajectory");
	virtual ~GMMTrajectory();
	typedef Eigen::VectorXd vec_3D;

public:
	Input<cp_type> cp_input;
	Output<cp_type> cp_output;
	Output<cv_type> cv_output;
protected:
	Output<cp_type>::Value* cp_out_val;
    Output<cv_type>::Value* cv_out_val;

private:
	GMM gmm;
	cp_type cp_tmp;
	cv_type cv_tmp;
	vec_3D cp;
	double Ts, gain;

protected:
	void operate();
public:
	void start(const cp_type& curr_cp);
	void stop();



private:
	DISALLOW_COPY_AND_ASSIGN(GMMTrajectory);
};
}
#include<Detail/GMMTrajectory-inl.h>
#endif /* GMMTRAJECTORY_H_ */
