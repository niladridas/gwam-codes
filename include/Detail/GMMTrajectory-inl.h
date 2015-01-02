/*
 * GMMTrajectory.cpp
 *
 *  Created on: 18-Dec-2014
 *      Author: mobman
 */

//#include <GMMTrajectory.h>
using namespace barrett;
namespace isl{

GMMTrajectory::GMMTrajectory(systems::ExecutionManager* em,double Rho_0, double K, const std::string& sysName):
		System(sysName), cp_input(this), cp_output(this, &cp_out_val), cv_output(this, &cv_out_val),	gmm(Rho_0,K), cp_tmp(0.0),  cv_tmp(0.0),cp(vec_3D()),Ts(0.0) {

//	if (em != NULL){
//		em->startManaging(*this);
//	}
	cp.resize(3,1);

	Ts = em->getPeriod();
}

GMMTrajectory::~GMMTrajectory() {
	this->mandatoryCleanUp();
	// TODO Auto-generated destructor stub
}

void GMMTrajectory::operate(){


	cp = cp_input.getValue();
	gmm.calculate_model_output(cp);

	cp=cp + Ts*gmm.output;
	cp_tmp = cp;
	cv_tmp = gmm.output;
	cp_out_val->setData(&cp_tmp);
	cv_out_val->setData(&cv_tmp);


}

}
