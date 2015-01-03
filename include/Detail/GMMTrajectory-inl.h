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
		System(sysName), cp_input(this), cp_output(this, &cp_out_val), cv_output(this, &cv_out_val), gmm_out(this, &gmm_out_val), diff_norm_out(this, &diff_norm_out_val),
		gmm(Rho_0,K),cvFilter(em->getPeriod()), cp_tmp(0.0), cv_tmp(0.0),gmm_cv_tmp(0.0),cp_ref(vecXD()), cp_real(vecXD()),Ts(0.0),gain(0.0), nrm(0.0), difnrm(0.0)
{

//	if (em != NULL){
//		em->startManaging(*this);
//	}
	cp_ref.resize(3,1);
	cp_real.resize(3,1);

	Ts = em->getPeriod();
	cvFilter.setLowPass(cv_type(50.0));
}

GMMTrajectory::~GMMTrajectory() {
	this->mandatoryCleanUp();
	// TODO Auto-generated destructor stub
}
void GMMTrajectory::start(const cp_type& curr_cp){
  BARRETT_SCOPED_LOCK(getEmMutex());
  cp_ref = curr_cp;
  gain=1.0;
}

void GMMTrajectory::stop(){
  BARRETT_SCOPED_LOCK(getEmMutex());
  gain=0.0;
}

void GMMTrajectory::setLowpass(double omega){
  cvFilter.setLowPass(cv_type(omega));
}

void GMMTrajectory::operate(){


	cp_real = cp_input.getValue();
	gmm.calculate_model_output(cp_ref);

	nrm = gmm.output.norm();
	for(size_t i=0; i<3; ++i){
	  gmm_cv_tmp[i]=gmm.output(i);
	}
	gmm_out_val->setData(&gmm_cv_tmp);

	if(nrm>MAX_VEL){
	  gmm.output = gmm.output*(MAX_VEL/nrm);
	}
//	cv_tmp = gmm.output;
	cv_tmp = cvFilter.eval(gmm.output);

	cp_ref=cp_ref + Ts*cv_tmp*gain;
	difnrm=(cp_ref-cp_real).norm();
	if(difnrm>0.02){
	  cp_ref = cp_real;
	}
	diff_norm_out_val->setData(&difnrm);
	cp_tmp = cp_ref;
//	cv_tmp = gmm.output;
	cp_out_val->setData(&cp_tmp);
	cv_out_val->setData(&cv_tmp);


}

}
