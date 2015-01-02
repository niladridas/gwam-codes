/*
 * myWAM.h
 *
 *  Created on: 02-Oct-2013
 *      Author: mobman
 */

#ifndef MYWAM_H_
#define MYWAM_H_
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()
#include <sys/time.h>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include<errDotGen.h>
#include <myController.h>
#include <compTorqueController.h>
#include<BmatrixGenerator.hpp>
#include <massMatrixHolder.hpp>

/* Auxiliary System (Computes M!) */
#include <dp_DynamicsBaseSystem.hpp>
//#include <virtualSensor.h>
namespace barrett {

template <typename InputType, typename controllerType, typename tgType, size_t DOF>
class myUtilWAMwrapper {

public:
	typedef typename math::Traits<InputType>::unitless_type unitless_type;
private:
	controllerType* controllerPtr;
	systems::Wam<DOF>* wamPtr;
	tgType* tgPtr;
	systems::Ramp* timePtr;

public:
	explicit myUtilWAMwrapper(controllerType* controller,  systems::Wam<DOF>* wam, tgType* tg, systems::Ramp* time);
	void myWamMoveTo(const InputType& desp, const double& mFact=3.0, const double& maxVel =1.5);
	void MoveWAMto(const InputType& desp, const double& mFact=3.0, const double& maxVel =1.5,const double TRANSITION_DURATION = 0.1);



};




}
#include<Detail/myUtilWAMwrapper-inl.h>
#endif /* MYWAM_H_ */
