/*
 * myUtilWAMwrapper-inl.h
 *
 *  Created on: 03-Oct-2013
 *      Author: mobman
 */

#ifndef MYUTILWAMWRAPPER_INL_H_
#define MYUTILWAMWRAPPER_INL_H_
#include<iostream>
#include <barrett/thread/abstract/mutex.h>
namespace barrett{
template<typename InputType, typename controllerType, typename tgType, size_t DOF>
myUtilWAMwrapper<InputType, controllerType, tgType, DOF>::myUtilWAMwrapper(controllerType* controller,  systems::Wam<DOF>* wam, tgType* tg, systems::Ramp* time)
{



	controllerPtr = controller;
	wamPtr = wam;
	tgPtr = tg;
	timePtr = time;

}



template<typename InputType, typename controllerType, typename tgType, size_t DOF>
void myUtilWAMwrapper<InputType, controllerType, tgType, DOF>::myWamMoveTo(const InputType& desp, const double& mFact, const double& maxVel){

		InputType sp(wamPtr->getJointPositions());


		systems::errDotGen<InputType> edg(sp, desp, mFact, maxVel);
		systems::connect(wamPtr->jpOutput, edg.feedbackInput);
		systems::connect(timePtr->output, edg.timeInput);

//
		systems::disconnect(controllerPtr->referenceInput);
		systems::connect(edg.refOut, controllerPtr->referenceInput);
		systems::connect(edg.errorOut, controllerPtr->posErrorIn);

		systems::connect(edg.errorDotOut, controllerPtr->posErrorDotIn);

		systems::connect(edg.refAcc, controllerPtr->refAccIn);
//		connect(edg.refAcc, tgPtr->template getInput<1>());
//		connect(edg1.timeOutput, tgPtr->template getInput<2>());
//		connect(edg.errorOut,tgPtr->template getInput<4>());
//		connect(edg.errorDotOut,tgPtr->template getInput<5>());
		controllerPtr->resetRcount();
		timePtr->start();
		wamPtr->moveTo(desp);
		systems::disconnect(edg.timeInput);


		systems::disconnect(controllerPtr->posErrorIn);
		systems::disconnect(controllerPtr->posErrorDotIn);
		systems::disconnect(controllerPtr->refAccIn);
//		disconnect(tgPtr->template getInput<1>());

	}




template<typename InputType, typename controllerType, typename tgType, size_t DOF>
void myUtilWAMwrapper<InputType, controllerType, tgType, DOF>::MoveWAMto(const InputType& desp, const double& mFact, const double& maxVel,const double TRANSITION_DURATION){


  timePtr->setOutput(0.0);
  InputType sp(wamPtr->getJointPositions());
  std::cout<<"current position of WAM: "<<sp<<std::endl;
  std::cout<<"Desired position of WAM: "<<desp<<std::endl;

  systems::errDotGen<InputType> edg(sp, desp, mFact, maxVel);
  systems::connect(wamPtr->jpOutput, edg.feedbackInput);
  systems::connect(timePtr->output, edg.timeInput);
  systems::disconnect(controllerPtr->referenceInput);
  systems::connect(edg.refOut, controllerPtr->referenceInput);
  systems::connect(edg.errorOut, controllerPtr->posErrorIn);
  systems::connect(edg.errorDotOut, controllerPtr->posErrorDotIn);
  systems::connect(edg.refAcc, controllerPtr->refAccIn);
  controllerPtr->resetRcount();
  wamPtr->supervisoryController.connectInputTo(controllerPtr->controlOutput);
//  wamPtr->trackReferenceSignal(controllerPtr->controlOutput);
//  systems::forceConnect(controllerPtr->controlOutput,wamPtr->input);
  timePtr->smoothStart(TRANSITION_DURATION);
//  timePtr->start();


  printf("Press [Enter] to stop.");
//  detail::waitForEnter();
  while((desp-wamPtr->getJointPositions()).norm()>0.01){
    btsleep(0.00001);
  }
  timePtr->smoothStop(TRANSITION_DURATION);
//  timePtr->stop();
  wamPtr->idle();

  systems::disconnect(edg.timeInput);
  systems:: disconnect(controllerPtr->posErrorIn);
  systems::disconnect(controllerPtr->posErrorDotIn);
  systems::disconnect(controllerPtr->refAccIn);
}

}



#endif /* MYUTILWAMWRAPPER_INL_H_ */
