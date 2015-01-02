#define COMP_CONTROLLER
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
#define LOCAL_MASSMAT
#include<myUtilWAMwrapper.h>
#include<BmatrixGenerator.hpp>

//#include "myDataProcessing.h"
#define BARRETT_SMF_VALIDATE_ARGS

#include <barrett/standard_main_function.h>
//#include "../include/myDataProcessing.h"

//#include <CriticController.h>
//#include<errDotGen.h>
//#include <myController.h>
//#include <compTorqueController.h>
//#include <virtualSensor.h>
using namespace barrett;
using systems::connect;
using detail::waitForEnter;
// This function template will accept a math::Matrix with any number of rows,
// any number of columns, and any units. In other words: it will accept any
// barrett::units type.
template<int R, int C, typename Units>
bool parseDoubles(math::Matrix<R,C, Units>* dest, const std::string& str) {
	const char* cur = str.c_str();
	const char* next = cur;

	for (int i = 0; i < dest->size(); ++i) {
		(*dest)[i] = strtod(cur, (char**) &next);
		if (cur == next) {
			return false;
		} else {
			cur = next;
		}
	}

	// Make sure there are no extra numbers in the string.
	double ignore = strtod(cur, (char**) &next);
	(void)ignore;  // Prevent unused variable warnings

	if (cur != next) {
		return false;
	}

	return true;
}

template<size_t DOF, int R, int C, typename Units>
void moveToStr(systems::Wam<DOF>& wam, math::Matrix<R,C, Units>* dest,
		const std::string& description, const std::string& str)
{
	if (parseDoubles(dest, str)) {
		std::cout << "Moving to " << description << ": " << *dest << std::endl;
		wam.moveTo(*dest);
	} else {
		printf("ERROR: Please enter exactly %d numbers separated by "
				"whitespace.\n", dest->size());
	}
}



bool validate_args(int argc, char** argv) {
	if (argc <2) {
		std::cout << "Usage: " << argv[0] << " <File Write>"<< " <Filter1 cutoff freq.>"<< std::endl;
		return false;
	}
	return true;
}





template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
//	const uint NS=4;
	const int NSTATES =4;
	float convertion = 0.0174532925;


	const int NR = 9;
	double omega;
	double mFact, mFactUsr;
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);


	typedef typename systems::compTorqueController<NSTATES, NR, jp_type, jv_type, jt_type>::compTorque_type compTorque_type;
	typedef typename systems::compTorqueController<NSTATES, NR, jp_type, jv_type, jt_type>::unitless_type_jp compTorUless_type;
	typedef systems::TupleGrouper<double, jp_type, jv_type,  jt_type, compTorque_type> tg_type;
	typedef systems::compTorqueController<NSTATES, NR, jp_type, jv_type, jt_type> compController_type;
	srand(time(NULL));
	float ThMax1, ThMax3;
	float ThMin1, ThMin3;
	ThMax1 = 1.5;
	ThMin1 = -1.5;
	ThMax3 = 3.1;
	ThMin3 = -0.9;


	ThMax1 = 8 * convertion;
	ThMin1 = -8*convertion;
	ThMax3 = 8 * convertion;;
	ThMin3 = -8 * convertion;
	double ts = pm.getExecutionManager()->getPeriod();
	std::cout<<"TS:  "<<ts<<std::endl;
	std::string line;
	//bool going = true;
	//bool homed = false;
	
	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

	if(argc>2){
		omega = boost::lexical_cast<double>(argv[2]);

	}
	else{
		omega =100;

	}
	if(argc>3){
		mFactUsr = boost::lexical_cast<double>(argv[3]);
		mFact = mFactUsr;
	}
	else{
		mFactUsr = 4;
		mFact = mFactUsr;
	}
	std::cout<<"Omega_P1: "<<omega<<std::endl;

	wam.gravityCompensate();
	std::cout<<"gravity is on"<<std::endl;
	jp_type jp(0.0), desp(0.0), zerop(0.0);

	// Create our system - this needs our default PID settings from the configuration file
	const libconfig::Setting& setting = pm.getConfig().lookup(pm.getWamDefaultConfigPath());

	// MassMatrix System
	systems::BmatrixGenerator<DOF, NSTATES/2> mg(pm.getExecutionManager(),pm.getConfig().lookup(pm.getWamDefaultConfigPath()),
	    wam.getJointVelocities(),wam.getJointPositions());
	systems::connect(wam.jpOutput, mg.jpInput);
	systems::connect(wam.jvFilter.output, mg.jvInput);

	compController_type compCon(setting["joint_position_control"]);
	systems::connect(wam.jpOutput, compCon.feedbackInput);
	systems::connect(wam.jvFilter.output, compCon.wamJVIn);
	//==============================================

	systems::Ramp time(pm.getExecutionManager(), 1.0); //slope 1/ts will give increment 1 for eash time step

	tg_type tg;
	// Connect our feedback

	// Register our controller as the default PID controller for Joint Position Control
	 wam.supervisoryController.registerConversion(systems::makeIOConversion(compCon.referenceInput, compCon.controlOutput));

	connect(time.output, tg.template getInput<0>());
//	connect(compCon.refOut, tg.template getInput<1>());
	connect(wam.jpOutput,tg.template getInput<1>());
	connect(wam.jvOutput,tg.template getInput<2>());
	connect(compCon.controlOutput,tg.template getInput<3>());
	connect(compCon.compConOp, tg.template getInput<4>());
//	connect(compCon.pEout, tg.template getInput<5>());
//	connect(compCon.vEout, tg.template getInput<6>());
//	connect(compCon.refOut, tg.template getInput<7>());

	typedef boost::tuple<double, jp_type, jv_type,  jt_type, compTorque_type> tuple_type;




waitForEnter();
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(pm.getExecutionManager(),new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
	    PERIOD_MULTIPLIER);

	jp_type jp1(0.0), jp2(0.0), jp3(0.0), jp4(0.0);

	for (size_t i=0; i<1; i++){
	  if(i==0){
	    jp[1]=ThMin1;
	    jp[3]=ThMin3;
	  }
	  if(i==1){
	    jp[1]=ThMin1;
	    jp[3]=ThMax3;
	  }
	  if(i==2){
	    jp[1]=ThMax1;
	    jp[3]=ThMin3;
	  }
	  if(i==3) {
	    jp[1]=ThMax1;
	    jp[3]=ThMax3;
	  }
	  std::cout<<"moving to: "<<jp<<std::endl;
	  time.setOutput(0.0);
	  jp_type sp(wam.getJointPositions()), hp(wam.getHomePosition());
	  connect(tg.output, logger.input);
	  wam.moveTo(jp);
	  myUtilWAMwrapper<jp_type, compController_type, tg_type, DOF> wamRap2(&compCon, &wam, &tg, &time);
	  wamRap2.MoveWAMto(desp,3.0, 0.8);
	  std::cout << "moved to desp" << std::endl;
	  systems::disconnect(logger.input);
	  time.stop();
	  usleep(3000000);
	}
	logger.closeLog();
	std::cout<<"Logger closed"<<std::endl;
	//=======================================
	wam.supervisoryController.registerConversion(systems::makeIOConversion(wam.jpController.referenceInput, wam.jpController.controlOutput));
	//==============================================
	std::cout<<"Default controller in action"<<std::endl;

	wam.moveHome();
	// Wait for the user to press Shift-idle
//	    		jp[1] = -0.7;
//			jp[3]= 2.5;
//			std::cout<<"moving to: "<<jp<<std::endl;
//			wam.moveTo(jp);
	std::cout<<"moved to Pos: "<<jp<<std::endl;
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	printf("Logging stopped.\n");

	log::Reader<tuple_type> lr(tmpFile);
	lr.exportCSV(argv[1]);
	printf("Output written to %s.\n", argv[1]);
	std::remove(tmpFile);
	
	return 0;
}


