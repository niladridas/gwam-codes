#include <iostream>
#include <string>

#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()

#include <boost/tuple/tuple.hpp>

#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
//#include "myDataProcessing.h"
//#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>
#include <myController.h>
//#include <myDataProcessing.h>
#include<inputMatrix.h>
//#include "myDataProcessing.h"


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

void printMenu() {
	printf("Commands:\n");
	printf("  t  start training \n");
	printf("  z  Move to the home position\n");
	printf("  h  Move to the home position\n");
	printf("  i  Idle (release position/orientation constraints)\n");
	printf("  q  Quit\n");
}

bool validate_args(int argc, char** argv) {
	if (argc != 2) {
		printf("Usage: %s <fileName>\n", argv[0]);
		return false;
	}
	return true;
}


template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
	const uint NS=4;
	const uint NR=5;

	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
//	typedef typename systems::myDataProcessing <DOF,NS>::vecRearranged_type vecRearranged_type;
//	typedef  typename systems::myDataProcessing <DOF,NS>::reducedState_type reducedState_type;
//	typedef typename systems::inputMatrix<DOF, NS, NR>::inputMatrix_type inputMatrix_type;

//	typedef typename systems::inputMatrix<DOF, NS, NR>::gx_type gx_type;
	srand(time(NULL));
	float ThMax1, ThMax3;
	float ThMin1, ThMin3;
	ThMax1 = 1.3;
	ThMin1 = -1.5;
	ThMax3 = 3.0;
	ThMin3 = -0.6;
	double ts = pm.getExecutionManager()->getPeriod();
	std::cout<<"TS:"<<std::endl;
	std::string line;
	//bool going = true;
	//bool homed = false;

	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

	wam.gravityCompensate();
	std::cout<<"gravity is on"<<std::endl;
	//myController registration
	// Create our system - this needs our default PID settings from the configuration file
	const libconfig::Setting& setting = pm.getConfig().lookup(pm.getWamDefaultConfigPath());
	const int NSTATES =4;
	systems::myController<NSTATES, jp_type, jv_type, jt_type> criticon(setting["joint_position_control"]);

	  // Connect our feedback
	        // Create a low level wam object
	        systems::connect(wam.jpOutput, criticon.feedbackInput);
	        systems::connect(wam.jvOutput, criticon.wamJVIn);
//	        // Register our controller as the default PID controller for Joint Position Control
	        wam.supervisoryController.registerConversion(
	                        systems::makeIOConversion(criticon.referenceInput,
	                                        criticon.controlOutput));



	jp_type jp(0.0), desp(0.0), zerop(0.0);
	printMenu();
//	systems::inputMatrix<DOF,NS,NR> im;
	systems::Ramp time(pm.getExecutionManager(), 1.0/ts);
	std::cout<<"call dp myData constructor"<<std::endl;
//	systems::myDataProcessing<DOF, NS> dp(ts);
//	dp.setQweight(10);
//	std::cout<<"Q: "<<dp.getQweight()<<std::endl;
//	std::cout<<"Dimension: "<<dp.getDim()<<std::endl;
//	connect(wam.jpOutput, dp.jpInput);
//	connect(wam.jvOutput, dp.jvInput);
//	connect(wam.jpController.controlOutput, dp.jtInput);
//	connect(wam.jtSum.output, dp.jtSumInput);
	//systems::TupleGrouper<double, jp_type, jv_type, jt_type> tg;
	//systems::TupleGrouper<double, jp_type> tg;

	std::cout<<jp.transpose()<<std::endl;


//	systems::FirstOrderFilter<jv_type> hp3;
//	double omega_p = 100.0;
//		wam.jvFilter.setLowPass(jv_type(omega_p));
//		hp3.setHighPass(jp_type(omega_p), jp_type(omega_p));
//			pm.getExecutionManager()->startManaging(hp3);
//			systems::Gain<jv_type, double, ja_type> changeUnits2(1.0);
//			connect(wam.jvOutput, hp3.input);
//				connect(hp3.output, changeUnits2.input);

//	systems::TupleGrouper<double, vecRearranged_type, reducedState_type, double > tg;
//	systems::TupleGrouper<double, jp_type, jv_type, ja_type, jt_type > tg;

//				connect(wam.jpOutput, im.jPositions);
//				connect(wam.jvOutput, im.jVelocities);
//				connect(wam.jpController.controlOutput, im.jTorques);

//				systems::TupleGrouper<double, jp_type, jv_type, jt_type, gx_type> tg;
				systems::TupleGrouper<double, jp_type, jv_type, jt_type> tg;
				connect(time.output, tg.template getInput<0>());
//	connect(im.gxrow, tg.template getImput<1>());
//	connect(dp.rearranged_StateOutput, tg.template getInput<1>());
//	connect(dp.stateOutput, tg.template getInput<2>());
//	connect(dp.costOutput, tg.template getInput<3>());
//	connect(dp.processed_jtSumOutput, tg.template getInput<3>());

				connect(wam.jpOutput, tg.template getInput<1>());

				connect(wam.jvOutput, tg.template getInput<2>());

				connect(criticon.controlOutput, tg.template getInput<3>());

//				connect(im.gxrow, tg.template getInput<4>());
//	connect(im.gx, tg.template getInput<4>());


//	typedef boost::tuple<double, vecRearranged_type, reducedState_type, double> tuple_type;
//	typedef boost::tuple<double, jp_type, jv_type, ja_type, jt_type> tuple_type;
//				typedef boost::tuple<double, jp_type, jv_type, jt_type, gx_type> tuple_type;
				typedef boost::tuple<double, jp_type, jv_type, jt_type> tuple_type;

	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);



			for (size_t i=0; i<1; i++){
				std::cout<<"joint position: "<<jp<<std::endl;
				jp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
				jp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);
//				desp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
//				desp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);
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
				wam.moveTo(jp);
				usleep(1000000);
				time.setOutput(0.0);
				time.start();
				connect(tg.output, logger.input);

				wam.moveTo(desp);
				systems::disconnect(logger.input);
				time.stop();
				usleep(1000000);
				//waitForEnter();
			}
			logger.closeLog();
			 wam.moveHome();
			// Wait for the user to press Shift-idle
			pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);



			printf("Logging stopped.\n");

			log::Reader<tuple_type> lr(tmpFile);
			lr.exportCSV(argv[1]);
			printf("Output written to %s.\n", argv[1]);
			std::remove(tmpFile);


	return 0;
}
