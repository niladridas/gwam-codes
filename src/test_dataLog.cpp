#include <iostream>
#include <string>

#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <barrett/math.h>
#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
//#include "myDataProcessing.h"
#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>
#include <myDataProcessing.h>
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
	if (argc != 4) {
		std::cout << "Usage: " << argv[0] << " <File Write>"<< " <Filter1 cutoff freq.>"<<" <Filter2 cutoff freq.>" << std::endl;
		return false;
	}
	return true;
}




template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {





	const uint NS=4;
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	typedef typename systems::myDataProcessing <DOF,NS>::vecRearranged_type vecRearranged_type;
	typedef  typename systems::myDataProcessing <DOF,NS>::reducedState_type reducedState_type;
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
	double omega1, omega2;
	if(argc>3){
	 omega1 = boost::lexical_cast<double>(argv[2]);
	omega2 = boost::lexical_cast<double>(argv[3]);
	}
	else{
		omega1 =10;
		omega2 =10;
	}
	std::cout<<"Omega_P1: "<<omega1<<"  Omega_P2: "<<omega2<<std::endl;


	wam.gravityCompensate();
	std::cout<<"gravity is on"<<std::endl;
	jp_type jp(0.0), desp(0.0), zerop(0.0), jcut1(50.0), jcut2(50.0);

	printMenu();

	systems::FirstOrderFilter<jp_type> filter1;
	filter1.setLowPass(jp_type(omega1));
	connect(wam.jpOutput, filter1.input);

	systems::FirstOrderFilter<jp_type> filter2;
		filter2.setLowPass(jp_type(omega2));
		connect(wam.jpOutput, filter2.input);


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
//	systems::TupleGrouper<double, vecRearranged_type, reducedState_type, double > tg;
	systems::TupleGrouper<double, jp_type, jp_type, jt_type, jp_type > tg;


	connect(time.output, tg.template getInput<0>());
//	connect(dp.rearranged_StateOutput, tg.template getInput<1>());
//	connect(dp.stateOutput, tg.template getInput<2>());
//	connect(dp.costOutput, tg.template getInput<3>());
//	connect(dp.processed_jtSumOutput, tg.template getInput<3>());

	connect(filter1.output, tg.template getInput<1>());
	connect(filter2.output, tg.template getInput<2>());

	connect(wam.jpController.controlOutput, tg.template getInput<3>());
	connect(wam.jpOutput, tg.template getInput<4>());
	//typedef boost::tuple<double, jp_type, jv_type, jt_type> tuple_type;
//	typedef boost::tuple<double, jp_type> tuple_type;
//	typedef boost::tuple<double, jp_type> tuple_type;
	
	
//	typedef boost::tuple<double, vecRearranged_type, reducedState_type, double> tuple_type;
	typedef boost::tuple<double, jp_type, jp_type, jt_type, jp_type> tuple_type;
	
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);



			for (size_t i=0; i<30; i++){
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
