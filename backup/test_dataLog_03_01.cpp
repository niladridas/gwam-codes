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
#include "../include/myDataProcessing.h"
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
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	typedef typename systems::myDataProcessing <DOF,NS>::vecRearranged_type vecRearranged_type;
	srand(time(NULL));
	float ThMax1, ThMax3;
	float ThMin1, ThMin3;
	ThMax1 = 1.5;
	ThMin1 = -1.5;
	ThMax3 = 3.1;
	ThMin3 = -0.9;
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
	jp_type jp(0.0), desp(0.0), zerop(0.0);
	printMenu();


	systems::Ramp time(pm.getExecutionManager(), 1.0/ts);
	std::cout<<"call dp myData constructor"<<std::endl;
	systems::myDataProcessing<DOF, NS> dp(ts);
	std::cout<<"R: "<<dp.getRweight()<<std::endl;
	std::cout<<"Dimension: "<<dp.getDim()<<std::endl;
	connect(wam.jpOutput, dp.jpInput);
	connect(wam.jvOutput, dp.jvInput);
	connect(wam.jpController.controlOutput, dp.jtInput);
//	connect(wam.jtSum.output, dp.jtSumInput);
	//systems::TupleGrouper<double, jp_type, jv_type, jt_type> tg;
	//systems::TupleGrouper<double, jp_type> tg;

	std::cout<<jp.transpose()<<std::endl;
	systems::TupleGrouper<double, vecRearranged_type, jp_type, jv_type, jt_type, double > tg;
//	systems::TupleGrouper<double, vecJp_type, jt_type, jt_type > tg;


	connect(time.output, tg.template getInput<0>());
	connect(dp.rearranged_StateOutput, tg.template getInput<1>());
	connect(dp.jpOutput, tg.template getInput<2>());
	connect(dp.jvOutput, tg.template getInput<3>());
	connect(dp.processed_jtOutput, tg.template getInput<4>());
	connect(dp.costOutput, tg.template getInput<5>());
//	connect(dp.processed_jtSumOutput, tg.template getInput<3>());

	//connect(wam.jvOutput, tg.template getInput<2>());
	//connect(wam.jtSum.output, tg.template getInput<3>());
	//connect(wam.jpController.controlOutput, tg.template getInput<3>());

	//typedef boost::tuple<double, jp_type, jv_type, jt_type> tuple_type;
//	typedef boost::tuple<double, jp_type> tuple_type;
//	typedef boost::tuple<double, jp_type> tuple_type;
	
	
	typedef boost::tuple<double, vecRearranged_type, jp_type, jv_type, jt_type, double> tuple_type;
//	typedef boost::tuple<double, vecJp_type, jt_type, jt_type> tuple_type;
	
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);



			for (size_t i=0; i<3; i++){
				std::cout<<"joint position: "<<jp<<std::endl;
				jp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
				jp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);
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
