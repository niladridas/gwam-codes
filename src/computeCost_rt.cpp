//#include"QuadTS.h"
//#include"system.h"
#define PI 3.14592653
#include <iostream>
#include <string>

#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/stl_utils.h>
#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()

#include <boost/tuple/tuple.hpp>

#include <barrett/log.h>

#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>


using namespace barrett;
using detail::waitForEnter;
using systems::connect;
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
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	srand(time(NULL));
	float ThMax1, ThMax3;
	float ThMin1, ThMin3;
	ThMax1 = 1.5;
	ThMin1 = -1.5;
	ThMax3 = 3.1;
	ThMin3 = -0.9;


	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
		}

	wam.gravityCompensate();
	jp_type jp(0.0), desp(0.0), zerop(0.0);

	systems::Ramp time(pm.getExecutionManager(), 1.0);

	printMenu();


	std::string line;
	bool going = true;
	bool homed = false;

	systems::TupleGrouper<double, jp_type, jv_type, jt_type> tg;
		connect(time.output, tg.template getInput<0>());
		connect(wam.jpOutput, tg.template getInput<1>());
		connect(wam.jvOutput, tg.template getInput<2>());
		//connect(wam.jtSum.output, tg.template getInput<3>());
		connect(wam.jpController.controlOutput, tg.template getInput<4>());
		typedef boost::tuple<double, jp_type, jv_type, jt_type> tuple_type;
			const size_t PERIOD_MULTIPLIER = 1;
			systems::PeriodicDataLogger<tuple_type> logger(
					pm.getExecutionManager(),
					new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
					PERIOD_MULTIPLIER);
			time.start();

	while (going) {

		printf(">>> ");
		std::getline(std::cin, line);
		switch (line[0]) {
				case 't':
					for (size_t i=0; i<2; i++){
						jp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
						jp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);
						std::cout<<"moving to the position: "<<jp<<std::endl;
						//waitForEnter();
						wam.moveTo(jp);
						std::cout<<"done moving"<<std::endl;
						connect(tg.output, logger.input);
						wam.moveTo(desp);
						logger.closeLog();
						}
					log::Reader<tuple_type> lr(tmpFile);
					lr.exportCSV(argv[1]);
					printf("Output written to %s.\n", argv[1]);
					std::remove(tmpFile);

					break;


				case 'z':
					std::cout<<"moving to the position: "<<zerop<<std::endl;
					waitForEnter();
					wam.moveTo(zerop);
					break;
//				case 'p':
//					moveToStr(wam, &cp, "tool position", line.substr(1));
//					break;

				case 'h':
					std::cout << "Moving to home position: "
							<< wam.getHomePosition() << std::endl;
					wam.moveHome();
					homed = true;
					break;

				case 'i':
					printf("WAM idled.\n");
					wam.idle();
					break;

				case 'q':
				case 'x':
					printf("Quitting.\n");
					going = false;
					break;

				default:
					if (line.size() != 0) {
						printf("Unrecognized option.\n");
						printMenu();
					}
					break;
				}
				if(!homed){
					std::cout<<"get back to desired position?"<<std::endl;
					waitForEnter();
					wam.moveTo(desp);
				}
			}


			wam.idle();
			pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
			return 0;


}
