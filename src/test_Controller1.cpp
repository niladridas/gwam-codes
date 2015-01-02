#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()
#include <sys/time.h>
#include <boost/tuple/tuple.hpp>

#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
//#include "myDataProcessing.h"
//#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>
//#include "../include/myDataProcessing.h"
#include <myDataProcessing.h>
//#include <CriticController.h>

#include <myController.h>
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
	typedef  typename  systems::myDataProcessing <DOF,NS>::reducedState_type reducedState_type;
	srand(time(NULL));
	float ThMax1, ThMax3;
	float ThMin1, ThMin3;
	ThMax1 = 1.5;
	ThMin1 = -1.5;
	ThMax3 = 3.1;
	ThMin3 = -0.9;
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



	wam.gravityCompensate();

	// Create our system - this needs our default PID settings from the configuration file
	const libconfig::Setting& setting = pm.getConfig().lookup(pm.getWamDefaultConfigPath());
	const int NSTATES =4;
	systems::myController<NSTATES, jp_type, jv_type, jt_type> criticon(setting["joint_position_control"]);

	  // Connect our feedback
	        // Create a low level wam object
	        systems::connect(wam.jpOutput, criticon.feedbackInput);
	        systems::connect(wam.jvOutput, criticon.wamJVIn);
	        // Register our controller as the default PID controller for Joint Position Control
	        wam.supervisoryController.registerConversion(
	                        systems::makeIOConversion(criticon.referenceInput,
	                                        criticon.controlOutput));
//	        wam.supervisoryController.registerConversion(
//	        	                        systems::makeIOConversion(criticon.wamJVref,
//	        	                                        criticon.controlOutput));

//	jv_type oldKd = criticon.getKd();

//	jv_type newKd = oldKd / 2.0;

//	criticon.setKd(newKd);


	std::cout<<"gravity is on"<<std::endl;
	jp_type jp(0.0), desp(0.0), zerop(0.0);
	printMenu();


	systems::Ramp time(pm.getExecutionManager(), 1.0/ts); //slope 1/ts will give increment 1 for eash time step
	std::cout<<"call dp myData constructor"<<std::endl;
	systems::myDataProcessing<DOF, NS> dp(ts);
	std::cout<<"R: "<<dp.getRweight()<<std::endl;
	std::cout<<"Dimension: "<<dp.getDim()<<std::endl;
	connect(wam.jpOutput, dp.jpInput);
	connect(wam.jvOutput, dp.jvInput);
	connect(wam.jpController.controlOutput, dp.jtInput);
//	connect(wam.jtSum.output, dp.jtSumInput);

//	std::cout<<jp.transpose()<<std::endl;
	systems::TupleGrouper<double, vecRearranged_type, reducedState_type, double> tg;
//	systems::TupleGrouper<double, vecJp_type, jt_type, jt_type > tg;


	connect(time.output, tg.template getInput<0>());
	connect(dp.rearranged_StateOutput, tg.template getInput<1>());
	connect(dp.stateOutput, tg.template getInput<2>());
	connect(dp.costOutput, tg.template getInput<3>());



	
	
	typedef boost::tuple<double, vecRearranged_type, reducedState_type, double> tuple_type;
//	typedef boost::tuple<double, vecJp_type, jt_type, jt_type> tuple_type;
	
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);

			jp_type jp1(0.0), jp2(0.0), jp3(0.0), jp4(0.0);

			for (size_t i=0; i<1; i++){
				std::cout<<"joint position: "<<jp<<std::endl;
				jp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
				jp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);
				jp[1] =-0.7;
				jp[3]= 2.5;
				 jp1[1] = -1.96354;
                                jp1[3]= 0.5;
				jp2[1] = 0.0;
				jp2[3] = 3.1358;
				jp3[1] = 1.5;
				jp3[3] = 3.1358;
				jp4[1] = 1.5;
				jp4[3]= 0.0;
				std::cout<<"moving to: "<<jp<<std::endl;
				wam.moveTo(jp1);
				std::cout<<"moved to Pos: "<<jp<<std::endl;
				usleep(1000000);
				time.setOutput(0.0);
				time.start();
				connect(tg.output, logger.input);

				//wam.moveTo(desp);
				//wam.moveHome();
				systems::disconnect(logger.input);
				time.stop();
				usleep(1000000);
				//waitForEnter();
			}
			logger.closeLog();
			 wam.moveHome();
			// Wait for the user to press Shift-idle
//			jp[1] = -0.7;
//			jp[3]= 2.5;
//			std::cout<<"moving to: "<<jp<<std::endl;
//			wam.moveTo(jp);
			std::cout<<"moved to Pos: "<<jp<<std::endl;
			pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);


/*
///////////////////////////////////////////////////////////////////////////////
			//Time Calc
			jv_type jv(0.2);
			typedef math::Traits<jp_type> MathTraits;
			typedef MathTraits MT;
			jp_type currentJp,error_jp;
			jv_type currentJv;
			Sam::Dynamics<jp_type, jv_type> dyn;
			Eigen::Matrix<double,NSTATES,NSTATES/2> gx_wam;
			QuadTS critic("param.txt", "Wmatrices.txt", "center.txt");
			unsigned int NU =2;
			double* Rinv;
			double* gTld;
			double* gth;
			double* delx;
			double* err_states_wam;
			double* states_wam;
			double* controlSignal_critic;
			Rinv = allocate<double>((NU)*(NU));
			eye(Rinv,NU);
			gTld = allocate<double>(NU);
			gth = allocate<double>(NSTATES*(NU));
			delx = allocate<double>(NSTATES);
			err_states_wam = allocate<double>(NSTATES);
			states_wam = allocate<double>(NSTATES);
			controlSignal_critic = allocate<double>(NU);

			 timeval t1, t2;
			    double elapsedTime;

				std::cout<<"Tyme Calc Starts"<<std::endl;
//
				gettimeofday(&t1, NULL);
			currentJp = jp;//criticon.feedbackInput.getValue();
			error_jp = MT::sub(zerop, currentJp);
			currentJv = jv;//criticon.wamJVIn.getValue();
			dyn.computeModel(currentJp, currentJv);
			gx_wam=dyn.getGx();
			gth[0] = 0;
			gth[1] = 0;

			gth[2] = gx_wam(2,0);            //2nd row of gx
			gth[3] = gx_wam(2,1);


			gth[4] = 0;
			gth[5] = 0;


			gth[6] = gx_wam(3,0);
			gth[7] = gx_wam(3,1);

			states_wam[0] = currentJp[1];
			states_wam[1] = currentJv[1];
			states_wam[2] = currentJp[3];
			states_wam[3] = currentJv[3];

			err_states_wam[0]= -error_jp[1];
			err_states_wam[1]= states_wam[1]; // As ref velocity is 0
			err_states_wam[2]= -error_jp[3];
			err_states_wam[3]= states_wam[3];

			critic.computeDelX(err_states_wam, delx, states_wam);

			multiply(gth, delx, gTld, NU,1,NSTATES, TRANS);
			multiply(Rinv, gTld, controlSignal_critic, NU,1,NU);

			multiply(controlSignal_critic, double(-0.5), NU);
//			endT=clock();
			gettimeofday(&t2, NULL);

			elapsedTime = (t2.tv_sec - t1.tv_sec);
			std::cout<<"Time taken : "<<elapsedTime<<std::endl;
			del(Rinv);
			del(gTld);
			del(gth);
			del(delx);
			del(err_states_wam);
			del(states_wam);

			del(controlSignal_critic);

///////////////////////////////  End of time Calc/////////////////////////////////////
*/

	
			printf("Logging stopped.\n");

//			log::Reader<tuple_type> lr(tmpFile);
//			lr.exportCSV(argv[1]);
			printf("Output written to %s.\n", argv[1]);
			std::remove(tmpFile);
	

	return 0;
}


