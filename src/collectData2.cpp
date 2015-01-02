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
#include<myUtilWAMwrapper.h>
#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>

using namespace barrett;
using systems::connect;
using detail::waitForEnter;



bool validate_args(int argc, char** argv) {
    if (argc <2) {
        std::cout << "Usage: " << argv[0] << " <File Write>"<< "<num. of Loop>"<<" <Filter1 cutoff freq.>"<< std::endl;
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
	float convertion = 0.0174532925;
	ThMax1 = 1.3;
	ThMin1 = -1.5;
	ThMax3 = 3.0;
	ThMin3 = -0.6;

	ThMax1 = 15 * convertion;
	ThMin1 = -15*convertion;
	ThMax3 = 15 * convertion;;
	ThMin3 = -15 * convertion;
	const int NSTATES =4;
	const int NR = 7;

	typedef systems::compTorqueController<NSTATES, NR, jp_type, jv_type, jt_type> compController_type;
	double ts = pm.getExecutionManager()->getPeriod();
	std::cout<<"TS:"<<std::endl;
	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

	double omega;
	size_t numLoop;
	if(argc==3){
	  numLoop = boost::lexical_cast<size_t>(argv[2]);

	}
	else{
	  numLoop =1;

	}

	if(argc==4){
	  omega = boost::lexical_cast<double>(argv[2]);

	        }
	else{
	  omega =100;

	}
	std::cout<<"Omega_P1: "<<omega<<"  num. of Loops: "<<numLoop<<std::endl;


	wam.gravityCompensate();
	jp_type jp(0.0), desp(0.0), zerop(0.0);

	//Define Controller
	const libconfig::Setting& setting = pm.getConfig().lookup(pm.getWamDefaultConfigPath());
	compController_type compCon(setting["joint_position_control"]);



	// group all joint variable
	typedef systems::TupleGrouper<double, jp_type, jv_type, ja_type, jt_type> tg_type;
	tg_type tg;

	//==============================================
	// Joint acc calculator
	systems::FirstOrderFilter<jv_type> filter;
	wam.jvFilter.setLowPass(jv_type(omega));
	filter.setHighPass(jp_type(omega), jp_type(omega));
	pm.getExecutionManager()->startManaging(filter);
	systems::Gain<jv_type, double, ja_type> changeUnits1(1.0);
	connect(wam.jvOutput, filter.input);
	connect(filter.output, changeUnits1.input);

	//=======================================================
	systems::connect(wam.jpOutput, compCon.feedbackInput);
	systems::connect(wam.jvFilter.output, compCon.wamJVIn);


	systems::Ramp time(pm.getExecutionManager(), 1.0/ts);

	connect(time.output, tg.template getInput<0>());

	connect(wam.jpOutput, tg.template getInput<1>());
	connect(wam.jvOutput, tg.template getInput<2>());
	connect(changeUnits1.output, tg.template getInput<3>());
	connect(compCon.controlOutput, tg.template getInput<4>());




	typedef boost::tuple<double, jp_type, jv_type, ja_type, jt_type> tuple_type;

	wam.supervisoryController.registerConversion(systems::makeIOConversion(compCon.referenceInput, compCon.controlOutput));


	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);
	std::string str(argv[3]);

	for (size_t i=0; i<2; i++){
					std::cout<<"joint position: "<<jp<<std::endl;
					jp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
					jp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);


//					if (str.compare("neq")==0){
//					if (i>8){
//						desp[1] = (ThMin1 + (ThMax1 - ThMin1) * float(rand())/RAND_MAX);
//						desp[3] = (ThMin3 + (ThMax3 - ThMin3) * float(rand())/RAND_MAX);
//					}
//					if(i>3 && i<8){
//						if(i==4){
//							jp[1]=ThMax1;
//							jp[3]=ThMax3;
//							desp[1]=ThMin1;
//							desp[3]=ThMin3;
//						}
//						if(i==5){
//							jp[1]=ThMax1;
//							jp[3]=ThMin3;
//							desp[1]=ThMin1;
//							desp[3]=ThMax3;
//						}
//						if(i==6){
//							jp[1]=ThMin1;
//							jp[3]=ThMax3;
//							desp[1]=ThMax1;
//							desp[3]=ThMin3;
//						}
//						if(i==7) {
//							jp[1]=ThMin1;
//							jp[3]=ThMin3;
//							desp[1]=ThMax1;
//							desp[3]=ThMax3;
//						}
//					}
//
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

				connect(tg.output, logger.input);
				time.start();

				myUtilWAMwrapper<jp_type, compController_type, tg_type,  DOF > wamRap2(&compCon, &wam, &tg,&time);
				wamRap2.myWamMoveTo(desp);
				systems::disconnect(logger.input);
				time.stop();
				usleep(1000000);
					//waitForEnter();
				}
				logger.closeLog();
				wam.supervisoryController.registerConversion(systems::makeIOConversion(wam.jpController.referenceInput, wam.jpController.controlOutput));
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
