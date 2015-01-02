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
#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>
#include <virtualSensor.h>
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
	ThMax1 = 1.3;
	ThMin1 = -1.5;
	ThMax3 = 3.0;
	ThMin3 = -0.6;
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
	systems::FirstOrderFilter<jp_type> filter;
	filter.setLowPass(jp_type(omega));
	connect(wam.jpOutput, filter.input);
	const size_t NSTATES = 4;

	typedef  typename systems::virtualSensor<DOF, NSTATES/2>::jp2DOF_type jp2DOF_type;
	typedef  typename systems::virtualSensor<DOF, NSTATES/2>::jv2DOF_type jv2DOF_type;
	typedef  typename systems::virtualSensor<DOF, NSTATES/2>::ja2DOF_type ja2DOF_type;
	systems::virtualSensor<DOF, NSTATES/2> stateEstimator;

	std::cout<<"TS FROM STATEESTIMATOR: "<<stateEstimator.getTS()<<std::endl;

	connect(filter.output, stateEstimator.jp7DOF);

	systems::Ramp time(pm.getExecutionManager(), 1.0/ts);
	systems::TupleGrouper<double, jp_type, jv_type, jp2DOF_type, jv2DOF_type, ja2DOF_type > tg;
	connect(time.output, tg.template getInput<0>());
//						connect(filter.output, tg.template getInput<1>()); // comented because collecting unfiltered data
	connect(wam.jpOutput, tg.template getInput<1>());
	connect(wam.jvOutput, tg.template getInput<2>());
	connect(stateEstimator.jp2DOF, tg.template getInput<3>());
	connect(stateEstimator.jv2DOF, tg.template getInput<4>());
	connect(stateEstimator.ja2DOF, tg.template getInput<5>());





	typedef boost::tuple<double, jp_type, jv_type, jp2DOF_type, jv2DOF_type, ja2DOF_type> tuple_type;
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);
	for (size_t i=0; i<numLoop; i++){
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
