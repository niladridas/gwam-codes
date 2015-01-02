#include <iostream>
#include <string>

#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/stl_utils.h>
#include <GMMTrajectory.h>
#include <barrett/log.h>
//#include <massMatrixHolder.hpp>

#include <barrett/standard_main_function.h>


using namespace barrett;


template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	systems::TupleGrouper<double, cp_type, cv_type> tg;
	char tmpFile[] = "/tmp/btXXXXXX";
	  if (mkstemp(tmpFile) == -1) {
	    printf("ERROR: Couldn't create temporary file!\n");
	    return 1;
	  }


	wam.gravityCompensate();

	jp_type startPos(0.0);
	startPos[1]= 0.32;
	startPos[3]= 3;
	cp_type xi;
	xi[0]=0.697;
	xi[1]= - 0.012;
	xi[2]= - 0.19;

	systems::Ramp time(pm.getExecutionManager(), 1.0);


	wam.moveTo(startPos);
	printf("Press [Enter] to follow traj.");
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<boost::tuple<double, cp_type, cv_type> > logger(pm.getExecutionManager(),
	    new log::RealTimeWriter<boost::tuple<double, cp_type, cv_type> > (tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
	    PERIOD_MULTIPLIER);
	Sam::waitForEnter();
	isl::GMMTrajectory gt(pm.getExecutionManager(), 0.01,0.01);
	systems::connect(wam.toolPosition.output,gt.cp_input);
	systems::connect(time.output, tg.template getInput<0>());
	systems::connect(gt.cp_output, tg.template getInput<1>());
	systems::connect(gt.cv_output, tg.template getInput<2>());
	time.start();
	pm.getExecutionManager()->startManaging(gt);

	systems::connect(tg.output, logger.input);
	btsleep(.01);
	wam.trackReferenceSignal(gt.cp_output);
	size_t cnt=0;
	while((xi-wam.getToolPosition()).norm()>0.1){
	    btsleep(0.00001);
	    cnt++;
	    if(cnt>2000){
	      break;
	    }
	}
	logger.closeLog();
	wam.moveHome();
	wam.idle();




//
//	printf("Press [Enter] to stop.");
//	Sam::waitForEnter();
//	wam.idle();

	// Wait for the user to press Shift-idle
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	log::Reader<boost::tuple<double, cp_type, cv_type> > lr(tmpFile);
	lr.exportCSV(argv[1]);
	printf("Output written to %s.\n", argv[1]);
	std::remove(tmpFile);
	return 0;
}
