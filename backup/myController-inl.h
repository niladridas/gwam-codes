#include <cassert>

#include <libconfig.h++>

#include <barrett/math/utils.h>
#include <barrett/thread/abstract/mutex.h>


namespace barrett {
namespace systems {
template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::myController(const std::string& sysName) :
	Controller<InputType1, OutputType>(sysName),
	wamJVIn(this), criticconOp(this, &criticconOpValue),/*gxrow(this, &gxrowVal), pEout(this, &pEoutVal), vEout(this, &vEoutVal), refout(this, &refoutVal),
	Rout1(this, &RoutVal1), Rout2(this, &RoutVal2),*/ ctorque(0.0), rCount(0), R(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()),
	Rinv(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), x(Eigen::Matrix<double, NSTATES, 1>()), err_states_wam(Eigen::Matrix<double, NSTATES, 1>()),
	u(Eigen::Matrix<double, NSTATES/2, 1>()), fuz("./FuzzParams/cons.txt","./FuzzParams/weight.txt", "./FuzzParams/model_gaussParams.txt"),

	gx_fuz(Eigen::Matrix<double,NSTATES,NSTATES/2>()),/* rout1(0.0),rout2(0.0), gx_matrow(0.0),*/T_s(0.0), error_jp(0.0), error_1_jp(0.0),currentJp(0.0), prevJp(0.0), error_jv(0.0), error_1_jv(0.0),currentJv(0.0),intError(0.0),
	intErrorLimit(0.0),	kp(0.0), ki(0.0), kd(0.0),controlSignal(0.0), controlSignalLimit(0.0)
	{
		getSamplePeriodFromEM();
		setR("./FuzzParams/Rmat.txt");
		Rinv=R.inverse();

		u.fill(0.0);
		x.fill(0.0);
		err_states_wam.fill(0.0);
		gx_fuz.fill(0.0);


	}



template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::myController(const libconfig::Setting& setting, const std::string& sysName) :
	Controller<InputType1, OutputType>(sysName), wamJVIn(this), criticconOp(this, &criticconOpValue),/*gxrow(this, &gxrowVal),
	pEout(this, &pEoutVal), vEout(this, &vEoutVal), refout(this, &refoutVal), Rout1(this, &RoutVal1), Rout2(this, &RoutVal2),*/
	ctorque(0.0),rCount(0),R(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()), Rinv(Eigen::Matrix<double, NSTATES/2, NSTATES/2>()),
	x(Eigen::Matrix<double, NSTATES, 1>()), err_states_wam(Eigen::Matrix<double, NSTATES, 1>()),u(Eigen::Matrix<double, NSTATES/2, 1>()),
	fuz("./FuzzParams/cons.txt","./FuzzParams/weight.txt", "./FuzzParams/model_gaussParams.txt"),
	gx_fuz(Eigen::Matrix<double,NSTATES,NSTATES/2>()),/* rout1(0.0),rout2(0.0),	gx_matrow(0.0),*/
	T_s(0.0), error_jp(0.0), error_1_jp(0.0), currentJp(0.0), prevJp(0.0), error_jv(0.0), error_1_jv(0.0),currentJv(0.0), intError(0.0),intErrorLimit(0.0),
	kp(0.0), ki(0.0), kd(0.0),controlSignal(0.0), controlSignalLimit(0.0)
{
	getSamplePeriodFromEM();
	setFromConfig(setting);
	setR("./FuzzParams/Rmat.txt");
	Rinv=R.inverse();
	u.fill(0.0);
	x.fill(0.0);
	err_states_wam.fill(0.0);
	gx_fuz.fill(0.0);
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setR(const char* RmatFile){

	std::vector<std::vector<double> > r;
	size_t tr, tc;

	Sam::readFile(RmatFile, r, tr, tc);
	if(!(tr==NSTATES/2 || tc== NSTATES/2)){

		std::cout<<"File "<<RmatFile<<" is not proper for R matrix"<<std::endl;
		exit(EXIT_FAILURE);
	}

	for (size_t i=0; i<NSTATES/2; i++){
		for (size_t j=0; j<NSTATES/2; j++){
			R(i,j)=r[i][j];
		}
	}

}



template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setFromConfig(const libconfig::Setting& setting)
{


  ///////////////////////////////////////////////////////////////
  // Test for myController
  //===========================================================
  if (setting.exists("kp")) {
		setKp(unitless_type_jp(setting["kp"]));
	}
	if (setting.exists("ki")) {
		setKi(unitless_type_jp(setting["ki"]));
	}
	if (setting.exists("kd")) {
		setKd(unitless_type_jp(setting["kd"]));
	}
	if (setting.exists("integrator_limit")) {
		setIntegratorLimit(unitless_type_jp(setting["integrator_limit"]));
	}
	
	//////////////////////////////////////////////////////////////////

	
	if (setting.exists("control_signal_limit")) {
		setControlSignalLimit(unitless_type_jp(setting["control_signal_limit"]));
	}
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setSamplePeriod(double timeStep)
{
	T_s = timeStep;
}


  
  ///////////////////////////////////////////////////////////////////////////////
  //Test for myController
  //============================================================================
  template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setKp(const unitless_type_jp& proportionalGains)
{
	std::cout << "Our previous proportional control gain: kp = " << kp << std::endl;
	kp = proportionalGains;
	std::cout << "Our new proportional control gain: kp = " << kp << std::endl;
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setKi(const unitless_type_jp& integralGains)
{
	ki = integralGains;
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setKd(const unitless_type_jp& derivitiveGains)
{
	kd = derivitiveGains;
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setIntegratorState(
		const unitless_type_jp& integratorState)
{
	// intError is written and read in operate(), so it needs to be locked.
	BARRETT_SCOPED_LOCK(this->getEmMutex());
	intError = integratorState;
}

template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setIntegratorLimit(
		const unitless_type_jp& intSaturations)
{
	intErrorLimit = intSaturations;
}


 template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
inline void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::resetIntegrator()
{
	setIntegratorState(unitless_type_jp(0.0));
}
  /////////////////////////////////////////////////////////////////////////////////////////



  

 template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::setControlSignalLimit(
		const unitless_type_jp& csSaturations)
{
	controlSignalLimit = csSaturations;
}


  
template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::operate()
{
	typedef MathTraits MT;
		currentJp = this->feedbackInput.getValue();
		error_jp = MT::sub(this->referenceInput.getValue(), currentJp);
		currentJv = wamJVIn.getValue();

		pE=error_jp;
		ref=this->referenceInput.getValue();
		vE=MT::div(MT::sub(error_jp, error_1_jp), T_s);
		//collect gx
//		gx_fuz = gx.getValue();
		///O/P gx row
//		gx_matrow(0,0)=gx_fuz(1,0);
//		gx_matrow(1,0)=gx_fuz(1,1);
//		gx_matrow(2,0)=gx_fuz(3,0);
//		gx_matrow(3,0)=gx_fuz(3,1);
//		gxrowVal->setData(&gx_matrow);

		// Original model

		x<<currentJp[1],currentJv[1],currentJp[3],currentJv[3];
		u<<controlSignal[1],controlSignal[3];

				fuz.computeFuzModel(x,u);





		gx_fuz=fuz.model_B;
//		rout1=(fuz.getR()).transpose();
//		RoutVal1->setData(&rout1);
//		size_t rc=0;
//		for(int i=0; i<NSTATES; i++){
//			for(int j=0; j<NSTATES/2; j++){
//				gx_matrow(rc,0)= gx_fuz(i,j);
//				rc++;
//			}
//		}

//		gxrowVal->setData(&gx_matrow);

	//////////////////////////////////////////////////////////////////
  // Critic Controller
	err_states_wam(0)=  /*currentJp[1];//*/-error_jp[1];
	err_states_wam(1)=  /*currentJv[1];//*/-vE[1]; // As e=curr-des
	err_states_wam(2)=  /*currentJp[3];//*/-error_jp[3];
	err_states_wam(3)= 	/*currentJv[3];//*/-vE[3];




	fuz.makeDelJDelx(err_states_wam,x);




//	rout2=(fuz.getR()).transpose();
//	RoutVal2->setData(&rout2);


  //=================================================================

	intError = MT::add(intError, MT::mult(ki, MT::mult(T_s, error_1_jp)));
	if (intErrorLimit != MT::zero()) {
		intError = math::saturate(intError, intErrorLimit);
	}

//	controlSignal = MT::add(MT::mult(kp, error_jp),
//							MT::add(intError,
//								MT::mult(kd, MT::div(MT::sub(error_jp, error_1_jp), T_s))));

	controlSignal = MT::add(MT::mult(kp, error_jp),
							MT::add(intError,
								MT::mult(kd, vE)));


	//===========================================================================================

	u = -0.5*Rinv*gx_fuz.transpose()*fuz.delVdelX;
//		u(0) = -kp[1]*err_states_wam(0) -kd[1]*err_states_wam(1) + intError[1];
//		u(1) = -kp[3]*err_states_wam(2) -kd[3]*err_states_wam(3) + intError[3];

		    ctorque[0] = u(0);
		    ctorque[1] =u(1);

	//=============================================================================================

//vE=intError;
//vEoutVal->setData(&vE);
//pEoutVal->setData(&pE);
//ref[0]=fuz.delVdelX(0);
//ref[1]=fuz.delVdelX(1);
//ref[2]=fuz.delVdelX(2);
//ref[3]=fuz.delVdelX(3);
//refoutVal->setData(&ref);
	prevJp=currentJp;
	error_1_jp = error_jp;

//	if(rCount>50){
//
//	controlSignal[1]=ctorque[0];
//	controlSignal[3]=ctorque[1];
//	}
	if (controlSignalLimit != MT::zero()) {
			controlSignal = math::saturate(controlSignal, controlSignalLimit);
		}
	this->controlOutputValue->setData(&controlSignal);
	criticconOpValue->setData(&ctorque);
	rCount++;

  /////////////////////////////////////////////////////////////////////////

  
}


  
  template<int NSTATES, int NR, typename InputType1, typename InputType2, typename OutputType, typename MathTraits>
  void myController<NSTATES, NR, InputType1, InputType2, OutputType, MathTraits>::getSamplePeriodFromEM()
{
	if (this->hasExecutionManager()) {
		assert(this->getExecutionManager()->getPeriod() > 0.0);
		setSamplePeriod(this->getExecutionManager()->getPeriod());
	} else {
		setSamplePeriod(0.0);
	}
}
}
}
