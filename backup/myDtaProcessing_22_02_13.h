#ifndef MY_DATA_PROCESSING
#define MY_DATA_PROCESSING

#include<barrett/systems.h>
#include <barrett/systems/abstract/system.h>
#include <iostream>

 namespace barrett{
 namespace systems{

template<size_t DOF, uint NSTATES>
class myDataProcessing : public System
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	static const size_t dim= NSTATES*(NSTATES+1)/2;



public:
	typedef  typename units::JointPositions<dim>::type vecRearranged_type;
	typedef  typename units::JointPositions<NSTATES>::type reducedState_type;


public:
	Input<jp_type> jpInput;  //Wam jp o/p would be connected to it
	Input<jv_type> jvInput;	 //Wam jv o/p would be connected to it
	Input<jt_type> jtInput;	// Controller o/p would be connected to it
	//Input<jt_type> jtSumInput;

	Output<vecRearranged_type > rearranged_StateOutput;  //o/p available after rearrange of I/p states
	Output<double > costOutput;
	Output<jt_type > processed_jtOutput;
	Output<jp_type> jpOutput;
	Output<jv_type> jvOutput;
	Output<reducedState_type> stateOutput;
protected:
    typename Output<vecRearranged_type >::Value* rearrStateOutputValue;
    typename Output<double>::Value* costOutputValue;
    typename Output<jt_type>::Value* prosJtOutputValue;
    typename Output<jp_type>::Value* prosJpOutputValue;
    typename Output<jv_type>::Value* prosJvOutputValue;
    typename Output<reducedState_type>::Value* stateOutputValue;
    //typename Output<jt_type>::Value* prosJtSumOutputValue;
private:

    jp_type wamJP;
    jv_type wamJV;
    jt_type jtCurrent;
    //jt_type jtSumCurrent;
    vecRearranged_type rearrStateOutput;
    double Cost;
    Eigen::Matrix4d Q;
    Eigen::Matrix2d R;
    reducedState_type systemStates;
    Eigen::Vector2d torque;
    double TS_;

public:
    myDataProcessing(double ts=1.0,const std::string& sysName = "myDataProcessing") :
        System(sysName), jpInput(this), jvInput(this), jtInput(this), rearranged_StateOutput(this, &rearrStateOutputValue),
        costOutput(this, &costOutputValue), processed_jtOutput(this, &prosJtOutputValue), jpOutput(this, &prosJpOutputValue),jvOutput(this, &prosJvOutputValue),
       stateOutput(this, &stateOutputValue), rearrStateOutput(0.0), Cost(0.0), Q(Eigen::Matrix4d()), R(Eigen::Matrix2d()), systemStates(0.0), torque(Eigen::Vector2d()),
       TS_(ts)

    {
    	Q.setIdentity();
    	R=(R.setIdentity());
    }

    void setRweight(double rw){
    	R=rw*R;
    }
    Eigen::Matrix2d& getRweight(){
    	return R;
    }
    void setQweight(double qw){
        	Q=qw*Q;
        }
        Eigen::Matrix2d& getQweight(){
        	return Q;
        }

    virtual ~myDataProcessing()
    {
      this->mandatoryCleanUp();
    }

 const size_t getDim() const {
return dim;
}
protected:

    virtual void operate()
    {
      wamJP = jpInput.getValue();
      wamJV = jvInput.getValue();
      jtCurrent = jtInput.getValue();
      systemStates[0]=wamJP[1];
      systemStates[1]=wamJV[1];
      systemStates[2]=wamJP[3];
      systemStates[3]=wamJV[3];
      torque(0)=jtCurrent[1];
      torque(1)=jtCurrent[3];
    //  tmp(0,0)= systemStates.transpose()*Q*systemStates [0];
      Cost =  TS_*((systemStates.transpose()*Q*systemStates + torque.transpose()*R*torque)(0)) ;
      size_t i=0;

            for(size_t j=0; j<NSTATES; j++){
          		  for(size_t k=j; k<NSTATES;k++){
          			  if(j==k){
          				  rearrStateOutput[i]=systemStates(j)*systemStates(k);

          				  i++;

          			  }
          			  else{
          				  rearrStateOutput[i]=2*systemStates(j)*systemStates(k);
          				  i++;


          			  }
              	//Amat1(0,0)=1;
          		  }
      }



      rearrStateOutputValue->setData(&rearrStateOutput);
      costOutputValue->setData(&Cost);
      prosJtOutputValue->setData(&jtCurrent);
      prosJpOutputValue->setData(&wamJP);
      prosJvOutputValue->setData(&wamJV);
      stateOutputValue->setData(&systemStates);
//      prosJtSumOutputValue->setData(&prosJtSumOutput);
    }




private:
    DISALLOW_COPY_AND_ASSIGN(myDataProcessing);
};


 }
 }




#endif
