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
protected:
    typename Output<vecRearranged_type >::Value* rearrStateOutputValue;
    typename Output<double>::Value* costOutputValue;
    typename Output<jt_type>::Value* prosJtOutputValue;
    typename Output<jp_type>::Value* prosJpOutputValue;
    typename Output<jv_type>::Value* prosJvOutputValue;
    //typename Output<jt_type>::Value* prosJtSumOutputValue;
private:

    jp_type wamJP;
    jv_type wamJV;
    jt_type jtCurrent;
    //jt_type jtSumCurrent;

    Eigen::Matrix4d Q;
    Eigen::Matrix2d R;
    Eigen::Vector4d systemStates;
    Eigen::Vector2d torque;
    double TS_;

public:
    myDataProcessing(double ts=1.0,const std::string& sysName = "myDataProcessing") :
        System(sysName), jpInput(this), jvInput(this), jtInput(this), rearranged_StateOutput(this, &rearrStateOutputValue),
        costOutput(this, &costOutputValue), processed_jtOutput(this, &prosJtOutputValue), jpOutput(this, &prosJpOutputValue),jvOutput(this, &prosJvOutputValue),
       /* processed_jtSumOutput(this, &prosJtSumOutputValue),*/Q(Eigen::Matrix4d()), R(Eigen::Matrix2d()), systemStates(Eigen::Vector4d()), torque(Eigen::Vector2d()),
       TS_(ts), rearrStateOutput(0.0),
       Cost(0.0), prosJtOutput(0.0), jpo(0.0), jvo(0.0)/*, prosJtSumOutput(0.0)*/

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
      jpo=wamJP;
      jvo=wamJV;
      jtCurrent = jtInput.getValue();
      prosJtOutput = jtCurrent;
      systemStates(0)=wamJP[1];
      systemStates(1)=wamJV[1];
      systemStates(2)=wamJP[3];
      systemStates(3)=wamJV[3];
      torque(0)=jtCurrent[1];
      torque(1)=jtCurrent[3];
    //  tmp(0,0)= systemStates.transpose()*Q*systemStates [0];
      Cost =  TS_*((systemStates.transpose()*Q*systemStates + torque.transpose()*R*torque)(0)) ;

      rearrStateOutputValue->setData(&rearrStateOutput);
      costOutputValue->setData(&Cost);
      prosJtOutputValue->setData(&prosJtOutput);
      prosJpOutputValue->setData(&jpo);
      prosJvOutputValue->setData(&jvo);
//      prosJtSumOutputValue->setData(&prosJtSumOutput);
    }

public:
    vecRearranged_type rearrStateOutput;
    double Cost;
    jt_type prosJtOutput;
    jp_type jpo;
    jv_type jvo;

  private:
    DISALLOW_COPY_AND_ASSIGN(myDataProcessing);
};


 }
 }




#endif
