/*
 * costComputeCritic.hpp
 *
 *  Created on: 27-Oct-2014
 *      Author: mobman
 */

#ifndef COSTCOMPUTECRITIC_HPP_
#define COSTCOMPUTECRITIC_HPP_
#include <barrett/products/product_manager.h>
#include <barrett/detail/ca_macro.h>
#include <barrett/systems/abstract/system.h>
#include <barrett/systems.h>
#include <tsk.h>
using namespace barrett;
using namespace systems;

namespace Sam{
  template <size_t DOF, size_t NJOINTS, size_t NR>
  class costComputeCritic : public systems::System
  {

    BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
  public:
    Input<jp_type> jPosInp;
    Input<jv_type> jVelInp;
    Output<double> cost1Output, cost2Output;
  protected:
    typename Output<double>::Value* cost1OutputVal, *cost2OutputVal;

  private:
    jp_type jPos;
    jv_type jVel;
    Sam::tsk<NJOINTS*2, NJOINTS, NR> fuz1, fuz2;
    Eigen::Matrix<double,NJOINTS*2, 1> x;
    double cost1, cost2;

  public:
    costComputeCritic(systems::ExecutionManager* em,const std::string& w1Fname,const std::string& w2Fname,
        const std::string& consFname="./FuzzParams/cons.txt",const std::string& gaussFname="./FuzzParams/gaussParams.txt",
        const std::string& sysName = "costComputeCritic"):
          System(sysName), jPosInp(this), jVelInp(this), cost1Output(this,&cost1OutputVal), cost2Output(this, &cost2OutputVal), jPos(0.0), jVel(0.0),
          fuz1(consFname.c_str(),w1Fname.c_str(),gaussFname.c_str()),fuz2(consFname.c_str(),w2Fname.c_str(),gaussFname.c_str()), x(Eigen::Matrix<double,NJOINTS*2, 1>()),
          cost1(0.0),cost2(0.0)
    {
      std::cout<<"***************************\n Using following files to create costComputeCritic Object:\n "<<consFname<<"\n"<<gaussFname<<"\n"<<w1Fname<<"\n"<<w2Fname<<"\n";
      if (em != NULL){
        em->startManaging(*this);
      }

    }

  protected:
    virtual void operate(){
      jPos=jPosInp.getValue();
      jVel=jVelInp.getValue();

      x<<jPos[1],jVel[1],jPos[3],jVel[3];

      cost1 = fuz1.output(x);
      cost2 = fuz2.output(x);

      cost1OutputVal->setData(&cost1);
      cost2OutputVal->setData(&cost2);


    }

  };
}
#endif /* COSTCOMPUTECRITIC_HPP_ */
