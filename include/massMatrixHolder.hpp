
/********************************************
 * This class outputs Mass Matrix using bt library;
 *
 * ******************************************
 */
#ifndef MASSMATRIXHOLDER_H_
#define MASSMATRIXHOLDER_H_

#include <iostream>
#include <libconfig.h++>
#include <barrett/math/matrix.h>
#include <barrett/math/kinematics.h>
#include <barrett/cdlbt/kinematics.h>
#include <barrett/cdlbt/dynamics.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/ca_macro.h>
#include <barrett/systems.h>
#include <barrett/cdlbt/gsl.h>

//#include <barrett/standard_main_function.h>

#include <barrett/products/product_manager.h>

#include <barrett/detail/ca_macro.h>
#include <barrett/detail/stl_utils.h>

#include <barrett/systems/kinematics_base.h>
#include <barrett/systems/abstract/system.h>

using namespace barrett;
using namespace systems;


namespace Sam{

template <size_t DOF, size_t NJOINTS>
class massMatrixHolder : public systems::System
{

  BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
public:

  typedef math::Matrix<NJOINTS,NJOINTS,double> mMat2D_type;
  typedef math::Matrix<DOF,DOF,double> mMat7D_type;

public:
  Input<jp_type> jpInput;
  Input<jv_type> jvInput;

  Output<mMat2D_type> mMat2DOutput, mInv2JOutput;
  Output<mMat7D_type> mMat7DOutput;


public:
  struct bt_dynamics *my_dynamics;
protected:
  mMat7D_type mMat7D,mInv7J;
  mMat2D_type mMat2D,mInv2J;
  math::Kinematics<DOF> my_kin;


protected:
  typename Output<mMat7D_type>::Value* OutputVal_mMat7D;
  typename Output<mMat2D_type>::Value* OutputVal_mMat2D;
  typename Output<mMat2D_type>::Value* mInv2JOutputVal;

public:
  massMatrixHolder(systems::ExecutionManager* em,const libconfig::Setting& misetting,  jp_type jp_current, jv_type jv_current, const std::string& sysName = "massMatrixHolder"):
    System(sysName), jpInput(this),jvInput(this),mMat2DOutput(this, &OutputVal_mMat2D),mInv2JOutput(this,&mInv2JOutputVal),mMat7DOutput(this, &OutputVal_mMat7D),
    mMat7D(mMat7D_type()),mInv7J(mMat7D_type()), mMat2D(mMat2D_type()),mInv2J(mMat2D_type()), my_kin(misetting["kinematics"])
  {
    bt_kinematics_eval(my_kin.impl,jp_current.asGslType(),jv_current.asGslType());

    bt_dynamics_create(&my_dynamics,misetting["dynamics"].getCSetting(), DOF);//, my_kin.impl);

    my_dynamics->jsim = gsl_matrix_alloc(my_dynamics->dof,my_dynamics->dof);

    if (em != NULL){
      em->startManaging(*this);
    }
  }

  virtual ~massMatrixHolder(){this->mandatoryCleanUp();}

protected:
  virtual void operate(){

    bt_kinematics_eval(my_kin.impl,this->jpInput.getValue().asGslType(),this->jvInput.getValue().asGslType());
    bt_dynamics_eval_jsim(my_dynamics,my_kin.impl);
    gsl_matrix_memcpy(mMat7D.asGslType(),my_dynamics->jsim);


    mMat2D(0,0)=mMat7D(1,1);
    mMat2D(0,1)=mMat7D(1,3);
    mMat2D(1,0)=mMat7D(3,1);
    mMat2D(1,1)=mMat7D(3,3);

    mInv7J = mMat7D.inverse();
    mInv2J(0,0)=mInv7J(1,1);
    mInv2J(0,1)=mInv7J(1,3);
    mInv2J(1,0)=mInv7J(3,1);
    mInv2J(1,1)=mInv7J(3,3);


    OutputVal_mMat7D->setData(&mMat7D);
    OutputVal_mMat2D->setData(&mMat2D);
    mInv2JOutputVal->setData(&mInv2J);


}

private:
    DISALLOW_COPY_AND_ASSIGN(massMatrixHolder);


  };


}






#endif
