
/********************************************
 * This class outputs B Matrix;
 * Needed for myController to generate critic output.
 * ******************************************
 */
#ifndef BMATRIXGENERATOR_H_
#define BMATRIXGENERATOR_H_

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
class BmatrixGenerator : public systems::System
{

  BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
public:

  typedef math::Matrix<DOF,DOF,double> Mat7D_type;

  typedef math::Matrix<1,DOF*DOF,double> BmatAsVec_type;
  typedef math::Matrix<1,DOF, double> Bu_type;
public:
  Input<Mat7D_type> inputMassMat;
  Input<jt_type> controlInp;

//  Output<Mat7D_type> BmatOutput7J;
  Output<BmatAsVec_type> BmatAsVecOutput;
  Output<Bu_type> BuOutput;

private:
  Mat7D_type mMatInv, mMat;

  BmatAsVec_type BmatAsVec;
  jt_type uInp;
  math::Matrix<DOF,1,double> u;
  Bu_type BuOut;
protected:

//  typename Output<Mat7D_type>::Value* BmatOutput7JVal;
  typename Output<BmatAsVec_type>::Value* BmatAsVecOutputVal;
  typename Output<Bu_type>::Value* BuOutputVal;
public:
  BmatrixGenerator(systems::ExecutionManager* em, const std::string& sysName = "massMatrixHolder"):
    System(sysName), inputMassMat(this), controlInp(this), /* BmatOutput7J(this, &BmatOutput7JVal),*/ BmatAsVecOutput(this, &BmatAsVecOutputVal),
    BuOutput(this, &BuOutputVal), mMatInv(Mat7D_type()),mMat(Mat7D_type()), BmatAsVec(BmatAsVec_type()), uInp(0.0), u(math::Matrix<DOF,1,double>()),
    BuOut(Bu_type())
  {


    if (em != NULL){
      em->startManaging(*this);
    }
  }

  virtual ~BmatrixGenerator(){this->mandatoryCleanUp();}

protected:
  virtual void operate(){

    mMat = inputMassMat.getValue();
    uInp=controlInp.getValue();
    for(size_t i=0; i<DOF; ++i){
      u(i)=uInp[i];
    }
    mMatInv = mMat.inverse();


    BuOut = mMatInv*u;

    for (size_t i=0; i<mMatInv.SIZE; ++i){
      BmatAsVec(i)= *(mMatInv.data()+i);
    }


//    Bmat2JAsVec<<BMat2J(0,0), BMat2J(1,0),BMat2J(0,1),BMat2J(1,1);

//    BmatOutput7JVal->setData(&BMat7J);
    BmatAsVecOutputVal->setData(&BmatAsVec);
    BuOutputVal->setData(&BuOut);
}

private:
    DISALLOW_COPY_AND_ASSIGN(BmatrixGenerator);
};

}


#endif
