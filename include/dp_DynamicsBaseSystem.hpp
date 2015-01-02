/*
 * dp_DynamicsBaseSystem.hpp
 *
 *  Last update: September 16, 2011
 *      Author: Diego Pardo
 *
 *  Long Term Goal : Inverse Dynamics Control
 *
 *  This file computes de JointSpace Inertial Matrix -> M-Matrix
 */


#ifndef IRI_INERTIA_MATRIX_H_
#define IRI_INERTIA_MATRIX_H_


#include <iostream>

#include <libconfig.h>
#include <libconfig.h++>
#include <barrett/cdlbt/gsl.h>

#include <barrett/math/kinematics.h>
#include <barrett/cdlbt/kinematics.h>
#include <barrett/cdlbt/dynamics.h>

//#include <barrett/standard_main_function.h>

#include <barrett/products/product_manager.h>

#include <barrett/detail/ca_macro.h>
#include <barrett/detail/stl_utils.h>

#include <barrett/systems/kinematics_base.h>
#include <barrett/systems/abstract/system.h>
#include <barrett/systems/abstract/single_io.h>
#include <barrett/systems.h>  // includes all non-abstract Systems



//using namespace barrett;

namespace barrett{

template <size_t DOF>
class DynamicsBaseSystem : public systems::System, public systems::SingleOutput<typename math::Matrix<DOF,DOF,double>::type> {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

public:	Input<jp_type> jpInput;
public:	Input<jv_type> jvInput;
  
	
public:
DynamicsBaseSystem(systems::ExecutionManager* em, const libconfig::Setting& misetting, jv_type jv_current, jp_type jp_current,	const std::string& sysName = "DynamicsBaseSystem"):
	systems::System(sysName), //systems::KinematicsInput<DOF>(this),
	systems::SingleOutput<math::Matrix<DOF,DOF,double> >(this),
	jpInput(this),jvInput(this),
	my_kin(misetting["kinematics"])
	{

	//Evaluar la kinematica antes de crear la dinamica : Digo Yo
	//bt_kinematics_eval(my_kin.impl,this->jpInput.getValue().asGslType(),this->jvInput.getValue().asGslType());
	bt_kinematics_eval(my_kin.impl,jp_current.asGslType(),jv_current.asGslType());

	bt_dynamics_create(&my_dynamics,misetting["dynamics"].getCSetting(), DOF);//, my_kin.impl);

	my_dynamics->jsim = gsl_matrix_alloc(my_dynamics->dof,my_dynamics->dof);

	std::cout << "Memoria para la matriz localizada" << std::endl;

	if (em != NULL)
		em->startManaging(*this);
	}
	//pm.getConfig().lookup(pm.getWamDefaultConfigPath())
virtual ~DynamicsBaseSystem(){ this->mandatoryCleanUp(); }


public:

struct bt_dynamics *my_dynamics;
math::Matrix<DOF,DOF,double> getJMatrix(){ return mimat;}
	
	
	
	
protected:
    
	math::Matrix<DOF,DOF,double> mimat;
	//bool flag;

private:
	math::Kinematics<DOF> my_kin;
	

protected:

	virtual void operate(){


		//if(inputsValid())
	    //{    
	  
		bt_kinematics_eval(my_kin.impl,this->jpInput.getValue().asGslType(),this->jvInput.getValue().asGslType());
		
		bt_dynamics_eval_jsim(my_dynamics,my_kin.impl);
		
		gsl_matrix_memcpy(mimat.asGslType(),my_dynamics->jsim);

		this->outputValue->setData(&mimat);
	    //}
	  
	    }
private:
	DISALLOW_COPY_AND_ASSIGN(DynamicsBaseSystem);
	

	    
};


}
#endif /*IRI_INERTIA_MATRIX_H_*/
