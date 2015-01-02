//#include <system.h>
//#include <QuadTS.h>
//#include"system.h"
//==========================================================
QuadTS::QuadTS()
{

  m_NI=0;
  m_NW=0;
  m_NF=0;

  m_NZ = 0;
  m_NC = 0;

  m_Share = 0;

  m_vW = 0;
  m_vdelW = 0;


  m_mC = 0;
  m_vMu = 0;
  m_vSigma = 0;
  m_vWnet = 0;

  m_Out = 0;
  m_Er = 0;

  m_vden = 0;


  m_vxB = 0;

  m_vMuW = 0;

  m_vDelMu = 0;


  tval=0.0;
  musum =0.0;


#ifdef DEBUG

  warning("QuadTS is created with zero initialization");

#endif


}
//==========================================================
QuadTS::~QuadTS()
{

  clearBuffers();

}
//==========================================================
void QuadTS::clearBuffers()
{

  del(m_vW,m_NC);
  del(m_vdelW,m_NC);
  del(m_mC,m_NC);
  del(m_NZ);

  del(m_vMu);
  del(m_normMemVal);
  del(m_vSigma);
  del(m_vWnet);


  del(m_vden);

  del(m_vMuW);


  del(m_vxB);


  del(m_vDelMu);


}
//==========================================================
void QuadTS::createBuffers()
{

  uint i=0;


  //==============================
  //Varaible Memory allocations

  m_NC = 1;

  //calculating the number of fuzzy zones
  for(i = 0;i < m_NF;i++)
    m_NC *= m_NZ[i];

  //weight matrix allocation
  m_vW = allocate<double>(m_NC, m_NW );

  //change in weight matrix allocation
  m_vdelW = allocate<double>(m_NC, m_NW );

  m_mC = allocate<double>(m_NC,m_NF);
  m_vMu = allocate<double>(m_NC);
  m_normMemVal = allocate<double>(m_NC);
  m_vDelMu = allocate<double>(m_NC*m_NF);

  m_vWnet = allocate<double>(m_NW);

  m_vMuW = allocate<double>(m_NW);


  m_vxB = allocate<double>(m_NW);

  //temp Variables needed in ::ComputeDelX
   tval=0.0;
   musum =0.0;

  //=============================
#ifdef DEBUG

  display("Number of Inputs::",m_NI);
  display("Number of Fuzzy Input::",m_NF);

  display("Number of fuzzy zones of inidividual fuzzy layer states :: ");
  for(i=0;i<m_NF;i++)
    std::cout<<i<<"::\t"<<m_NZ[i]<<std::endl;

  newline();

  display("Number of fuzzy zones ::",m_NC);
  newline();



#endif

  //==================================

  return ;
}//==========================================================
void QuadTS::initFuzzyZone(const double* const& xMax, const double* const& xMin)
{ 

  uint i=0;
  uint j=0;


  m_vSigma = allocate<double>(m_NF);
  m_vden = allocate<double>(m_NF);


  //==============================
  //computing width of fuzzy zones

  double* width= NULL;
  width = allocate<double>(m_NF);

  uint nf = m_NF;
  sub(xMax,xMin, width,nf);
#ifdef DEBUG
  display("Range of Operation :: ");
  display(width,m_NF);
#endif
  //===============================
  //Initializing Fuzzy centers

  
  //Linear System
  if(m_NC == 1)
    {
      add(xMax, xMin, m_mC[0] , m_NF);
      divide(m_mC[0], double(2.0), m_NF);

    }
  else
    {

      //=====================================
      //computing the width of fuzzy zones
      for(i=0;i<m_NF;i++)
	{
	  if(m_NZ[i]!=1)
	    width[i] = width[i]/(m_NZ[i]-1);
	  else
	    width[i] = ( xMax[i] + xMin[i] )/2;
	}
      //=====================================
      uint *nb=0;

      //temporary buffer to store the number of repetition terms
      nb = allocate<uint>(m_NF);
      //======================================
      //number of occurences of each zone of inidividual variables
      for(i=0;i < m_NF;i++)
	{
	  if(i==0)
	    nb[i] = 1;
	  else
	    nb[i] = nb[i-1] * m_NZ[i-1];
	}
      //=======================================

#ifdef DEBUG
      display("Width of fuzzy zones::");
      display(width,m_NF);
      display("zone centers created::");
      display(nb,m_NF);
#endif

      //=====================================
      //Arranging the centers
      for(i=0; i<m_NC; i++)
	{
	  for(j=0; j<m_NF; j++)
	    {
	      if(m_NZ[j] != 1)
		m_mC[i][j]= xMin[j] + ( (i/nb[j]) % m_NZ[j] ) * width[j];
	      else
		m_mC[i][j] = xMin[j] + width[j];

	    }

	}

      delete []nb;
    } 

  //=================================
#ifdef DEBUG
  display("Centers of Fuzzy Zone::");
  for(uint i=0;i<m_NC;i++)
    display(m_mC[i], m_NF);
#endif
  //==============================
  //Initializing the variance

  initFuzzySigma(width);

  //================================

  del(width);
  return ;

}
//==========================================================
QuadTS::QuadTS(const uint& ni, const uint* const& nz, const double* const &xMax, const double* const &xMin, const double &share, const double* const &w)
{

  uint i = 0;
  srand(time(NULL));

  m_NI = ni;
  m_NW = m_NI* (m_NI+1)/2;
  m_NF = ni;

  //==============================
  //Varaible Memory allocations

  m_NZ = allocate<uint>(m_NF);

  //calculating the number of fuzzy zones
  copy(nz,m_NZ,m_NF);


  createBuffers();

  //===============================
  //Initializing the weights
  for(i=0;i<m_NC;i++)
    {
      if(w != 0)
	{
	  copy(w, m_vW[i], m_NW );
	}
      else
	{
	  randomize(m_vW[i], m_NW, double(0.1));
#ifdef DEBUG
	  set(m_vW[i],double(i+1),m_NW );
#endif
	}


#ifdef DEBUG
      display("zone::  ",i);
      display(m_vW[i],m_NW);
#endif


    }


  //==============================
  //computing width of fuzzy zones
  m_Share = share;
  initFuzzyZone(xMax,xMin);

#ifdef DEBUG
  display("Fuzzy zone created");
  display("Range of fuzzy layer states ::");
  display("Maximum::",xMax, m_NF);
  display("Minimum::",xMin,m_NF);
#endif

  //==============================


}
//==========================================================
QuadTS::QuadTS(const char* const pname, const char* const wname, const char* const cname)
{


	//uint  i = 0;

  readParam(pname);

#ifdef DEBUG
  display("Parameter file read");
#endif

  createBuffers();

  readWeight(wname);

#ifdef DEBUG
  display("Weight file read");
#endif

#ifdef DEBUG
  for(i=0;i<m_NC;i++)
    {
      display("weights::  ",i);
      display(m_vW[i],m_NW);
    }
#endif

  readCenter(cname);

#ifdef DEBUG
  for(i=0;i<m_NC;i++)
    {
      display("zone::  ",i);
      display(m_mC[i],m_NF);
    }
#endif



#ifdef DEBUG
  display("QuadTS created from file successfully");
#endif


}

//==========================================================
QuadTS::QuadTS(const uint& ni, const uint& nf, const uint* const& nz, const double* const &xMax, const double* const &xMin, const double &share, const double* const &w)
{

  uint i = 0;
  srand(time(NULL));

  m_NI = ni;
  m_NW = m_NI * (m_NI + 1) /2;
  m_NF = nf;

  //==============================
  //Varaible Memory allocations
  m_NZ = allocate<uint>(m_NF);


  //calculating the number of fuzzy zones
  copy(nz,m_NZ,m_NF);

  createBuffers();


  //===============================
  //Initializing the weights
  for(i=0;i<m_NC;i++)
    {
      if(w != 0)
	{
	  copy(w, m_vW[i], m_NW);
	}
      else
	{
	  randomize(m_vW[i], m_NW, double(0.1));
#ifdef DEBUG
	  set(m_vW[i],double(i+1),m_NW );
#endif
	}


#ifdef DEBUG
      display("zone::  ",i);
      display(m_vW[i],m_NW);
#endif


    }

  //==============================
  //computing width of fuzzy zones
  m_Share = share;
  initFuzzyZone(xMax,xMin);

#ifdef DEBUG
  display("Fuzzy zone created");
  display("Range of fuzzy layer states ::");
  display("Maximum::",xMax, m_NF);
  display("Minimum::",xMin,m_NF);
#endif


  //==============================


}
//==========================================================
void QuadTS::readWeight(const char* const fname)
{
  uint i=0;
  uint j=0;

  std::ifstream fw;
  open2read(fw,fname);



  for(i=0;i<m_NC;i++)
    {
      for(j=0;j<m_NW;j++)
	{
	  fw>>m_vW[i][j];
	}
    }


  fw.close();


}
//==========================================================
void QuadTS::readCenter(const char* const fname)
{
  uint i=0;
  uint j=0;

  std::ifstream fw;
  open2read(fw,fname);



  for(i=0;i<m_NC;i++)
    {
      for(j=0;j<m_NF;j++)
	{
	  fw>>m_mC[i][j];
	}
    }


  fw.close();


}
//==========================================================
void QuadTS::readParam(const char* const fname)
{

  uint i=0;

  std::ifstream fw;
  open2read(fw,fname);

  fw>>m_NI;
  m_NW = m_NI * (m_NI+1)/2;
  fw>>m_NF;


  m_NZ=allocate<uint>(m_NF);
 
  for(i=0;i<m_NF;i++)
    fw>>m_NZ[i];

  m_vSigma = allocate<double>(m_NF);
  m_vden = allocate<double>(m_NF);

  double* sigma = allocate<double>(m_NF);

  for(i=0;i<m_NF;i++)
    fw>>sigma[i];

  setFuzzySigma(sigma);

  fw>>m_Share;

  del(sigma);
  fw.close();

  return;
}

//==========================================================
void QuadTS::storeWeight(const char* const fname)const
{

  uint i=0;

  std::ofstream fw;

  if(fname !=0)
    {
      open2write(fw,fname);
    }
  else
    {
      open2write(fw,"weight.txt");
    }



  for( i = 0;i<m_NC;i++)
    {

      print(fw,m_vW[i],m_NW);
    }

  fw<<std::endl;

  fw.close();

#ifdef DEBUG
  display("\nWeights are stored ");
#endif


}
//==========================================================
void QuadTS::storeCenter(const char* const fname)const
{

  uint i=0;

  std::ofstream fw;

  if(fname !=0)
    {
      open2write(fw,fname);
    }
  else
    {
      open2write(fw,"center.txt");
    }


  for( i = 0;i<m_NC;i++)
    {

      print(fw,m_mC[i],m_NF);
    }

  fw<<std::endl;

  fw.close();

#ifdef DEBUG
  display("\nCenters are stored ");
#endif


}
//==========================================================
void QuadTS::storeParam(const char* const fname)const
{
//  uint i=0;

  std::ofstream fw;

  if(fname !=0)
    {
      open2write(fw,fname);
    }
  else
    {
      open2write(fw,"param.txt");
    }


  fw<<m_NI<<std::endl;
  fw<<m_NF<<std::endl;
  fw<<std::endl;

  print(fw,m_NZ,m_NF);
  fw<<std::endl;

  print(fw,m_vSigma,m_NF);
  fw<<std::endl;

  fw<<m_Share;
  fw<<std::endl;


#ifdef DEBUG
  display("\nParameters are stored");
#endif


  fw.close();
  return;

}
//==========================================================
void QuadTS::getParam(uint& ni, uint& nf, uint* &nz, double& share)const
{

  ni = m_NI;
  nf = m_NF;

  copy(m_NZ,nz,nf);

  share = m_Share;
  return;

}

//==========================================================
void QuadTS::getWeight(double** const &w)const
{

  uint i = 0;


  for(i=0; i<m_NC;i++)
    {
      copy(m_vW[i], w[i], m_NW);
    }

  return;
}

//==========================================================
void QuadTS::getFuzzyCenter(double** const &c)const
{

  uint i = 0;

  for(i=0; i<m_NC;i++)
    {
      copy(m_mC[i], c[i], m_NF);
    }

  return;
}

//==========================================================
void QuadTS::getFuzzySigma(double* const& sigma)const
{

//  uint i=0;

  copy(m_vSigma,sigma,m_NF);

  return;
}
//==========================================================
void QuadTS::computeMembership(const double* const &x)
{


  uint i=0;
  uint j=0;

//  double norm;

  set(m_vMu, double(1.0), m_NC);

  //========================
  for(i = 0;i < m_NC;i++)
    {
      for(j = 0;j<m_NF;j++)
	{
//    	  if((x[j]>ThMax_[j] && abs((ThMax_[j]-m_mC[i][j]))<0.01) || (x[j]<ThMin_[j] && abs((ThMin_[j]-m_mC[i][j]))<0.01)){
//    		  norm=1;
//    		  std::cout<<"entered IF::"<<std::endl<<"center: "<<m_mC[i][j]<<"  Xval: "<<x[j]<<"  ThMx and min:m "<<ThMax_[j]<<"  "<<ThMin_[j]<<"  and MU: "<<norm<<std::endl;
//
//
//    	  }
//    	  else{
    		  norm = exp( -pow((m_mC[i][j]-x[j] ),2)/m_vden[j] );
//    		  std::cout<<"entered IF::"<<std::endl<<"center: "<<m_mC[i][j]<<"  Xval: "<<x[j]<<"  ThMx and min:m "<<ThMax_[j]<<"  "<<ThMin_[j]<<"  and MU: "
//    		      				  <<exp( -pow((m_mC[i][j]-x[j] ),2)/m_vden[j] )<<std::endl;

//    	  }

	  m_vMu[i] = std::min(norm ,m_vMu[i]);
//	  std::cout<<i<<"th mu val:  "<<m_vMu[i]<<std::endl;

	}
    }

  //=========================
#ifdef DEBUG
  display("Membership Values & centers::");
  for(i=0;i<m_NC;i++)
    {
      std::cout<<i<<"\t"<<m_vMu[i]<<"\t";
      display(m_mC[i], m_NF);
    }
#endif

  double musum = sum(m_vMu,m_NC);
  if(musum< 0.001)
    {
      display("given x::",x,m_NF);
      display("Sum of membership", musum );
      exit(EXIT_FAILURE);
//       display(m_vMu,m_NC);
    }


  //==========================
  return ;
}
///==========================================================
void QuadTS::computeDelW(const double& er, const double* const &In)
{

  uint i=0;

  //==========================
#ifdef DEBUG
  std::cout<<"_____________Computing DelW________________\n";
  display("eta * Er::", er);
  display("In::",In,m_NI);
#endif
  //==========================
 

  //Computing delW
  //for i^{th} fuzzy zone delW_i = eta * Er^T * (-1) * mu_i*  xBasis
  for(i = 0;i<m_NC;i++)
    {

      //copying xBasis vector
      copy(m_vxB, m_vdelW[i], m_NW);      


      //delW_i = -1* eta * er * mu_i* basis (-1 is multiplied to include negative
      //sign associated with y in E^2 term)
      multiply(m_vdelW[i], -(er* m_vMu[i]), m_NW);

      //========================
#ifdef DEBUG
      display("Membership value and zone ::  ",m_vMu[i]);
      display(m_mC[i], m_NF);
      display("DelW::");
      display(m_vdelW[i],m_NW);
#endif
      //========================
    }

  //============================
  return ;

}
//==========================================================
void QuadTS::computeBasis(const double* const& x, double* const& xB)
{
  //computes the Basis vector for the fuzzy system
  //The vector is the Quadratic polynomials of the input
  // V = \frac{1}{2} x^T W x

  //==========================

  uint i=0;
  uint j=0;
  uint n=0;

  for(i = 0; i < m_NI; i++)
    {
      for(j=i; j< m_NI; j++)
	{


          if( i == j)
	    xB[n] = 0.5 * pow(x[i],2.0);
	  else
	    xB[n] = x[i] * x[j];

	  n = n+1;
	}
    }
  //==========================

#ifdef DEBUG
  display("Computing the gradient w.r.t W");
  display("Input vector::", x, m_NI);
  display("Gradient Vector::", xB, m_NW);
#endif


  //==========================
  return;
}
//==========================================================
const double& QuadTS::Output(const double* const &x, const double* const &fx)
{

  double musum=0;
  uint i=0;

  zero(m_vWnet, m_NW );

  //=======================
  //Membership Computation
  if(fx !=0)
    {
#ifdef DEBUG
      display("computing membership with fx");
#endif
      computeMembership(fx);
    }
  else
    {
#ifdef DEBUG
      display("computing membership with x ");
#endif
      computeMembership(x);
    }
  musum = sum(m_vMu, m_NC);
  divide(m_vMu, musum, m_NC);
  //==========================
#ifdef DEBUG
  display("musum ::  ",musum);
  display("Normalized membership and Fuzzy zones");
  for(i=0;i<m_NC;i++)
    {
      std::cout<<m_vMu[i]<<"\t";
      display(m_mC[i],m_NF);
    }
  newline();
#endif
  //==========================
  //Net weight Computation

  for(i = 0;i<m_NC;i++)
    {

      multiply(m_vW[i], m_vMu[i], m_vMuW, m_NW);

      add(m_vWnet, m_vMuW, m_NW);
    }
  //==========================
#ifdef DEBUG
  display("Net Weight matrix ::");
  display(m_vWnet,m_NW);
#endif

  //==========================
  //Basis computation

  computeBasis(x,m_vxB);

  //==========================
  //Output Computation	  
  m_Out = dotproduct(m_vWnet,m_vxB, m_NW );
  //==========================
#ifdef DEBUG
  display("Output::", m_Out);
#endif

  //==========================
  return(m_Out);

}
//==========================================================
void QuadTS::update(const double* const &In, const double& yd, const double &eta, const double* const &fx)
{

  uint i=0;

  //  //Compute Output for the given input
  if(fx !=0) // if fuzzy zone is different from input
    {
      m_Out = Output(In,fx);
    }
  else //If input is fuzzified
    {
      m_Out=Output(In);
    }

  //compute the error
  m_Er = yd - m_Out;

#ifdef DEBUG
  display("Error :: ",m_Er);
#endif

  m_Er = m_Er * eta;

  //==============================
#ifdef DEBUG
  display("Error x eta :: ", m_Er);
#endif
  //==============================

  //computing Del W
  computeDelW(m_Er,In);
  //===============================

  //Updating the weight

#ifdef DEBUG
  display("Updated Weight ::");

#endif

  // w(new)= w(old) - delw
  // w(ew) = w(old) - eta* (-1)* Er^T * In
  for(i=0;i<m_NC;i++)
    {
      sub(m_vW[i], m_vdelW[i], m_NW);
#ifdef DEBUG
      display(m_vW[i],m_NW);
      newline();
#endif
    }
  //  return ;
}
//==========================================================
QuadTS::QuadTS(const QuadTS &b)
{


  //==============================

  uint* nz=allocate<uint>(MAX_DATA_INPUT);
  double* sigma=0;

//  uint i=0;

  //==============================

  b.getParam(m_NI, m_NF, nz, m_Share);
  m_NW = m_NI * (m_NI + 1) /2;

  m_NZ = allocate<uint>(m_NF);

  //calculating the number of fuzzy zones
  copy(nz,m_NZ,m_NF);
 
  //==============================
  createBuffers();

  b.getWeight(m_vW);

  b.getFuzzyCenter(m_mC);

  //==============================


  sigma = allocate<double>(m_NF);
  m_vSigma = allocate<double>(m_NF);
  m_vden = allocate<double>(m_NF);

  b.getFuzzySigma(sigma);
  setFuzzySigma(sigma);


  //==============================
 
  del(sigma);
  del(nz);


  //==============================

#ifdef DEBUG
  display("QuadTS Copied for constructor successfully\n");
#endif

  //==============================

  return ;
}

//==========================================================
QuadTS& QuadTS::operator=(const QuadTS& b)
{

  uint* nz=allocate<uint>(MAX_DATA_INPUT);
  double* sigma = 0;
//  uint i=0;

  //==============================
  clearBuffers();


#ifdef DEBUG

  display("Existing network is cleared");
  display(" Empty network is ready for assignment");

#endif

  //==============================
  b.getParam(m_NI, m_NF, nz, m_Share);

  m_NW = m_NI * (m_NI+1)/2;

  m_NZ = allocate<uint>(m_NF);

  //==============================
  //calculating the number of fuzzy zones
  copy(nz,m_NZ,m_NF);
 
  //==============================

  createBuffers();

  b.getWeight(m_vW);
 
  b.getFuzzyCenter(m_mC);


  //==============================

  sigma = allocate<double>(m_NF);
  m_vSigma = allocate<double>(m_NF);
  m_vden = allocate<double>(m_NF);


  b.getFuzzySigma(sigma);
  setFuzzySigma(sigma);


  //==============================
 
  del(sigma);
  del(nz);

  //==============================

#ifdef DEBUG
  display("Assignment Successful\n");
#endif

  //==============================

  return (*this);

}
//==========================================================
void QuadTS::setWeight(double** &w)
{

  for(uint i=0;i<m_NC;i++)
    {
      copy(w[i],m_vW[i],m_NW);
    }

  return;
}
//==========================================================
void QuadTS::setFuzzyCenter(double** &c)
{

  for(uint i=0;i<m_NC;i++)
    {
      copy(c[i],m_mC[i],m_NF);
    }

  return;
}
//==========================================================
void QuadTS::initFuzzySigma(const double* const& width)
{

  uint i=0;


  for(i = 0;i < m_NF;i++)
    {
      m_vSigma[i] = sqrt( -pow(width[i],2)/(2*log(m_Share)) );
    }



  for(i = 0;i < m_NF;i++)
    m_vden[i] = 2 * pow(m_vSigma[i],2);


  //========================

#ifdef DEBUG 
  display("Received Width",width,m_NF);
  display("Standard Deviation of fuzzy centers",m_vSigma,m_NF);
  display("Denominator of fuzzy centers",m_vden,m_NF);
#endif

  //========================

  return;
}
//==========================================================
void QuadTS::setFuzzySigma(const double* const& sigma)
{

  uint i=0;

  copy(sigma,m_vSigma,m_NF);

  for(i = 0;i < m_NF;i++)
    m_vden[i] = 2 * pow(m_vSigma[i],2);



  //========================
#ifdef DEBUG 
  display("Standard Deviation of fuzzy centers");
  display(m_vSigma,m_NF);
#endif
  //========================

  return;
}
//==========================================================
void QuadTS::computeDelX(const double* const& x, double* const& delx, const double* const& fx/*=0*/){

////  double delV1[m_NI];
////  double delV2[m_NI];
////  double P[m_NI * m_NI];
//
  tval= 0.0;
  musum = 0.0;
//  //=======================
//  //Membership Computation
  if(fx !=0)
    {

      computeMembership(fx);
      computeDelMu(fx);

    }
  else
    {

      computeMembership(x);
      computeDelMu(x);
    }
//
  for(uint nc = 0; nc < m_NC; nc++){
    musum = musum + m_vMu[nc];
  }
//
//  //==========================
//
  for(uint index = 0; index < m_NI; index ++){
    delV1[index] = 0.0f;
    delV2[index] = 0.0f;
  }
//
////  double delmuSum[m_NF];
//
  for(uint nf=0; nf < m_NF; nf++){
    delmuSum[nf] = 0.0;
    for(uint nc = 0; nc < m_NC; nc++){
      delmuSum[nf] += m_vDelMu[nc*m_NF + nf];
    }
  }
//
//
  for(uint nc = 0; nc < m_NC; nc++){


    //=============================
    uint count = 0;
    for(uint row = 0; row < m_NI; row++){
      for(uint col = row; col < m_NI; col++){
	  P[row*m_NI + col] =  m_vW[nc][count];
	  P[col*m_NI + row] =  m_vW[nc][count];
	count ++;

      }
    }//for loop - forming individual P matrix


    //=============================
    //computation delV1

    for(uint row = 0; row < m_NI; row++){
      tval = 0;
      for(uint col = 0; col < m_NI; col++){
	tval = tval +  P[row*m_NI + col] * x[col];
      }
      delV1[row] = delV1[row] + m_vMu[nc]/musum * tval;
    }
//
//    //=============================
//    //computation of delV2
//
//
    tval = 0;
    for(uint row = 0; row < m_NI; row++){
      for(uint col = 0; col < m_NI; col++){
	tval = tval + 0.5* x[row] * P[row*m_NI + col] * x[col];
      }
    }
//
//
//

    for(uint nf = 0; nf < m_NF; nf++)
      delV2[nf] = delV2[nf] + tval * (m_vDelMu[nc*m_NF + nf] - m_vMu[nc]/musum * delmuSum[nf]) / musum;

  }
//
//
//  //==========================
////     delV2[1] = 0.0f;
//
//
  for(uint index = 0; index < m_NI; index ++){
    delx[index] = delV1[index] + delV2[index];
  }
//
//
//  //==========================
  return;
}
//==========================================================
void QuadTS::computeDelMu(const double* const& x){

  for(uint nc = 0; nc < m_NC; nc++){
    for(uint nf = 0; nf < m_NF; nf++){
      //       m_vDelMu[nc*m_NF + nf] = m_vMu[nc] * (-1/pow(m_vSigma[nf], 2.0) ) * (x[nf] - m_mC[nc][nf]);
      m_vDelMu[nc*m_NF + nf] = m_vMu[nc] * (-2/m_vden[nf] ) * (x[nf] - m_mC[nc][nf]);
    }
  }


#ifdef DEBUG
  display("DelMu::");
  display(m_vDelMu, m_NC, m_NF);
#endif
  return;
}
//==========================================================
//void QuadTS::compNormalizedMemVal(const double* const& x, double* const& normMemVal){
//	  computeMembership(x);
//	  double musum = sum(m_vMu,m_NC);
//	  for (size_t i=0; i<m_NC; i++){
//		  normMemVal[i]=m_vMu[i]/musum;
//	  }
//  }
 double* QuadTS::compNormalizedMemVal(const double* const& x){
	  computeMembership(x);
	  double musum = sum(m_vMu,m_NC);
//	 std::cout<<"Musum:  "<< musum<<std::endl;
	  for (size_t i=0; i<m_NC; i++){
		  m_normMemVal[i]=m_vMu[i]/musum;
	  }
	  return m_normMemVal;
  }
//==========================================================
//EOF

