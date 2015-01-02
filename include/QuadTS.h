/*
============================================================
created on 04.08.2010 (PKP)

===========================================================
*/


#ifndef _QUADTS_
#define _QUADTS_

#include <constants.h>


#define MAX_DATA_INPUT (20)


//=========================================================
class QuadTS
{

 private:


  //No. of Inputs
  uint m_NI;                          


  //No. of weights (NI * (NI+1)/2)  
  uint m_NW;                          

  //Array storing the number of fuzzy zones of individual inputs
  uint *m_NZ;                         

  //Total number of fuzzy zones in workspace
  uint m_NC;                          

  //No. of fuzzy inputs
  uint m_NF;                          


  //contribution of neighbour guassian functions at the center 
  double m_Share;

  //Weight matrix at each fuzzy zone
  double **m_vW;

  //Centers of fuzzy zone
  double **m_mC;

  //Membership Value
  double* m_vMu;

  //normalized membership value
  double* m_normMemVal;
  //Width of guassian function for each inputs
  double* m_vSigma;
  

  //Net weight matrix (sum of product of guassian function and weight)
  double* m_vWnet;

  //Output
  double m_Out;

  //Error 
  double m_Er;

   //Varinace of individual fuzzy input zone
  double* m_vden;

  //change in  weight vector
  double** m_vdelW;

  //W_{\mu_i} = \sigma_i * W_i
  double* m_vMuW;

  //Input Basis vector to the network
  double* m_vxB;

//  const double* const& getMemValues(){
//	  return m_vMuW;
//  }
  //creates members of appropriate size
  void createBuffers();

  //create NC fuzzy zones by computing the centers
  void initFuzzyZone(const double* const& xMax, const double* const& xMin);


  //computing the change in weight
  void computeDelW(const double& Er, const double* const &In);


  //Clears the dynamically allocated memories
  void clearBuffers();


 //Computing the fuzzy membership of given fuzzy input
  void computeMembership(const double* const &x);



  //computes the Basis vector for the fuzzy system
  //The vector is the Quadratic polynomials of the input
  void computeBasis(const double* const& x, double* const& gradW);

  //===============================
  //critic portion
  double* m_vDelMu;

  void computeDelMu(const double* const& x);


 private: //variables which were defined inside the function::ComputeDelX
  double delV1[4];
    double delV2[4];
    double P[16];
    double delmuSum[4];
    double tval;
    double musum;
    double norm;


 public:




  //=====================================================
  //CONSTRUCTOR AND DESTRUCTOR



  //default constructor
  QuadTS();

  ~QuadTS();


  //constructor for pre-defined initialization with same fuzzy and system input
  QuadTS(const uint& ni, const uint* const& nz, const double* const &xMax, const double* const &xMin, const double& share, const double* const &w=0);


  //Constructor for pre-defined weight initialization with separate fuzzy and system input
  QuadTS(const uint& ni, const uint& nf, const uint* const& nz, const double* const& xMax, const double* const& xMin, const double& share, const double* const& w=0);


  //Fuzzy Initialized from file
  QuadTS(const char* const pname, const char* const wname, const char* const cname);


  
  //copy constructor
  QuadTS(const QuadTS& b);

  //===================================================
  //ACCESSING THE ELEMENTS OF QuadTS-FUZZY MODEL



  //get Parameters of fuzzy model
  void getParam(uint& ni, uint& nf, uint* &nz, double& share)const;

  //get Weight Matrix
  void getWeight(double** const &w)const;


  //get Weight Matrix
  void getFuzzyCenter(double** const &c)const;



  //get the variance of fuzzy zones
  void getFuzzySigma(double* const& sigma)const;

  //store weights in pre-defined file
  void storeWeight(const char* const fname=0)const;

  //store Parameters in give file
  void storeParam(const char* const fname=0)const;


  //store Parameters in give file
  void storeCenter(const char* const fname=0)const;



  //read weights from given file 
  void readWeight(const char* const fname);
 
  //read Parameters from given file
  void readParam(const char* const fname);

  //read Parameters from given file
  void readCenter(const char* const fname);

 


  //==================================================
  //DATA MANIPULATION


  //computing output for given input with separate fuzzy zone activation
  const double& Output(const double* const &x, const double* const& fx=0);


  //updating the weight with separate fuzzy zone activation
  void update(const double* const &In,  const double& yd, const double& eta, const double* const &fx=0);


  //===================================================
  //assignment operator overloading
  QuadTS& operator=(const QuadTS& b);


  //===================================================
  //Setting the TSFuzzy parameters manually

  //initialzing the individual weight matrix
  void setWeight(double** & w);


  //initialzing the individual fuzzy centers
  void setFuzzyCenter(double** & c);


  //initialzing the individual fuzzy width 
  void setFuzzySigma(const double* const& sigma);


  //initialzing the individual fuzzy sigma based on width
  void initFuzzySigma(const double* const& width);



  //===================================================
  //portion of critic

  void computeDelX(const double* const& x, double* const& delx, const double* const& fx=0);


//   void compNormalizedMemVal(const double* const& x, double* const& normMemVal);
  double* compNormalizedMemVal(const double* const& x);


};
//=========================================================
#include <Detail/QuadTS-inl.h>
#endif
//=========================================================

//EOF


