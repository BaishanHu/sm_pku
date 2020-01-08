#ifndef _SM_SOLVER_H_
#define _SM_SOLVER_H_
#include<bitset>
#include<eigen3/Eigen/Eigen>//matrix library
#include<eigen3/Eigen/Sparse>//matrix library

#include<boost/mpi.hpp>//boost mpi library

#include"sm_system.h"
#include"lanczos_mpi.h"

using std::size_t;
using std::bitset;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::SelfAdjointEigenSolver;
using Eigen::ComplexEigenSolver;

using Eigen::SparseMatrix;
using Eigen::Triplet;

using namespace boost::mpi;

namespace SM
{
  class SMSolver
  {
  public:
    typedef Matrix<DataType,Dynamic,Dynamic> Mat;
    typedef Matrix<DataType,Dynamic,1> Vec;
    typedef SparseMatrix<DataType> SpMat;//sparse matrix

  SMSolver(const communicator & _world,const System_Tz_SM_M * Ptr):world(_world),configs(Ptr),numStates(0){}

    void addConfigs(int vN,int vZ,int MM,int Par,const vector<int> &restriction,double Emax);
    void clearConfigs(){configs.clear();}
    void diag(bool calOper=false);//diagonalize directly
    void diag_lanczos(int num,bool savevecs=false,bool calOper=false,int step=1);//lanczos algorithm
    void calOccuNums(int i);//cal. the occupied num of every orbital for state i

    void printOccuNums(int num=-1);
    void printStates(int num, ostream& fout);

    int numStates;
    Vec eigVals;
    Mat eigVecs;
    Vec J2s;
    vector<DataType> OccuNums;//only meaningful on process 0


    void printOperInfo(int num, ostream& fout);
    Mat Oper;//only meaningful on process 0
    int ConfigMM;
    
    Configs configs;


    //some info after lanczos
    int totaldim;
    int numOfparts;
    int partdim;
    communicator diagcom;
    
  private:
    const communicator & world;
  };

  //****************************************************************

  
}

#endif
