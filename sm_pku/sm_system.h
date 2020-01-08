#ifndef _SM_SYSTEM_H_
#define _SM_SYSTEM_H_
#include<string>
#include<unordered_map>
#include"help.h"
#include"base.h"
#include"orbital.h"
#include"twobodystate.h"
#include"system_ho.h"
#include<eigen3/Eigen/Sparse>//matrix library

using std::string;
using std::unordered_map;


namespace SM
{
  typedef complexd DataType;
  //  typedef double DataType;
  ///orbitals with l,jj,tz
  class Orbital_SM_Tz:public Orbital_ljjtz
  {
  public:
    Orbital_SM_Tz(){};    
  Orbital_SM_Tz(int _l,int _jj,int _tz,DataType _e,int min=-1,int max=-1,int _type=0):Orbital_ljjtz(_l,_jj,_tz),e(_e),low_limit(min),up_limit(max),type(_type){}
    //used when generating partition
    int low_limit;//minimum particles in this orbital
    int up_limit;//maximum  particles in this orbital
    DataType e;
    int type;
  };

  ///orbitals with l,jj,mm,tz
  class Orbital_SM_Tz_M:public Orbital_ljjmmtz
  {
  public:
    Orbital_SM_Tz_M(){}
  Orbital_SM_Tz_M(int _l,int _jj,int _mm,int _tz,DataType _e):Orbital_ljjmmtz(_l,_jj,_mm,_tz),e(_e){}

    DataType e;
  };

  //****************************************************************
  
  ///system hold interaction
  class System_Tz_SM:public System_Tz<Orbital_SM_Tz,DataType>
    {
    public:
      typedef Orbital_SM_Tz OrbitalType;
      typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
      const static int MatEledim;
      //    System_Tz_SM(const string & _OrbFile,const string &_MatFile):OrbFile(_OrbFile),MatFile(_MatFile){}
    System_Tz_SM():lambda1B(-1),lambda2B(-1){}

      void setInputFiles(const string & _OrbFile,const string & _MatFile)
      {
	OrbFile=_OrbFile;
	MatFile=_MatFile;
      }
	
      void setup(const string & _OrbFile,const string & _MatFile)
      {
	using std::cout;
	setInputFiles(_OrbFile,_MatFile);
	setupOrbitals();
      	setupTwoBodyStates();
	setupIndiceIn2BStates();
	setupTwoBodyMat();
      }

      DataType get1B(int bra,int ket)
      {
	//return a==b?Orbitals[a].e:0;
        if(bra==ket) return Orbitals[bra].e;
        if(Orbitals[bra].isInAGroupWith(Orbitals[ket]))
          {
            int group=GroupIndice[bra];
            int a=OrbIndice[bra];
            int b=OrbIndice[ket];
            int dim=Groups[group].size();
            if(a>b) swap(a,b);
            return OneBodyMat[group][mapab_e(a,b,dim)];
          }
        else
          return 0;
      }
      
      //reduced matrix of operator
      DataType get1BOper(int a,int b) const;
      DataType get2BOper(int a,int b,int c,int d,int Jab,int Jcd) const;
      int lambda1B;
      int lambda2B;      
      void setupOper1B(const string & file);
      void setupOper2B(const string & file);
      typedef unordered_map<int,DataType> Table;
      Table Oper1B,Oper2B;

      
      void setupOrbitals();
      void setupTwoBodyMat();

      void printOrbitals() const;
      void printTwoBodyStates() const;
    private:
      string OrbFile;
      string MatFile;
    };



  
  //****************************************************************




  
  ///system hold interaction in M scheme, from J scheme to M scheme
  class System_Tz_SM_M:public System<Orbital_SM_Tz_M,TwoBodyState_TzParM,DataType>
    {
    public:
      typedef System_Tz_SM VSystem;//interaction system
      typedef Orbital_SM_Tz_M OrbitalType;
      typedef TwoBodyState_TzParM TwoBodyStateType;
      typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
      const static int MatEledim;

      void setup(VSystem *Ptr)
      {
	pSystem=Ptr;
	setupOrbitals();
      	setupTwoBodyStates();
	//setupTwoBodyMat();
      }
      void setupOrbitals();
      void setupTwoBodyStates();
      void setupTwoBodyMat();

      void setupIndiceIn2BStates();
      vector<int> IndiceIn2BStates;

      void get2BmatAt(int a,int b,int c,int d,TwoBodyMatEle &MatEle) const;
      void get2BmatOnNeed(int a,int b,int c,int d,TwoBodyMatEle &MatEle) const;


      DataType get1B(int a,int b) const
      {
	if(Orbitals[a].mm==Orbitals[b].mm)
	  return pSystem->get1B(Orders[a],Orders[b]);
	else
	  return 0;
      }
      
      DataType get2B(int a,int b,int c,int d) const
      {
	TwoBodyMatEle MatEle;
	//get2BmatAt(a,b,c,d,MatEle);
	get2BmatOnNeed(a,b,c,d,MatEle);
	if(MatEle.empty()) return 0;
	return MatEle[0];
      }

      //<a|\lambda 0|b>
      DataType get1BOper(int a,int b) const;
      //<ab |\lambda 0|cd >
      DataType get2BOper(int a,int b,int c,int d) const;

      
      double jplus(int a,int b) const;///< J+ matrix element between two orbital a and b
      double jminus(int a,int b) const;///< J+ matrix element between two orbital a and b
      
      void printOrbitals() const;
      void printTwoBodyStates() const;

      VSystem * pSystem;
      vector<int> Orders;//split from which orbital in pSystem->Orbitals, the main shell of each orbital

      int maxProtons;
      int maxNeutrons;

    private:
      int mapab(int a,int b) const//a<b,a map from a<b matrix to one dimension vector
      {
	int totalOrbitals=Orbitals.size();
	return a*totalOrbitals-(a*(a+3))/2+b-1;
      }
    };




  
  //****************************************************************

  const size_t numBits=2048;//number of bits for config., must be greater than totalOrbitals
  ///many body basis, configurations
  class Configs
  {
  public:
    typedef vector<int> Partition;//a partition is a way to distribute particles to orbitals.
    typedef bitset<numBits> Config;
  Configs(const System_Tz_SM_M *Ptr):pSystem(Ptr),maxProtons(pSystem->maxProtons),\
      maxNeutrons(pSystem->maxNeutrons),Orbitals(Ptr->Orbitals)
    {
      //check if numBit is enough
      if(pSystem->Orbitals.size() > numBits)
	{
	  std::cerr<<"numBits("<<numBits<<") "<< "must not be less than totalOrbitals("<<pSystem->Orbitals.size()<<"), please set numBits larger!!!"<<std::endl;
	  exit(0);
	}
    }

    void setvNZ(int _vN,int _vZ)
    {
      //check if A,Z is properly set.
      if(_vZ > maxProtons || _vN > maxNeutrons)
	{
	  std::cerr<<"number of valence protons or neutrons exceed number of orbitals(maxProtons:"<<maxProtons<<",maxNeutrons:"<<maxNeutrons<<")!!!\n";
	  exit(0);
	}
      vZ=_vZ;vN=_vN;      
    }
    
    void addConfigs(int vN,int vZ,int MM,int Par,const vector<int>& restriction,double Emax);///< add configs which has total angular momentum projection MM/2 and parity Par. Par=0 for +, 1 for -.
    void genPartitions(int Par,const vector<int> &restriction,double Emax);//generate partitions with valence neutron(proton) vN(vZ)
    vector<Partition> partitions;//partitions
    void addConfigs(int MM);//generate configs after generating partition
    void clear(){partitions.clear();configs.clear();}
    
    void setup_simple(int MM,int Par);//simple version, don't generate partitions first
    
    vector< Config > configs;///< different configurations; A basis contain a sequence of bits, 0-maxProtons-1 for protons, this is in accordance with the orders of orbitals in pSystem.

    DataType getME(int i,int j) const;///return interaction matrix element between config i and j.
    double getJ2(int i,int j) const;//return J^2 matrix element between config. i and j.

    //operator matrix element between config. i and j
    DataType getOper(int i,int j) const;
    
    const System_Tz_SM_M * pSystem;
    const vector<System_Tz_SM_M::OrbitalType> & Orbitals;
    int maxNeutrons,maxProtons;

    int vN,vZ;///< number of valence nuclei.
    int dim;
  private:
    int getMM(const Config & config) const;///< cal the total angular momentum projection of a config.
    int getPar(const Config & config) const;///< cal the parity of a config.
    DataType get1B(int a,int b) const
    {
      //return a==b?Orbitals[a].e:0;
      return pSystem->get1B(a,b);
    }
    DataType get2B(int a,int b,int c,int d) const
    {
      return pSystem->get2B(a,b,c,d);
    }
    
    //recursion algorithm to part 'num' particles into 'numorbs' parts.
    void part(int num,int numorbs,vector<Partition> & Partitions,int *low_limits,int *up_limits);
  };

  //****************************************************************



}
#endif
