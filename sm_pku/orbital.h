#ifndef _ORBITAL_H_
#define _ORBITAL_H_

#include<limits>
#include<cmath>
#include<complex>

#include"base.h"
#include"help.h"

using std::abs;
using std::numeric_limits;
using std::complex;
typedef complex<double> complexd;

///orbital with quantum number (l,j,tz)
class Orbital_ljjtz:public Orbital
{
 public:
  Orbital_ljjtz(){}
  Orbital_ljjtz(int _l,int _jj,int _tz):l(_l),jj(_jj),tz(_tz){}
  int l,jj,tz;///< tz=-1 for proton,1 for neutron
  bool isInAGroupWith(const Orbital&) const;
};

//sepration line
//*******************************************************
///orbital with quantum number (par,m,tz)
class Orbital_parmmtz:public Orbital
{
 public:
  Orbital_parmmtz(){}
 Orbital_parmmtz(int _par,int _mm,int _tz):par(_par),mm(_mm),tz(_tz){}
  int par,mm,tz;///< tz=-1 for proton,1 for neutron, par=0 plus,1 negative
  bool isInAGroupWith(const Orbital&) const;///< is used in system_m to group orbitals according to m,par,tz .
};

///orbital with quantum number (l,j,m,tz)
class Orbital_ljjmmtz:public Orbital_parmmtz
{
 public:
  Orbital_ljjmmtz(){}
 Orbital_ljjmmtz(int _l,int _jj,int _mm,int _tz):l(_l),jj(_jj),Orbital_parmmtz(_l%2,_mm,_tz){}
  int l,jj;///< tz=-1 for proton,1 for neutron
};

//sepration line
//*******************************************************


///harmonic oscillator orbital, including proton and neutron
class HO_Orbital_Tz:public Orbital_ljjtz
{
 public:
  HO_Orbital_Tz(){}
 HO_Orbital_Tz(int _n,int _l,int _jj,int _tz,double _e):Orbital_ljjtz(_l,_jj,_tz),n(_n),e(_e){}
  int n;///< HO. quantum number n, single particle energy is (2n+l+3/2)hbaromega
  double e;///< single particle energy
};

//*******************************************************
///orbital with quantum number (l,j)
class Orbital_ljj:public Orbital
{
 public:
  Orbital_ljj(){}
 Orbital_ljj(int _l,int _jj):l(_l),jj(_jj){}
  int l,jj;
  bool isInAGroupWith(const Orbital&) const;
};
//**********************************************************
///harmonic oscillator orbital, not distinguish proton and neutron
class HO_Orbital_T:public Orbital_ljj
{
 public:
  HO_Orbital_T(){}
 HO_Orbital_T(int _n,int _l,int _jj,double _e):Orbital_ljj(_l,_jj),n(_n),e(_e){}
  int n;///< HO. quantum number n, single particle energy is (2n+l+3/2)hbaromega
  double e;///< single particle energy
};

//**********************************************************

///Hartree-Fock orbital.
template<class DataType>
class SHF_Orbital:public Orbital_ljjtz
{
 public:
  SHF_Orbital(){}
 SHF_Orbital(int _l,int _jj,int _tz,DataType _e):Orbital_ljjtz(_l,_jj,_tz),e(_e){}

  void set(int _GroupIndex,int _l,int _jj,int _tz,DataType _e)
 {
   GroupIndex=_GroupIndex;
   l=_l;
   jj=_jj;
   tz=_tz;
   e=_e;
 }
  
  ///wave function under certain bases(eg. under harmonic bases).
  ///GroupIndex: index of the group in the groups of bases 
  ///coeff: coefficients under bases in that group
  int GroupIndex;
  vector<DataType> Coeff;
  bool state;///< 0 for hole, 1 for particle
  double OccuP;///< occupied possibility

  bool operator < (const SHF_Orbital<DataType> & Orb) const
  {
    //    return e<Orb.e;
    //in case that energy is too close
    if( abs(e-Orb.e) > numeric_limits<double>::epsilon() )
      {
    	return (e < Orb.e);
      }
    else
      {
    	return false;
      }
  }
  DataType e;///< single particle energy.
};




///WS orbital, including proton and neutron
class WS_Orbital_Tz:public Orbital_ljjtz
{
 public:
  WS_Orbital_Tz(){}
  WS_Orbital_Tz(int _l,int _jj,int _tz,complexd _e):Orbital_ljjtz(_l,_jj,_tz),e(_e){}
  complexd e;///< single particle energy
};


///orbital with the value of momentum as good quantum number
class K_Orbital_Tz:public Orbital_ljjtz
{
 public:
  K_Orbital_Tz(){}
  K_Orbital_Tz(int _l,int _jj,int _tz,double _e,double _k,double _wk):Orbital_ljjtz(_l,_jj,_tz),e(_e),k(_k),wk(_wk){}
  double e;
  double k;
  double wk;
};

#endif
