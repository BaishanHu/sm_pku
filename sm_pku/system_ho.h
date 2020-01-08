#ifndef _SYSTEM_HO_H_
#define _SYSTEM_HO_H_
#include<algorithm>

#include"help.h"
#include"base.h"
#include"orbital.h"
#include"twobodystate.h"


using std::swap;
using std::min;
using std::max;

///system, twobody state is in Tz representation
/*
  OrbitalType must be or derived from Orbital_ljjtz
 */
template<class OrbitalType,class DataType>
class System_Tz:public System<OrbitalType,TwoBodyState_TzParJ,DataType>
{
 public:
  typedef TwoBodyState_TzParJ TwoBodyStateType;
  typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
  
  void setupTwoBodyStates();

  /// setup three dimension array IndiceIn2BStates.
  void setupIndiceIn2BStates();
  vector< vector< vector<int> > >IndiceIn2BStates;///<IndiceIn2BStates(a,b,J-Jmin) return the indice of |ab,J) in TwoBodyStates. Jmin is abs(ja-jb),noting a<=b


  /// return <ab,J|2B|cd,J>. Note that the MatEle may be empty if you pass a empty TwoBodyMateEle.
  void get2BmatAt(int a,int b,int c,int d,int J,TwoBodyMatEle & MatEle) const;
 
  void set2BmatAt(int a,int b,int c,int d,int J,const TwoBodyMatEle & MatEle);///< set matrix element <ab,J|2B|cd,J>, called by setupTwoBodyMat to write value into the 2Bmat;

};
template<class OrbitalType,class DataType>
void System_Tz<OrbitalType,DataType>::setupTwoBodyStates()
{
  this->TwoBodyStates.clear();
  for(int i=0;i<this->Orbitals.size();i++)
    {
      for(int j=i;j<this->Orbitals.size();j++)
	{
	  int Jmin=abs(this->Orbitals[i].jj - this->Orbitals[j].jj)/2;
	  int Jmax=(this->Orbitals[i].jj + this->Orbitals[j].jj)/2;
	  for(int J=Jmin;J<=Jmax;J++)
	    {
	      int Tz=(this->Orbitals[i].tz + this->Orbitals[j].tz)/2;
	      int Par=(this->Orbitals[i].l + this->Orbitals[j].l)%2;
	      //if i==j, J must be even. |ab,J) = (-1)^{ja+jb-J+1} |ba,J)
	      if( (i==j) && J%2) continue;
	      this->TwoBodyStates.push_back( TwoBodyStateType(i,j,Tz,Par,J) );
	    }
	}
    }
  this->setupTwoBodyChannels();//noting that the twobody states are sorted after this calling.
}

template<class OrbitalType,class DataType>
void System_Tz<OrbitalType,DataType>::setupIndiceIn2BStates()
{
  int totalOrbitals=this->Orbitals.size();
  IndiceIn2BStates.resize(totalOrbitals);
  for(int i=0;i<totalOrbitals;i++) IndiceIn2BStates[i].resize(totalOrbitals);
  for(int a=0;a<totalOrbitals;a++)
    for(int b=a;b<totalOrbitals;b++)
      {
	int Jmin=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
	int Jmax=(this->Orbitals[a].jj+this->Orbitals[b].jj)/2;
	IndiceIn2BStates[a][b].resize(Jmax-Jmin+1);
	for(int J=Jmin;J<=Jmax;J++)
	  IndiceIn2BStates[a][b][J-Jmin]=-1;//initialize as -1
      }
  for(int i=0;i<this->TwoBodyStates.size();i++)
    {
      int a=this->TwoBodyStates[i].a;
      int b=this->TwoBodyStates[i].b;
      int J=this->TwoBodyStates[i].J;
      int Jmin=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
      IndiceIn2BStates[a][b][J-Jmin]=i;
    }
}

template<class OrbitalType,class DataType>
  void System_Tz<OrbitalType,DataType>::set2BmatAt(int a,int b,int c,int d,int J,const TwoBodyMatEle & MatEle)
{
  //make sure a<=b,c<=d,bra<=ket
  int phase=1;
  if(a>b)
    {
      swap(a,b);
      phase*=( (this->Orbitals[a].jj+this->Orbitals[b].jj)/2 -J + 1 )%2?-1:1;
    }
  if(c>d)
    {
      swap(c,d);
      phase*=( (this->Orbitals[c].jj+this->Orbitals[d].jj)/2 -J + 1 )%2?-1:1;
    }
  int Jminab=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
  int Jmincd=abs(this->Orbitals[c].jj-this->Orbitals[d].jj)/2;
  int ab=IndiceIn2BStates[a][b][J-Jminab];
  int cd=IndiceIn2BStates[c][d][J-Jmincd];
  int channelab=this->ChannelIndice[ab];
  int channelcd=this->ChannelIndice[cd];
  //channelab should = channelcd,we should add this error check when debuging
  int bra=this->TwoBodyStateIndice[ab];
  int ket=this->TwoBodyStateIndice[cd];
  //make sure bra<=ket
  if(bra>ket)
    swap(bra,ket);
  TwoBodyMatEle MatEleTemp=MatEle;
  MatEleTemp*=phase;
  this->set2Bmat(channelab,bra,ket,MatEleTemp);
}

template<class OrbitalType,class DataType>
  void System_Tz<OrbitalType,DataType>::get2BmatAt(int a,int b,int c,int d,int J,TwoBodyMatEle & MatEle) const
{
  MatEle.setZero();
  int phase=1;
  //make sure a<=b,c<=d,(ab)<=(cd)
  if(a>b)
    {
      swap(a,b);
      phase*=( (this->Orbitals[a].jj+this->Orbitals[b].jj)/2-J + 1 )%2?-1:1;
    }
  if(c>d)
    {
      swap(c,d);
      phase*=( (this->Orbitals[c].jj+this->Orbitals[d].jj)/2-J + 1 )%2?-1:1;
    }

  int Jminab=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
  int Jmincd=abs(this->Orbitals[c].jj-this->Orbitals[d].jj)/2;
  int Jmaxab=(this->Orbitals[a].jj+this->Orbitals[b].jj)/2;
  int Jmaxcd=(this->Orbitals[c].jj+this->Orbitals[d].jj)/2;

  if( J<max(Jminab,Jmincd) ) return;
  if( J>min(Jmaxab,Jmaxcd) ) return;
  int ab=IndiceIn2BStates[a][b][J-Jminab];
  int cd=IndiceIn2BStates[c][d][J-Jmincd];
  if(ab==-1||cd==-1) return;
  int channelab=this->ChannelIndice[ab];
  int channelcd=this->ChannelIndice[cd];
  //channelab should = channelcd
  if(channelab!=channelcd) return;
  int bra=this->TwoBodyStateIndice[ab];
  int ket=this->TwoBodyStateIndice[cd];
  if(bra>ket)
    swap(bra,ket);
  this->get2Bmat(channelab,bra,ket,MatEle);
  MatEle*=phase;
}





//sepration line
//****************************************************************************************


///System under HO. basis, twobody states is in Tz representation
/*
  dim of MatEle=4,the 4 double are:
      0  interaction
      1  (p1 \cdot p2)/(m hbar omega) + m omega r1 \cdot r2 / hbar
      2  (r1-r2)^2 m omega/(2 hbar)
      3  (p1 \cdot p2)/(m hbar omega)
*/
class System_Tz_HO:public System_Tz<HO_Orbital_Tz,double>
{
 public:
  typedef HO_Orbital_Tz OrbitalType;
  typedef double DataType;
  typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
  const static int MatEledim;//< dimension of TwoBodyMatEle.
  System_Tz_HO(){setQuad(100);};
 System_Tz_HO(const string & _OrbFile,const string & _MatFile):OrbFile(_OrbFile),MatFile(_MatFile){setQuad(100);}

  void setInputFiles(const string & _OrbFile,const string & _MatFile)
  {
    OrbFile=_OrbFile;
    MatFile=_MatFile;
  }
  /// setup orbitals,2B states and 2B mat.
  void setup()
  {
    setupOrbitals();
    setupTwoBodyStates();
    setupTwoBodyMat();
  }
  void setupOrbitals();
  void setupTwoBodyMat();


  ///interface for HF cal.
  double get1B(int A,int bra,int ket) const;
  ///interface for HF cal.
  double get2B(int A,int a,int b,int c,int d,int J) const;

  ///interface for nocore shell model
  double get1B_p2(int bra,int ket) const;
  double get1B_p2_plus_r2(int bra,int ket) const
  {
    return bra==ket?Orbitals[bra].e:0;
  }

  
  ///nuclear two body interaction
  double getV(int a,int b,int c,int d,int J) const;

  ///return <bra|k^2|ket>, dimensionless: k2 b^2
  double k2(int bra,int ket) const;
  ///return <bra|r^2|ket>, dimensionless: r^2/b^2
  double r2(int bra,int ket) const;
  
  
  ///set quadrature points and weights for integration. An option for rL, kL cal.
  void setQuad(const int QuadPoints)
  {
    x.resize(QuadPoints);
    w.resize(QuadPoints);
    gauleg(0,1,x,w);
    //use t=1/(1+x) map [0,infinity] to [0,1] 
    for(int i=0;i<QuadPoints;i++)
      {
    	double xx=x[i];
    	x[i]=(1-xx)/xx;
    	double ww=w[i];
    	w[i]=ww/pow(xx,2);
      }
  }
  
  //other maybe useful operators
  ///return <bra|r^L|ket>, note that this is just integration of radial wavefunction.
  double rL(int bra,int ket,int L) const;
  ///return <bra|k^L|ket>, the result should be multiplied by a phase (i)^lbra * (-i)^lket
  ///ref. "ho_k" in file "system.h". We multiply the phase in the final expresion to avoid handling with complex
  ///called by "k1k2"
  double kL(int bra,int ket,int L) const;
  
  ///return \f$ <ab J|k1 \cdot k2|cd J> \f$
  ///final result is dimesionless: k1 \cdot k2 b^2 = p1 \cdot p2/(m hbar omega)
  double k1k2(int a,int b,int c,int d,int J) const;
  ///return \f$ <ab J|r1 \cdot r2|cd J>\f$
  /// final result is dimesionless: r1 \cdot r2 m omega/hbar= r1 \cdot r2 m omega/hbar
  double r1r2(int a,int b,int c,int d,int J) const;
  ///return <ab J|(r1-r2)^2|cd J>
  ///final result is dimesionless: (r1-r2)^2 / b^2 = (r1- r2)^2 m omega/hbar
  double relr2(int a,int b,int c,int d,int J) const;
  ///return <ab J|(k1-k2)^2|cd J>
  ///final result is dimesionless: (k1-k2)^2 b^2 = (p1 -p2)^2/(m hbar omega)
  double relk2(int a,int b,int c,int d,int J) const;


  void printOrbitals() const;

  
  double hbar_omega;
  double lengthb;  
 private:
  string OrbFile;
  string MatFile;
  ///for integration
  vector<double> x;
  vector<double> w;
};




//****************************************************************************************

///System, twobody state is in T representation
/*
  OrbitalType must be or derived form Orbital_ljj
*/
template<class OrbitalType,class DataType>
class System_T:public System<OrbitalType,TwoBodyState_TParJ,DataType>
{
 public:
  typedef TwoBodyState_TParJ TwoBodyStateType;
  typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
  
  void setupTwoBodyStates();

  /// setup four dimension array IndiceIn2BStates.
  void setupIndiceIn2BStates();
  vector< vector< vector< vector<int> > > >IndiceIn2BStates;///<IndiceIn2BStates(a,b,T,J-Jmin) return the indice of |ab,JT) in TwoBodyStates. Jmin is abs(ja-jb),noting a<=b


  /// return <ab,JT|2B|cd,JT>. Note MatEle may be empty, we need check when we use.
  void get2BmatAt(int a,int b,int c,int d,int J,int T,TwoBodyMatEle & MatEle) const;

  void set2BmatAt(int a,int b,int c,int d,int J,int T,const TwoBodyMatEle & MatEle);///< set matrix element <ab,J|2B|cd,J>, called by setupTwoBodyMat to write value into the 2Bmat;
};
template<class OrbitalType,class DataType>
void System_T<OrbitalType,DataType>::setupTwoBodyStates()
{
  this->TwoBodyStates.clear();
  for(int i=0;i<this->Orbitals.size();i++)
    {
      for(int j=i;j<this->Orbitals.size();j++)
	{
	  int Jmin=abs(this->Orbitals[i].jj - this->Orbitals[j].jj)/2;
	  int Jmax=(this->Orbitals[i].jj + this->Orbitals[j].jj)/2;	  
	  for(int J=Jmin;J<=Jmax;J++)
	    {
	      for(int T=0;T<=1;T++)
		{
		  int Par=(this->Orbitals[i].l + this->Orbitals[j].l)%2;
		  //if i==j, J+T must be odd. |ab,JT) = (-1)^{ja+jb-J+T} |ba,J)
		  if( (i==j) && (J+T+1)%2) continue;
		  this->TwoBodyStates.push_back( TwoBodyStateType(i,j,T,Par,J) );
		}
	    }
	}
    }
  this->setupTwoBodyChannels();//noting that the twobody states are sorted after this calling.
}

template<class OrbitalType,class DataType>
void System_T<OrbitalType,DataType>::setupIndiceIn2BStates()
{
  int totalOrbitals=this->Orbitals.size();
  IndiceIn2BStates.resize(totalOrbitals);
  for(int i=0;i<totalOrbitals;i++) IndiceIn2BStates[i].resize(totalOrbitals);
  for(int a=0;a<totalOrbitals;a++)
    for(int b=a;b<totalOrbitals;b++)
      {
	IndiceIn2BStates[a][b].resize(2);
	for(int T=0;T<=1;T++)
	  {
	    int Jmin=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
	    int Jmax=(this->Orbitals[a].jj+this->Orbitals[b].jj)/2;
	    IndiceIn2BStates[a][b][T].resize(Jmax-Jmin+1);
	    for(int J=Jmin;J<=Jmax;J++)
	      IndiceIn2BStates[a][b][T][J-Jmin]=-1;//initialize as -1
	  }
      }
  for(int i=0;i<this->TwoBodyStates.size();i++)
    {
      int a=this->TwoBodyStates[i].a;
      int b=this->TwoBodyStates[i].b;
      int J=this->TwoBodyStates[i].J;
      int T=this->TwoBodyStates[i].T;
      int Jmin=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
      IndiceIn2BStates[a][b][T][J-Jmin]=i;
    }
}

template<class OrbitalType,class DataType>
  void System_T<OrbitalType,DataType>::get2BmatAt(int a,int b,int c,int d,int J,int T,TwoBodyMatEle & MatEle) const
{
  MatEle.setZero();
  int phase=1;
  //make sure a<=b,c<=d,(ab)<=(cd)
  if(a>b)
    {
      swap(a,b);
      phase*=( (this->Orbitals[a].jj+this->Orbitals[b].jj)/2-J + T )%2?-1:1;
    }
  if(c>d)
    {
      swap(c,d);
      phase*=( (this->Orbitals[c].jj+this->Orbitals[d].jj)/2-J + T )%2?-1:1;
    }

  int Jminab=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
  int Jmincd=abs(this->Orbitals[c].jj-this->Orbitals[d].jj)/2;
  int Jmaxab=(this->Orbitals[a].jj+this->Orbitals[b].jj)/2;
  int Jmaxcd=(this->Orbitals[c].jj+this->Orbitals[d].jj)/2;

  if( J<max(Jminab,Jmincd) ) return;
  if( J>min(Jmaxab,Jmaxcd) ) return;
  //T belong to (0,1); array bound check when debugging
  int ab=IndiceIn2BStates[a][b][T][J-Jminab];
  int cd=IndiceIn2BStates[c][d][T][J-Jmincd];
  if(ab==-1||cd==-1) return;
  int channelab=this->ChannelIndice[ab];
  int channelcd=this->ChannelIndice[cd];
  //channelab should = channelcd
  if(channelab!=channelcd) return;
  int bra=this->TwoBodyStateIndice[ab];
  int ket=this->TwoBodyStateIndice[cd];
  if(bra>ket)
    swap(bra,ket);
  this->get2Bmat(channelab,bra,ket,MatEle);
  MatEle*=phase;
}

template<class OrbitalType,class DataType>
  void System_T<OrbitalType,DataType>::set2BmatAt(int a,int b,int c,int d,int J,int T,const TwoBodyMatEle & MatEle)
{
  //make sure a<=b,c<=d,bra<=ket
  int phase=1;
  if(a>b)
    {
      swap(a,b);
      phase*=( (this->Orbitals[a].jj+this->Orbitals[b].jj)/2 -J + T )%2?-1:1;
    }
  if(c>d)
    {
      swap(c,d);
      phase*=( (this->Orbitals[c].jj+this->Orbitals[d].jj)/2 -J + T )%2?-1:1;
    }
  int Jminab=abs(this->Orbitals[a].jj-this->Orbitals[b].jj)/2;
  int Jmincd=abs(this->Orbitals[c].jj-this->Orbitals[d].jj)/2;
  int ab=IndiceIn2BStates[a][b][T][J-Jminab];
  int cd=IndiceIn2BStates[c][d][T][J-Jmincd];
  int channelab=this->ChannelIndice[ab];
  int channelcd=this->ChannelIndice[cd];
  //channelab should = channelcd,we should add this error check when debuging
  int bra=this->TwoBodyStateIndice[ab];
  int ket=this->TwoBodyStateIndice[cd];
  if(bra>ket)
    {
      swap(bra,ket);
    }
  TwoBodyMatEle MatEleTemp=MatEle;
  MatEleTemp*=phase;
  this->set2Bmat(channelab,bra,ket,MatEleTemp);
}





//****************************************************************************************

///System under HO. basis, twobody states is in T representation
/*
  dim of TwoBodyMatEle=6,the 6 double are:
      0  (p1-p2)^2/(4m hbar omega)
      1  (p1-p2)^2/(4m habr omega)+ (r1-r2)^2 m omega/(4 hbar), relative hamiltion
      2  coloumb interaction
      3  nuclear interaction for Tz=0
      4  nuclear interaction for Tz=1
      5  nuclear interaction for Tz=-1
*/
class System_T_HO:public System_T<HO_Orbital_T,double>
{
 public:
  typedef HO_Orbital_T OrbitalType;
  typedef double DataType;
  typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
  const static int MatEledim;//< dimension of TwoBodyMatEle.
  System_T_HO(){};
  System_T_HO(const string & _OrbFile,const string & _MatFile):OrbFile(_OrbFile),MatFile(_MatFile){}
  System_T_HO(const string & _OrbFile,const string & _MatFile,double _hbar_omega):OrbFile(_OrbFile),MatFile(_MatFile),hbar_omega(_hbar_omega){}

  void setInputFiles(const string & _OrbFile,const string & _MatFile)
  {
    OrbFile=_OrbFile;
    MatFile=_MatFile;
  }
  void setInput(const string & _OrbFile,const string & _MatFile,const double _hbar_omega)
  {
    OrbFile=_OrbFile;
    MatFile=_MatFile;
    hbar_omega=_hbar_omega;
  }
  /// setup orbitals,2B states and 2B mat.
  void setup()
  {
    setupOrbitals();
    setupTwoBodyStates();
    setupTwoBodyMat();
  }
  
  void setupOrbitals();
  
  void setupTwoBodyMat();
  
  double hbar_omega;
  double lengthb;  
 private:
  string OrbFile;
  string MatFile;
};


//****************************************************************************************
///transform form System_HO_T to System_HO_Tz
void TtoTz(const System_T_HO&,System_Tz_HO&);

///transform form System_HO_Tz to System_HO_T
void TztoT(const System_Tz_HO&,System_T_HO&);

#endif
