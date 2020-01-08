#ifndef _BASE_H_
#define _BASE_H_
#include<vector>
#include<string>
#include<algorithm>
#include<iostream>

using std::vector;
using std::string;
using std::swap;
using std::sort;
using std::ostream;
using std::endl;

/// base class for orbital
class Orbital
{
 public:
  virtual ~Orbital(){};
  ///whether it has same quantum number with
  virtual bool isInAGroupWith(const Orbital&) const=0;
};




//sepration line
/*********************************************************************/


///  base class for 2B states
class TwoBodyState
{
 public:
  TwoBodyState(){}
  TwoBodyState(int _a,int _b):a(_a),b(_b){}
  virtual ~TwoBodyState(){};
  int a,b;

  virtual bool isInAChannelWith(const TwoBodyState&) const=0;
  
  /// used for sort ...
  virtual bool operator<(const TwoBodyState&) const=0;
};


 
/*********************************************************************/



///twobody matrix element, may including many twobody operators, thus inherit form class vector.
template <class DataType>
class TwoBodyMatElement:public vector<DataType>
{
 public:
 TwoBodyMatElement(int length=0):vector<DataType>(length){}
  void setZero()
  {
    for(int i=0;i<this->size();i++)
      {
	this->at(i)=0;
      }
  }
  void operator *=(DataType fac)
  {
    for(int i=0;i<this->size();i++)
      {
	this->at(i)*=fac;
      }
  }
  TwoBodyMatElement  operator *(const DataType & fac)
    {
      TwoBodyMatElement v(*this);
      v*=fac;
      return v;
    }
  TwoBodyMatElement & operator +=(const TwoBodyMatElement &TwoB)
    {
      for(int i=0;i<this->size();i++)
	{
	  this->at(i)+=TwoB.at(i);
	}
      return *this;
    }


  friend ostream & operator <<(ostream & os,const TwoBodyMatElement<DataType> & MatEle)
  {
    for(int i=0;i<MatEle.size();i++)
      os<<MatEle[i]<<"\t";
    os<<endl;
    return os;
  }
};


/*********************************************************************/


/// base system class for many body system
template <class OrbitalType, class TwoBodyStateType,class DataType>
class System
{
 public:
  typedef TwoBodyMatElement<DataType> TwoBodyMatEle;
  
  virtual ~System(){}
  virtual void setupOrbitals()=0;
  vector<OrbitalType> Orbitals;
  void setupGroups();
  vector< vector<int> > Groups;///<first indice for different groups,second indice for different orbitals in a group.
  vector<int> GroupIndice;///< GroupIndice[i] return the first indice of ith Orbital in "Groups"
  vector<int> OrbIndice;///< OrbIndice[i] return the second indice of ith Orbital in "Groups", Groups[GroupIndice[i]][OrbIndice[i]]=i
  vector< vector<DataType> > OneBodyMat;


  virtual void setupTwoBodyStates()=0;
  vector<TwoBodyStateType> TwoBodyStates;
  void setupTwoBodyChannels();
  vector< vector<int> > TwoBodyChannels;///< TwoBodyChannels[a][b] return the indice in TwoBodyStates of a th channel b th state

  vector<int> ChannelIndice;///< ChannelIndice[i] return the channel number corresponding to the i th twobodystate in array TwoBodyStates.
  vector<int> TwoBodyStateIndice;///< TwoBodyStateIndice[i] return the  second  indice in TwoBodyChannels corresponding to the i th twobodystate in array TwoBodyStates.

  virtual void setupTwoBodyMat()=0;
  void set2Bmat(int channel,int bra,int ket,const TwoBodyMatEle & MatEle);///< set matrix element TwoBodymat[chan][bra*n+ket-bra*(bra+1)/2], make sure bra<=ket
  void get2Bmat(int channel,int bra,int ket,TwoBodyMatEle & MatEle) const;///< get matrix element TwoBodymat[chan][bra*n+ket-bra*(bra+1)/2], make sure bra<=ket
  vector< vector< TwoBodyMatEle > > TwoBodyMat;///< TwoBodymat[chan][bra*n+ket-bra*(bra+1)/2] return (bra,ket) matrix element of chan th channel. Noting that bra<=ket. n is total number of 2B states in this channel.

};
template <class OrbitalType, class TwoBodyStateType,class DataType>
void System<OrbitalType,TwoBodyStateType,DataType>::setupGroups()
{
  Groups.clear();
  int NumOrbitals=Orbitals.size();
  GroupIndice.resize(NumOrbitals);
  OrbIndice.resize(NumOrbitals);
  for(int i=0;i<NumOrbitals;i++)
    {
      int j=0;
      for(j=0;j<Groups.size();++j)
	{
	  if( Groups[j].empty() ||
	      ( Orbitals[i].isInAGroupWith(Orbitals[ Groups[j][0] ]) )
	      )
	    {
	      Groups[j].push_back(i);
	      GroupIndice[i]=j;
	      OrbIndice[i]=Groups[j].size()-1;
	      break;
	    }
	}
      if(j>=Groups.size())
	{
	  Groups.push_back(vector<int>() );
	  Groups[Groups.size()-1].push_back(i);
	  GroupIndice[i]=Groups.size()-1;
	  OrbIndice[i]=0;
	}
    } 
}

template <class OrbitalType,class TwoBodyStateType,class DataType>
  void System<OrbitalType,TwoBodyStateType,DataType>::setupTwoBodyChannels()
{
  sort(TwoBodyStates.begin(),TwoBodyStates.end());
  ChannelIndice.resize(TwoBodyStates.size());
  TwoBodyStateIndice.resize(TwoBodyStates.size());
  TwoBodyChannels.clear();
  int NumTwoBodyStates=TwoBodyStates.size();
  for(int i=0;i<NumTwoBodyStates;i++)
    {
      int LastChannel=TwoBodyChannels.size()-1;
      //TwoBodyChannels is empty, so push_back a new channel.
      if(LastChannel<0)
	{
	  TwoBodyChannels.push_back( vector<int>(1,i) );
	  ChannelIndice[i]=0;
	  TwoBodyStateIndice[i]=0;
	}
      else
	{
	  int indice=TwoBodyChannels[LastChannel][0];
	  TwoBodyState & Temp=TwoBodyStates[indice];
	  //the twobody state belong to LastChannel, so push_back the twobody state in this channel
	  if( TwoBodyStates[i].isInAChannelWith(Temp) )
	    {
	      TwoBodyChannels[LastChannel].push_back(i);
	      ChannelIndice[i]=LastChannel;
	      TwoBodyStateIndice[i]=TwoBodyChannels[LastChannel].size()-1;
	    }
	  else//don't belong to, push_back a new channel
	    {
	      TwoBodyChannels.push_back( vector<int>(1,i) );
	      ChannelIndice[i]=LastChannel+1;
	      TwoBodyStateIndice[i]=0;
	    }
	}
    }
}

template <class OrbitalType,class TwoBodyStateType,class DataType>
void System<OrbitalType,TwoBodyStateType,DataType>::set2Bmat(int channel,int bra,int ket,const TwoBodyMatEle & MatEle)
{
  int n=TwoBodyChannels[channel].size();
  TwoBodyMat[channel][bra*n+ket-(bra*(bra+1))/2]=MatEle;
}

template <class OrbitalType,class TwoBodyStateType,class DataType>
  void System<OrbitalType,TwoBodyStateType,DataType>::get2Bmat(int channel,int bra,int ket,TwoBodyMatEle & MatEle) const
{
  int n=TwoBodyChannels[channel].size();
  MatEle=TwoBodyMat[channel][bra*n+ket-(bra*(bra+1))/2];
}

#endif
