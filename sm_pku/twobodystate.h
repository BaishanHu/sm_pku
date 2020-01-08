#ifndef _TWOBODYSTATE_H_
#define _TWOBODYSTATE_H_
#include"base.h"
///twobody state with quantum number isospin projection, parity and angular momentum
class TwoBodyState_TzParJ:public TwoBodyState
{
 public:
  TwoBodyState_TzParJ(){}
  TwoBodyState_TzParJ(int _a,int _b,int _Tz,int _Par,int _J):TwoBodyState(_a,_b),
                                                         Tz(_Tz),Par(_Par),J(_J){};
  int Tz,Par,J;///< Tz=-1,0,1 for pp,pn,nn; Par=0 for positive,1 for negative
  bool operator < (const TwoBodyState&) const;
  bool isInAChannelWith(const TwoBodyState&) const;
};
//*********************************************************

///twobody state with quantum number isospin projection, parity, angular momentum projection
class TwoBodyState_TzParM:public TwoBodyState
{
 public:
  TwoBodyState_TzParM(){}
 TwoBodyState_TzParM(int _a,int _b,int _Tz,int _Par,int _M):TwoBodyState(_a,_b),
    Tz(_Tz),Par(_Par),M(_M){};
  int Tz,Par,M;///< Tz=-1,0,1 for pp,pn,nn; Par=0 for positive,1 for negative
  bool operator < (const TwoBodyState&) const;
  bool isInAChannelWith(const TwoBodyState&) const;
};

//*********************************************************

class TwoBodyState_TParJ:public TwoBodyState
{
 public:
  TwoBodyState_TParJ(){}
  TwoBodyState_TParJ(int _a,int _b,int _T,int _Par,int _J):TwoBodyState(_a,_b),T(_T),Par(_Par),J(_J){}
  bool operator < (const TwoBodyState&) const;
  bool isInAChannelWith(const TwoBodyState&) const;
  int T,Par,J;///< Par=0 for positive,1 for negative
};

#endif
