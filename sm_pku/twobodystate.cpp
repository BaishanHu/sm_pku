#include"twobodystate.h"
bool TwoBodyState_TzParJ::operator < (const TwoBodyState&state) const
{
  const TwoBodyState_TzParJ* Pointer=(TwoBodyState_TzParJ*)&state;
  if(Tz != Pointer->Tz)
    return Tz < Pointer->Tz;
  else if(Par != Pointer->Par)
    return Par < Pointer->Par;
  else if(J != Pointer->J)
    return J < Pointer->J;
  else if(a != Pointer->a)
    return a < Pointer->a;
  else if(b != Pointer->b)
    return b < Pointer->b;
  return false;
}

bool TwoBodyState_TzParJ::isInAChannelWith(const TwoBodyState&state) const
{
  const TwoBodyState_TzParJ* Pointer=(TwoBodyState_TzParJ*)&state;
  return  Tz == Pointer->Tz&&
    Par == Pointer->Par&&
    J == Pointer->J;
}

//**********************************************************************
bool TwoBodyState_TzParM::operator < (const TwoBodyState&state) const
{
  const TwoBodyState_TzParM* Pointer=(TwoBodyState_TzParM*)&state;
  if(Tz != Pointer->Tz)
    return Tz < Pointer->Tz;
  else if(Par != Pointer->Par)
    return Par < Pointer->Par;
  else if(M != Pointer->M)
    return M < Pointer->M;
  else if(a != Pointer->a)
    return a < Pointer->a;
  else if(b != Pointer->b)
    return b < Pointer->b;
  return false;
}

bool TwoBodyState_TzParM::isInAChannelWith(const TwoBodyState&state) const
{
  const TwoBodyState_TzParM* Pointer=(TwoBodyState_TzParM*)&state;
  return  Tz == Pointer->Tz&&
    Par == Pointer->Par&&
    M == Pointer->M;
}

//**********************************************************************
bool TwoBodyState_TParJ::operator < (const TwoBodyState&state) const
{
  const TwoBodyState_TParJ* Pointer=(TwoBodyState_TParJ*)&state;
  if(T != Pointer->T)
    return T < Pointer->T;
  else if(Par != Pointer->Par)
    return Par < Pointer->Par;
  else if(J != Pointer->J)
    return J < Pointer->J;
  else if(a != Pointer->a)
    return a < Pointer->a;
  else if(b != Pointer->b)
    return b < Pointer->b;
  return false;
}

bool TwoBodyState_TParJ::isInAChannelWith(const TwoBodyState&state) const
{
  const TwoBodyState_TParJ* Pointer=(TwoBodyState_TParJ*)&state;
  return  T == Pointer->T&&
    Par == Pointer->Par&&
    J == Pointer->J;
}
