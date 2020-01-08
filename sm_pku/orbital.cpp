#include"orbital.h"
bool Orbital_ljjtz::isInAGroupWith(const Orbital&orb) const
{
 const Orbital_ljjtz * OrbPointer = (Orbital_ljjtz*)&orb;
  return (l == OrbPointer->l)&&
         (jj == OrbPointer->jj)&&
         (tz == OrbPointer->tz);
}

//*******************************************************
bool Orbital_parmmtz::isInAGroupWith(const Orbital&orb) const
{
  const Orbital_parmmtz * OrbPointer = (Orbital_parmmtz*)&orb;
  return (par == OrbPointer->par)&&
         (mm == OrbPointer->mm)&&
         (tz == OrbPointer->tz);
}

//*******************************************************

bool Orbital_ljj::isInAGroupWith(const Orbital&orb) const
{
 const Orbital_ljj * OrbPointer = (Orbital_ljj*)&orb;
  return (l == OrbPointer->l)&&
         (jj == OrbPointer->jj);
}

//*******************************************************
