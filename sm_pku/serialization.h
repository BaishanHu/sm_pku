//serialization of SM classes for mpi communication use
//boost library is used for serialization and mpi parallel
#ifndef _SERIALIZATION_H_ 
#define _SERIALIZATION_H_
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/complex.hpp>

//#include <boost/serialization/bitset.hpp>
#include<boost/mpi.hpp>

#include"sm_system.h"
#include"base.h"
#include"eigen.hpp"
using namespace SM;

namespace boost {

  namespace serialization {

    template<class Archive>
      void serialize(Archive & ar, Orbital_SM_Tz & Orb, const unsigned int version)
      {
	ar & Orb.l;
	ar & Orb.jj;
	ar & Orb.tz;
	ar & Orb.e;	
      }

    template<class Archive>
      void serialize(Archive & ar, Orbital_SM_Tz_M & Orbm, const unsigned int version)
      {
	ar & Orbm.par;
	ar & Orbm.l;
	ar & Orbm.jj;
	ar & Orbm.mm;
	ar & Orbm.tz;
	ar & Orbm.e;
      }

    
    template<class Archive,class DataType>
      void serialize(Archive & ar, TwoBodyMatElement<DataType> & ME, const unsigned int version)
    {
      ar & boost::serialization::base_object< vector<DataType> >(ME);
    }
       

    template<class Archive>
      void serialize(Archive & ar, TwoBodyState_TzParJ & TwoBState, const unsigned int version)
      {
	ar & TwoBState.Tz;
	ar & TwoBState.Par;
	ar & TwoBState.J;
	ar & TwoBState.a;
	ar & TwoBState.b;	
      }

    template<class Archive>
      void serialize(Archive & ar, TwoBodyState_TzParM & TwoBState, const unsigned int version)
      {
	ar & TwoBState.Tz;
	ar & TwoBState.Par;
	ar & TwoBState.M;
	ar & TwoBState.a;
	ar & TwoBState.b;	
      }


    ///a self-implemention of bitset serialization
    template<size_t N>
      void bitset_to_bytes(const std::bitset<N>& bs,std::vector<unsigned char> & bits )
      {
	bits.resize((N + 7) >> 3);
	for (int j=0; j<int(N); j++)
	  bits[j>>3] |= (bs[j] << (j & 7));
      }

    template<size_t N>
      void bitset_from_bytes(const std::vector<unsigned char>& buf,std::bitset<N>& bs)
      {
	assert(buf.size() == ((N + 7) >> 3));
	for (int j=0; j<int(N); j++)
	  bs[j] = ((buf[j>>3] >> (j & 7)) & 1);
      }
    
    template <class Archive, std::size_t size>
      inline void save(
		       Archive & ar,
		       std::bitset<size> const & t,
		       const unsigned int  version 
		       ){
      std::vector<unsigned char> bits;
      bitset_to_bytes(t,bits);
      ar << BOOST_SERIALIZATION_NVP( bits );
    }

    template <class Archive, std::size_t size>
      inline void load(
		       Archive & ar,
		       std::bitset<size> & t,
		       const unsigned int  version 
		       ){
      std::vector<unsigned char> bits;
      ar >> BOOST_SERIALIZATION_NVP( bits );
      bitset_from_bytes(bits,t);
    }

    template <class Archive, std::size_t size>
      inline void serialize(
			    Archive & ar,
			    std::bitset<size> & t,
			    const unsigned int version
			    ){
      boost::serialization::split_free( ar, t, version );
    }

    // don't track bitsets since that would trigger tracking
    // all over the program - which probably would be a surprise.
    // also, tracking would be hard to implement since, we're
    // serialization a representation of the data rather than
    // the data itself.
    template <std::size_t size>
      struct tracking_level<std::bitset<size> >
      : mpl::int_<track_never> {} ;
    ///finish bieset serialization





    
    template<class Archive>
      void serialize(Archive & ar, Configs & configs, const unsigned int version)
      {
	ar & configs.configs;
	ar & configs.maxNeutrons;
	ar & configs.maxProtons;
	ar & configs.vN;
	ar & configs.vZ;
	ar & configs.dim;
      }
    
       
    template<class Archive>
      void serialize(Archive & ar, System_Tz_SM_M & Vm, const unsigned int version)
      {
	ar & Vm.pSystem;
	ar & Vm.maxProtons;
	ar & Vm.maxNeutrons;
	ar & Vm.IndiceIn2BStates;
	ar & Vm.Orders;
	ar & Vm.Orbitals;
	ar & Vm.Groups;
	ar & Vm.GroupIndice;
	ar & Vm.OrbIndice;
	ar & Vm.TwoBodyStates;
	ar & Vm.TwoBodyChannels;
	ar & Vm.ChannelIndice;
	ar & Vm.TwoBodyStateIndice;
	ar & Vm.TwoBodyMat;
      }

    template<class Archive>
      void serialize(Archive & ar, System_Tz_SM & V, const unsigned int version)
      {
	ar & V.IndiceIn2BStates;
	ar & V.Orbitals;
	ar & V.Groups;
	ar & V.GroupIndice;
	ar & V.OrbIndice;
	ar & V.TwoBodyStates;
	ar & V.TwoBodyChannels;
	ar & V.ChannelIndice;
	ar & V.TwoBodyStateIndice;
	ar & V.TwoBodyMat;
        ar & V.OneBodyMat;
	ar & V.lambda1B;
	ar & V.lambda2B;
	ar & V.Oper1B;
	ar & V.Oper2B;
      }

    
  }
}

#endif
