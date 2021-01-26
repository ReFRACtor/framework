#include "fp_serialize_support.h"
#include "boost_rational_serialize_support.h"
#include "generic_object.h"
#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive, class T>
void boost::serialization::save(Archive& ar, const boost::rational<T>& R, 
				const unsigned int UNUSED(version))
{
  int n = R.numerator(), d = R.denominator();
  ar & FP_NVP2("numerator", n) & FP_NVP2("denominator", d);
}

template<class Archive, class T>
void boost::serialization::load(Archive& ar, boost::rational<T>& R,
				const unsigned int UNUSED(version))
{
  int n, d;
  ar & FP_NVP2("numerator", n) & FP_NVP2("denominator", d);
  R = boost::rational<T>(n, d);
}

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 boost::rational<int>& R, 
					 const unsigned int version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const boost::rational<int>& R, 
					 const unsigned int version);
#endif
