#ifndef BLITZ_ARRAY_SERIALIZE_SUPPORT_H
#define BLITZ_ARRAY_SERIALIZE_SUPPORT_H
#include <blitz/array.h>
#include "geocal_time.h"

// This contains various support routines for *implementing* the boost
// serialization (as opposed to *using*  the serialization)
#include "geocal_config.h"
// This adds serialization of blitz::Array.
#ifdef GEOCAL_HAVE_BOOST_SERIALIZATION

// We can extend the size and types supported as needed.
#include <boost/serialization/split_free.hpp>
namespace boost {
  namespace serialization {
    template<class Archive, class T>
    void save(Archive& ar, const blitz::Array<T, 1>& A, 
	      const unsigned version);
    template<typename Archive, class T>
    void load(Archive& ar, blitz::Array<T, 1>& A, 
	      const unsigned version);
    template<class Archive, class T>
    void save(Archive& ar, const blitz::Array<T, 2>& A, 
	      const unsigned version);
    template<typename Archive, class T>
    void load(Archive& ar, blitz::Array<T, 2>& A, 
	      const unsigned version);
    template<class Archive, class T>
    void save(Archive& ar, const blitz::Array<T, 3>& A, 
	      const unsigned version);
    template<typename Archive, class T>
    void load(Archive& ar, blitz::Array<T, 3>& A, 
	      const unsigned version);
    template<class Archive, class T>
    void save(Archive& ar, const blitz::Array<T, 4>& A, 
	      const unsigned version);
    template<typename Archive, class T>
    void load(Archive& ar, blitz::Array<T, 4>& A, 
	      const unsigned version);
    template<class Archive, class T>
    void save(Archive& ar, const blitz::Array<T, 5>& A, 
	      const unsigned version);
    template<typename Archive, class T>
    void load(Archive& ar, blitz::Array<T, 5>& A, 
	      const unsigned version);
  }
}
typedef blitz::Array<double, 1> blitz_double_array_1d;
typedef blitz::Array<double, 2> blitz_double_array_2d;
typedef blitz::Array<double, 3> blitz_double_array_3d;
typedef blitz::Array<double, 4> blitz_double_array_4d;
typedef blitz::Array<double, 5> blitz_double_array_5d;
typedef blitz::Array<float, 1> blitz_float_array_1d;
typedef blitz::Array<float, 2> blitz_float_array_2d;
typedef blitz::Array<float, 3> blitz_float_array_3d;
typedef blitz::Array<float, 4> blitz_float_array_4d;
typedef blitz::Array<float, 5> blitz_float_array_5d;
typedef blitz::Array<int, 1> blitz_int_array_1d;
typedef blitz::Array<int, 2> blitz_int_array_2d;
typedef blitz::Array<int, 3> blitz_int_array_3d;
typedef blitz::Array<int, 4> blitz_int_array_4d;
typedef blitz::Array<int, 5> blitz_int_array_5d;
typedef blitz::Array<bool, 1> blitz_bool_array_1d;
typedef blitz::Array<bool, 2> blitz_bool_array_2d;
typedef blitz::Array<bool, 3> blitz_bool_array_3d;
typedef blitz::Array<bool, 4> blitz_bool_array_4d;
typedef blitz::Array<bool, 5> blitz_bool_array_5d;
typedef blitz::Array<unsigned int, 1> blitz_uint_array_1d;
typedef blitz::Array<unsigned int, 2> blitz_uint_array_2d;
typedef blitz::Array<unsigned int, 3> blitz_uint_array_3d;
typedef blitz::Array<unsigned int, 4> blitz_uint_array_4d;
typedef blitz::Array<unsigned int, 5> blitz_uint_array_5d;
typedef blitz::Array<unsigned short int, 1> blitz_ushort_array_1d;
typedef blitz::Array<unsigned short int, 2> blitz_ushort_array_2d;
typedef blitz::Array<unsigned short int, 3> blitz_ushort_array_3d;
typedef blitz::Array<unsigned short int, 4> blitz_ushort_array_4d;
typedef blitz::Array<unsigned short int, 5> blitz_ushort_array_5d;
typedef blitz::Array<short int, 1> blitz_short_array_1d;
typedef blitz::Array<short int, 2> blitz_short_array_2d;
typedef blitz::Array<short int, 3> blitz_short_array_3d;
typedef blitz::Array<short int, 4> blitz_short_array_4d;
typedef blitz::Array<short int, 5> blitz_short_array_5d;
typedef blitz::Array<GeoCal::Time, 1> blitz_time_array_1d;
typedef blitz::Array<std::string, 1> blitz_string_array_1d;
typedef blitz::Array<std::string, 2> blitz_string_array_2d;
BOOST_SERIALIZATION_SPLIT_FREE(blitz_double_array_1d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_double_array_2d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_double_array_3d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_double_array_4d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_double_array_5d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_string_array_1d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_string_array_2d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_int_array_1d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_int_array_2d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_int_array_3d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_int_array_4d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_int_array_5d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_bool_array_1d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_bool_array_2d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_bool_array_3d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_bool_array_4d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_bool_array_5d);
BOOST_SERIALIZATION_SPLIT_FREE(blitz_time_array_1d);

#endif

#endif
