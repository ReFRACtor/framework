#include "blitz_array_serialize_support.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"

#ifdef FP_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/array.hpp>

template<class Archive, class T>
void boost::serialization::save(Archive& ar, const blitz::Array<T, 1>& A, 
				const unsigned UNUSED(version))
{
  using boost::serialization::make_array;
  if(A.size() > 0 && !A.isStorageContiguous())
    throw FullPhysics::Exception("We can only save contiguous matrix data");
  int size = A.rows();
  ar << FP_NVP(size);
  if(A.size() > 0)
    ar << FP_NVP2("data", make_array(A.data(), A.size()));
}

template<typename Archive, class T>
void boost::serialization::load(Archive& ar, blitz::Array<T, 1>& A, 
				const unsigned UNUSED(version))
{
  using boost::serialization::make_array;
  int size;
  ar >> FP_NVP(size);
  A.resize(size);
  if(A.size() > 0)
    ar >> FP_NVP2("data", make_array(A.data(), A.size()));
}

template<class Archive, class T>
void boost::serialization::save(Archive& ar, const blitz::Array<T, 2>& A, 
				const unsigned UNUSED(version))
{
  using boost::serialization::make_array;
  if(A.size() > 0 && !A.isStorageContiguous())
    throw FullPhysics::Exception("We can only save contiguous matrix data");
  int rows = A.rows();
  int cols = A.cols();
  ar << FP_NVP(rows) << FP_NVP(cols);
  if(A.size() > 0)
    ar << FP_NVP2("data", make_array(A.data(), A.size()));
}

template<typename Archive, class T>
void boost::serialization::load(Archive& ar, blitz::Array<T, 2>& A, 
				const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  int rows, cols;
  ar >> FP_NVP(rows) >> FP_NVP(cols);
  A.resize(rows, cols);
  if(A.size() > 0)
    ar >> FP_NVP2("data", make_array(A.data(), A.size()));
}

template<class Archive, class T>
void boost::serialization::save(Archive& ar, const blitz::Array<T, 3>& A, 
			      const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  if(A.size() > 0 && !A.isStorageContiguous())
    throw FullPhysics::Exception("We can only save contiguous matrix data");
  int rows = A.rows();
  int cols = A.cols();
  int depth = A.depth();
  ar << FP_NVP(rows) << FP_NVP(cols) << FP_NVP(depth);
  if(A.size() > 0)
    ar << FP_NVP2("data", make_array(A.data(), A.size()));
}

template<typename Archive, class T>
void boost::serialization::load(Archive& ar, blitz::Array<T, 3>& A, 
				const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  int rows, cols, depth;
  ar >> FP_NVP(rows) >> FP_NVP(cols) >> FP_NVP(depth);
  A.resize(rows, cols, depth);
  if(A.size() > 0)
    ar >> FP_NVP2("data", make_array(A.data(), A.size()));
}

template<class Archive, class T>
void boost::serialization::save(Archive& ar, const blitz::Array<T, 4>& A, 
			      const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  if(A.size() > 0 && !A.isStorageContiguous())
    throw FullPhysics::Exception("We can only save contiguous matrix data");
  int shp0 = A.shape()(0);
  int shp1 = A.shape()(1);
  int shp2 = A.shape()(2);
  int shp3 = A.shape()(3);
  ar << FP_NVP(shp0) << FP_NVP(shp1) << FP_NVP(shp2)
     << FP_NVP(shp3);
  if(A.size() > 0)
    ar << FP_NVP2("data", make_array(A.data(), A.size()));
}

template<typename Archive, class T>
void boost::serialization::load(Archive& ar, blitz::Array<T, 4>& A, 
				const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  blitz::TinyVector<int, 4> shp;
  ar >> FP_NVP2("shp0", shp(0)) 
     >> FP_NVP2("shp1", shp(1)) 
     >> FP_NVP2("shp2", shp(2)) 
     >> FP_NVP2("shp3", shp(3));
  A.resize(shp);
  if(A.size() > 0)
    ar >> FP_NVP2("data", make_array(A.data(), A.size()));
}

template<class Archive, class T>
void boost::serialization::save(Archive& ar, const blitz::Array<T, 5>& A, 
			      const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  if(A.size() > 0 && !A.isStorageContiguous())
    throw FullPhysics::Exception("We can only save contiguous matrix data");
  int shp0 = A.shape()(0);
  int shp1 = A.shape()(1);
  int shp2 = A.shape()(2);
  int shp3 = A.shape()(3);
  int shp4 = A.shape()(4);
  ar << FP_NVP(shp0) << FP_NVP(shp1) << FP_NVP(shp2)
     << FP_NVP(shp3) << FP_NVP(shp4);
  if(A.size() > 0)
    ar << FP_NVP2("data", make_array(A.data(), A.size()));
}

template<typename Archive, class T>
void boost::serialization::load(Archive& ar, blitz::Array<T, 5>& A, 
				const unsigned UNUSED(version)) 
{
  using boost::serialization::make_array;
  blitz::TinyVector<int, 5> shp;
  ar >> FP_NVP2("shp0", shp(0)) 
     >> FP_NVP2("shp1", shp(1)) 
     >> FP_NVP2("shp2", shp(2)) 
     >> FP_NVP2("shp3", shp(3))
     >> FP_NVP2("shp4", shp(4));
  A.resize(shp);
  if(A.size() > 0)
    ar >> FP_NVP2("data", make_array(A.data(), A.size()));
}

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<double, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<double, 1>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<double, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<double, 2>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<double, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<double, 3>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<double, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<double, 4>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<double, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<double, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<FullPhysics::AutoDerivative<double>, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<FullPhysics::AutoDerivative<double>, 1>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<FullPhysics::AutoDerivative<double>, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<FullPhysics::AutoDerivative<double>, 2>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<FullPhysics::AutoDerivative<double>, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<FullPhysics::AutoDerivative<double>, 3>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<FullPhysics::AutoDerivative<double>, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<FullPhysics::AutoDerivative<double>, 4>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<FullPhysics::AutoDerivative<double>, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<FullPhysics::AutoDerivative<double>, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<float, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<float, 1>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<float, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<float, 2>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<float, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<float, 3>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<float, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<float, 4>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<float, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<float, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<std::string, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<std::string, 1>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<std::string, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<std::string, 2>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<int, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<int, 1>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<int, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<int, 2>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<int, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<int, 3>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<int, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<int, 4>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<int, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<int, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<bool, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<bool, 1>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<bool, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<bool, 2>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<bool, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<bool, 3>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<bool, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<bool, 4>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<bool, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<bool, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned int, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned int, 1>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned int, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned int, 2>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned int, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned int, 3>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned int, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned int, 4>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned int, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned int, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<short int, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<short int, 1>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<short int, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<short int, 2>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<short int, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<short int, 3>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<short int, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<short int, 4>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<short int, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<short int, 5>& A, 
					 const unsigned version);

template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned short int, 1>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned short int, 1>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned short int, 2>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned short int, 2>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned short int, 3>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned short int, 3>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned short int, 4>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned short int, 4>& A, 
					 const unsigned version);
template void boost::serialization::load(polymorphic_iarchive& ar, 
					 blitz::Array<unsigned short int, 5>& A, 
					 const unsigned version);
template void boost::serialization::save(polymorphic_oarchive& ar, 
					 const blitz::Array<unsigned short int, 5>& A, 
					 const unsigned version);

#endif
