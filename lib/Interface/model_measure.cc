#include "model_measure.h"
#include "fp_exception.h"
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ModelMeasure::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelState)
    & FP_NVP(msrmnt_does_not_change) 
    & FP_NVP(msrmnt) & FP_NVP(msrmnt_jacobian) & FP_NVP(Se);
}

FP_IMPLEMENT(ModelMeasure);
#endif

void ModelMeasure::set_measurement( const blitz::Array<double, 1>& measurement, 
                                    const blitz::Array<double, 1>& measurement_error_cov )
{
  if(measurement.rows() <= 0)
    throw Exception("The size of measurement data array is zero.");
  if(measurement_error_cov.rows() != measurement.rows() )
    throw Exception("Measured data and the related diagonal error covariance matrix need to be the same size.");

  msrmnt.resize(measurement.shape());
  msrmnt = measurement;

  msrmnt_jacobian.free();
  
  Se.resize(measurement_error_cov.shape());
  Se = measurement_error_cov;
}


void ModelMeasure::assert_model_correct(const blitz::Array<double, 1>& m) const
{
  if( m.rows() != msrmnt.rows() )
    throw Exception("Model data and measured data need to be the same size.");
}


void ModelMeasure::assert_jacobian_correct(const blitz::Array<double, 2>& k) const
{
  if( k.rows() != msrmnt.rows() )
    throw Exception("Model Jacobian and measurement error covariance matrix must have the same number of rows.");
  if( k.cols() != expected_parameter_size() )
    throw Exception("The number of model Jacobian columns must equal the number of parameters.");
}



Array<double, 2> ModelMeasure::uncert_weighted_jacobian()
{ 
  firstIndex i1; secondIndex i2;
  Array<double, 2> result(measurement_size(), parameter_size());
  Array<double, 2> mjac(measurement_jacobian());
  Array<double, 2> jac(jacobian());
  if(mjac.size() > 0) {
    if(jac.rows() != mjac.rows() ||
       jac.cols() != mjac.cols())
      throw Exception("Model jacobian and measurement jacobian need to be the same size");
    result = (jac(i1, i2) - mjac(i1, i2)) / sqrt(Se(i1));
  } else {
    result = jac(i1, i2) / sqrt(Se(i1));
  }
  return result;
}


Array<double, 2> ModelMeasure::uncert_weighted_jac_inner_product()
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> K(jacobian());
  Array<double, 2> result(K.cols(), K.cols());
  result = sum(K(i3,i1) / Se(i3) * K(i3, i2), i3);
  return result;
}
