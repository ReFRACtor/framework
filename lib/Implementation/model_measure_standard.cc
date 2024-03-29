#include <model_measure_standard.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ModelMeasureStandard::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasure)
    & FP_NVP(fm) & FP_NVP(obs) & FP_NVP(sv);
}

FP_IMPLEMENT(ModelMeasureStandard);
#endif


// Function to append a second array at the bottom of an existing
// array.

template<class T> void append_array(blitz::Array<T, 1>& A, const blitz::Array<T, 1>& A_to_append)
{
  Range r_to(A.rows(), A.rows() + A_to_append.rows()-1);
  A.resizeAndPreserve(A.rows() + A_to_append.rows());
  A(r_to) = A_to_append(Range::all());
}

template<class T> void append_array(blitz::Array<T, 2>& A, const blitz::Array<T, 2>& A_to_append)
{
  if(A.cols() != A_to_append.cols())
    throw Exception("Matrices need to have the same number of columns");
  Range r_to(A.rows(), A.rows() + A_to_append.rows()-1);
  A.resizeAndPreserve(A.rows() + A_to_append.rows(), A.cols());
  A(r_to, Range::all()) = A_to_append(Range::all(), Range::all());
}


ModelMeasureStandard::ModelMeasureStandard
(const boost::shared_ptr<ForwardModel>& forward_model, 
 const boost::shared_ptr<Observation>& observation, 
 const boost::shared_ptr<StateVector>& state_vector)
  : sv(state_vector)
{
  fm.push_back(forward_model);
  obs.push_back(observation);
  SpectralRange s0 = obs[0]->radiance_all().spectral_range();
  set_measurement(s0.data(),
		  Array<double, 1>(sqr(s0.uncertainty())));
}

ModelMeasureStandard::ModelMeasureStandard
(const std::vector<boost::shared_ptr<ForwardModel> >& forward_model, 
 const std::vector<boost::shared_ptr<Observation> >& observation, 
 const boost::shared_ptr<StateVector>& state_vector)
  : fm(forward_model), obs(observation), sv(state_vector)
{
  if(forward_model.size() < 1)
    throw Exception("Need to have at least one forward model.");
  if(forward_model.size() != observation.size())
    throw Exception("forward_model and observation vectors need to be the same size");
  
  SpectralRange s0 = obs[0]->radiance_all().spectral_range();
  Array<double, 1> a1(s0.data());
  Array<double, 1> a2(sqr(s0.uncertainty()));
  for(int i = 1; i < (int) obs.size(); ++i) {
    SpectralRange s1 = obs[i]->radiance_all().spectral_range();
    append_array(a1, s1.data());
    Array<double, 1> a3(sqr(s1.uncertainty()));
    append_array(a2, a3);
  }
  set_measurement(a1, a2);
}


void ModelMeasureStandard::model_eval()
{
  if(M.size() <= 0)
    radiance_from_fm();
}


void ModelMeasureStandard::jacobian_eval()
{
  if(K.size() <= 0)
    radiance_from_fm();
}


void ModelMeasureStandard::model_jacobian_eval()
{
  if((K.size() <= 0) or (M.size() <= 0))
    radiance_from_fm();
}

int ModelMeasureStandard::expected_parameter_size() const
{
  return sv->observer_claimed_size();
}


void ModelMeasureStandard::radiance_from_fm()
{
  assert_parameter_set_correctly();

  SpectralRange s0 = obs[0]->radiance_all().spectral_range();
  Array<double, 1> temp_msrmnt(s0.data());
  Spectrum rad_spec = fm[0]->radiance_all(false);
  SpectralRange rad_mod = rad_spec.spectral_range().convert(s0.units());
  M.reference(rad_mod.data_ad().value());
  K.reference(rad_mod.data_ad().jacobian());
  for(int i = 1; i < (int) obs.size(); ++i) {
    SpectralRange s1 = obs[i]->radiance_all().spectral_range();
    append_array(temp_msrmnt, s1.data());
    Spectrum s2 = fm[i]->radiance_all(false);
    SpectralRange s3 = s2.spectral_range().convert(s1.units());
    append_array(M, s3.data_ad().value());
    append_array(K, s3.data_ad().jacobian());
  }

  //  TEMPORARY
  //
  //  There should be no need for this error check block.
  //  The problem is that with our current model code 
  //  measurement may change.
  //
  if(msrmnt.rows() != temp_msrmnt.rows())
    throw Exception("Measurement has changed during the retrieval. :( ");
  double msrmnt_L1_norm = sum(abs(msrmnt));
  if(msrmnt_L1_norm <= 0.0)
    throw Exception("Measurement has no signal (just 0). :( ");
  if( sum(abs(temp_msrmnt-msrmnt))/msrmnt_L1_norm > 0.0000001 )
    throw Exception("Measurement has changed during the retrieval. :( ");

  assert_model_correct(M);
  M.makeUnique();
  assert_jacobian_correct(K);
  K.makeUnique();
}

void ModelMeasureStandard::parameters(const blitz::Array<double, 1>& x)
{
  //  This check is needed because it is not obvious 
  //  whether or not the Forward Model does caching.
  //
  if(!parameters_different(x)) return;
  ModelMeasure::parameters(x);
  sv->update_state(x);
}
