#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>

#include "example_level_1b.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "example_observation_id.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ExampleLevel1b::serialize(Archive & ar,
			       const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Level1bSampleCoefficient)
    & FP_NVP_(input) & FP_NVP_(data_index);
}

FP_IMPLEMENT(ExampleLevel1b);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ExampleLevel1b, Level1bSampleCoefficient)
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&, const std::string&>())
REGISTER_LUA_END()
#endif

ExampleLevel1b::ExampleLevel1b
(const boost::shared_ptr<HdfFile>& input_file,
 const std::string& observation_id)
:input_(input_file)
{
  ExampleObservationId<std::string> obs_id(input_file, "observation_ids",
					   observation_id);
  data_index_ = obs_id.data_index();
}

ExampleLevel1b::ExampleLevel1b
(const std::string& input_filename,
 const std::string& observation_id)
: input_(boost::make_shared<HdfFile>(input_filename))
{
  ExampleObservationId<std::string> obs_id(input_, "observation_ids",
					   observation_id);
  data_index_ = obs_id.data_index();
}

int ExampleLevel1b::number_spectrometer() const
{
  TinyVector<int, 2> lat_shape =
    input_->read_shape<2>(group_name + "/latitude");
  return lat_shape(1);
}

int ExampleLevel1b::number_sample(int Spec_index) const
{
  std::string rad_ds_name = group_name + "/radiance_" +
    boost::lexical_cast<std::string>(Spec_index + 1);
  TinyVector<int, 2> rad_shape = input_->read_shape<2>(rad_ds_name);
  return rad_shape(1);
}

blitz::Array<double, 1> ExampleLevel1b::spectral_variable
(int sensor_index) const
{ 
  blitz::Array<double, 1> var_vals(number_sample(sensor_index));
  blitz::firstIndex i1; 
  var_vals = i1 + 1;
  return var_vals;
}

SpectralRange ExampleLevel1b::radiance(int Spec_index) const
{
  std::string rad_ds_name = group_name + "/radiance_" +
    boost::lexical_cast<std::string>(Spec_index + 1);
  std::string uncert_ds_name = group_name + "/uncertainty_" +
    boost::lexical_cast<std::string>(Spec_index + 1);

  TinyVector<int, 2> rad_ds_shape = input_->read_shape<2>(rad_ds_name);
  TinyVector<int, 2> rad_start, rad_size;
  rad_start = data_index_, 0;
  rad_size  = 1, rad_ds_shape(1);
  ArrayWithUnit<double, 2> radiance =
    input_->read_field_with_unit<double, 2>(rad_ds_name, 
        Unit("Ph sec^{-1} m^{-2} sr^{-1} um^{-1}"),
        rad_start, rad_size);

  if (input_->has_object(uncert_ds_name)) {
    TinyVector<int, 2> uncert_ds_shape = input_->read_shape<2>(uncert_ds_name);
    TinyVector<int, 2> uncert_start, uncert_size;
    uncert_start = data_index_, 0;
    uncert_size  = 1, uncert_ds_shape(1);
    Array<double, 2> uncertainty =
      input_->read_field<double, 2>(uncert_ds_name, uncert_start, uncert_size);

    return SpectralRange(radiance(0, Range::all()).value, radiance.units,
			 uncertainty(0, Range::all()));
  } else
    return SpectralRange(radiance(0, Range::all()).value, radiance.units);
}

DoubleWithUnit ExampleLevel1b::read_scalar_with_unit
(const std::string& dataset_name, int i, const Unit& default_unit) const
{
  range_check(i, 0, number_spectrometer());
  TinyVector<int, 2> start, size;
  start = data_index_, i;
  size = 1, 1;
  ArrayWithUnit<double, 2> val_arr =
    input_->read_field_with_unit<double, 2>(dataset_name, default_unit,
					    start, size);
  return val_arr(0, 0);
}

double ExampleLevel1b::read_scalar
(const std::string& dataset_name, int i) const
{
  range_check(i, 0, number_spectrometer());
  TinyVector<int, 2> start, size;
  start = data_index_, i;
  size = 1, 1;
  Array<double, 2> val_arr = input_->read_field<double, 2>(dataset_name,
							   start, size);
  return val_arr(0, 0);
}

ArrayWithUnit<double, 1> ExampleLevel1b::read_array_with_unit
(const std::string& dataset_name, int i, const Unit& default_unit) const
{
  range_check(i, 0, number_spectrometer());
  
  TinyVector<int, 3> ds_shape = input_->read_shape<3>(dataset_name);
  TinyVector<int, 3> start, size;
  start = data_index_, i, 0;
  size = 1, 1, ds_shape(2);

  ArrayWithUnit<double, 3> val_arr =
    input_->read_field_with_unit<double, 3>(dataset_name, default_unit,
					    start, size);

  ArrayWithUnit<double, 1> ret_val;
  ret_val.units = val_arr.units;
  ret_val.value.resize(ds_shape(2));
  ret_val.value = val_arr.value(0, 0, Range::all());

  return ret_val;
}

Array<double, 1> ExampleLevel1b::read_array
(const std::string& dataset_name, int i) const
{
  range_check(i, 0, number_spectrometer());
  TinyVector<int, 3> ds_shape = input_->read_shape<3>(dataset_name);
  TinyVector<int, 3> start, size;
  start = data_index_, i, 0;
  size = 1, 1, ds_shape(2);
  Array<double, 3> val_arr = input_->read_field<double, 3>(dataset_name,
							   start, size);
  return val_arr(0, 0, Range::all());
}
