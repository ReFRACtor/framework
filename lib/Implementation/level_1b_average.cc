#include "level_1b_average.h"
#include <boost/foreach.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Level1bAverage, Level1b)
.def(luabind::constructor<const std::vector<boost::shared_ptr<Level1b> > &>())
REGISTER_LUA_END()
#endif


int Level1bAverage::number_spectrometer() const
{
    assert_field_equal(&Level1b::number_spectrometer);
    return  l1b[0]->number_spectrometer();
}

DoubleWithUnit Level1bAverage::relative_velocity(int Spec_index) const
{
    assert_field_equal(&Level1b::relative_velocity, Spec_index);
    return l1b[0]->relative_velocity(Spec_index);
}

Time Level1bAverage::time(int Spec_index) const
{
    // TODO: Enable assert
    // assert_field_equal(&Level1b::time, Spec_index);
    return l1b[0]->time(Spec_index);
}

SpectralDomain Level1bAverage::sample_grid(int Spec_index) const
{
    // TODO: Enable assert
    // assert_field_equal(&Level1b::sample_grid, Spec_index);
    return l1b[0]->sample_grid(Spec_index);
}

DoubleWithUnit Level1bAverage::latitude(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->latitude(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->latitude(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::longitude(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->longitude(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->longitude(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::sounding_zenith(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->sounding_zenith(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->sounding_zenith(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::sounding_azimuth(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->sounding_azimuth(1).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->sounding_azimuth(i);
  return sum / ((int) l1b.size());
}

Array<double, 1> Level1bAverage::stokes_coefficient(int i) const
{
  Array<double, 1> res(l1b[0]->stokes_coefficient(i).rows());
  res = 0;
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    res += f->stokes_coefficient(i);
  res /= (int) l1b.size();
  return res;
}

DoubleWithUnit Level1bAverage::solar_zenith(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->solar_zenith(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->solar_zenith(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::solar_azimuth(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->solar_azimuth(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->solar_azimuth(i);
  return sum / ((int) l1b.size());
}

DoubleWithUnit Level1bAverage::altitude(int i) const
{
  DoubleWithUnit sum(0, l1b[0]->altitude(i).units);
  BOOST_FOREACH(const boost::shared_ptr<Level1b>& f, l1b)
    sum += f->altitude(i);
  return sum / ((int) l1b.size());
}

SpectralRange Level1bAverage::radiance(int Spec_index) const
{
  SpectralRange t = l1b[0]->radiance(Spec_index);
  Array<double, 1> rad(t.data());
  Array<double, 1> sum(rad.shape());
  sum = 0;
  bool have_uncertainty = (t.uncertainty().rows() > 0);
  if(have_uncertainty)
    sum += sqr(t.uncertainty());
  for(int i = 1; i < (int) l1b.size(); ++i) {
    SpectralRange t2 = l1b[i]->radiance(Spec_index);
    rad += t2.data() * FullPhysics::conversion(t2.units(), t.units());
    if(have_uncertainty)
      sum += sqr(t2.uncertainty() * FullPhysics::conversion(t2.units(), t.units()));
  }
  rad /= ((int) l1b.size());
  Array<double, 1> uncer;
  if(have_uncertainty)
    uncer.reference(Array<double,1>(sqrt(sum) / (int) l1b.size()));
  return SpectralRange(rad, t.units(), uncer);
}

template <typename T>
bool Level1bAverage::check_field_equal(T && check_field) const {
    auto first_result = std::bind(check_field, l1b[0])();
    return std::all_of(l1b.begin() + 1, l1b.end(), [first_result, check_field](boost::shared_ptr<Level1b> l1b_i){return std::bind(check_field, l1b_i)() == first_result;});
}

template <typename T>
bool Level1bAverage::check_field_equal(T && check_field, int arg1) const {
    auto first_result = std::bind(check_field, l1b[0], arg1)();
    return std::all_of(l1b.begin() + 1, l1b.end(), [first_result, check_field, arg1](boost::shared_ptr<Level1b> l1b_i){return std::bind(check_field, l1b_i, arg1)() == first_result;});
}

/* TODO: change signature to accept string name of function to provide better error output? */
template <typename T>
void Level1bAverage::assert_field_equal(T && check_field) const {
    bool field_equal = check_field_equal(check_field);
    if(!field_equal) {
        Exception e;
        e << "All instances of checked field not equal.\n";
        throw e;
    }
}

/* TODO: change signature to accept string name of function to provide better error output? */
template <typename T>
void Level1bAverage::assert_field_equal(T && check_field, int arg1) const {
    bool field_equal = check_field_equal(check_field, arg1);
    if(!field_equal) {
        Exception e;
        e << "All instances of checked field not equal.\n";
        throw e;
    }
}


void Level1bAverage::print(std::ostream& Os) const
{
  Os << "Level1bAverage";
}
