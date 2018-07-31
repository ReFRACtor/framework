#ifndef LEVEL_1B_AVERAGE_H
#define LEVEL_1B_AVERAGE_H
#include "level_1b.h"
#include <boost/shared_ptr.hpp>
#include <functional>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This reads averages a set of Level1b classes to get the various
  values. This is used for example on Gosat, where we average the S
  and P data.
*******************************************************************/
class Level1bAverage: public Level1b {
public:
  Level1bAverage(const std::vector<boost::shared_ptr<Level1b> >& Data)
    : l1b(Data) {}
  virtual ~Level1bAverage() {}
  virtual int number_spectrometer() const;
  virtual DoubleWithUnit relative_velocity(int Spec_index) const;
  virtual Time time(int Spec_index) const;
  virtual SpectralDomain sample_grid(int Spec_index) const;

  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual DoubleWithUnit altitude(int i) const;

  virtual void print(std::ostream& Os) const;
  virtual SpectralRange radiance(int Spec_index) const;
private:
  template <typename T>
  bool check_field_equal(T && check_field) const;
  template <typename T>
  bool check_field_equal(T && check_field, int arg1) const;
  template <typename T>
  void assert_field_equal(T && check_field) const;
  template <typename T>
  void assert_field_equal(T && check_field, int arg1) const;
  std::vector<boost::shared_ptr<Level1b> > l1b;
};
}
#endif
