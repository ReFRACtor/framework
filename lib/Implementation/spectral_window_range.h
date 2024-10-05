#ifndef SPECTRAL_WINDOW_RANGE_H
#define SPECTRAL_WINDOW_RANGE_H

#include "spectral_window.h"
#include "array_with_unit.h"
#include "sample_grid.h"

namespace FullPhysics {

/****************************************************************//**
  This is an implementation of a SpectralWindow that covers a fixed
  window. The window can be made up of multiple microwindow ranges if
  desired. It can also include a bad pixel mask.

  Note that there are a few closely related classes, with similar
  sounding names. See \ref spectrumdoxygen for a description of each
  of these.
*******************************************************************/

class SpectralWindowRange : public SpectralWindow {
public:
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges);

  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges,
		      const blitz::Array<bool, 2>& Bad_sample_mask);

  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges,
		      const std::vector<blitz::Array<bool, 1> >& Bad_sample_mask);

  virtual ~SpectralWindowRange() {}

  //-----------------------------------------------------------------------
  /// The Dispersion is optional. If supplied, we can convert from
  /// sample_index to wavelength/wavenumber.
  //-----------------------------------------------------------------------

  const std::vector<boost::shared_ptr<SampleGrid> >& dispersion() const
  {
    return disp_;
  }
  void dispersion(const std::vector<boost::shared_ptr<SampleGrid> >& D)
  {
    disp_ = D;
  }

  virtual int number_spectrometer() const
  {
    return range_.value.rows();
  }

  virtual SpectralBound spectral_bound() const;

  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, int sensor_index) const;

  void print(std::ostream& Os) const;

  const ArrayWithUnit<double, 3>& range_array() const
  {
    return range_;
  }

  void range_array(const ArrayWithUnit<double, 3>& Ran)
  {
    range_ = Ran;
  }

  const blitz::Array<bool, 1>& bad_sample_mask(int sensor_index) const
  {
    if (sensor_index < 0 || sensor_index >= (int) bad_sample_mask_.size()) {
      Exception err;
      err << "Bad sensor index " << sensor_index << ", out of range of bad sample mask size: "
	  << bad_sample_mask_.size();
      throw err;
    }

    return bad_sample_mask_[sensor_index];
  }

  void bad_sample_mask(const blitz::Array<bool, 1>& M, int sensor_index)
  {
    // Populate with empty arrays for all channels
    if (bad_sample_mask_.size() == 0) {
      if (sensor_index < 0 || sensor_index >= (int) range_.rows()) {
	Exception err;
	err << "Bad sensor index " << sensor_index << ", out of range of bad sample mask size: "
	    << range_.rows();
	throw err;
      }
      for(int chan_idx = 0; chan_idx < range_.rows(); chan_idx++) {
	bad_sample_mask_.push_back( blitz::Array<bool, 1>() );
      }
    }

    if (sensor_index < 0 || sensor_index >= (int) bad_sample_mask_.size()) {
      Exception err;
      err << "Bad sensor index " << sensor_index << ", out of range of bad sample mask size: "
	  << bad_sample_mask_.size();
      throw err;
    }
    bad_sample_mask_[sensor_index].reference(M.copy());
  }

private:
  ArrayWithUnit<double, 3> range_;
    
  // Mask of bad samples, True for a bad sample, False for a good one
  // One array per channel
  std::vector<blitz::Array<bool, 1> > bad_sample_mask_;

  std::vector<boost::shared_ptr<SampleGrid> > disp_;
  SpectralWindowRange() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(SpectralWindowRange);
#endif
