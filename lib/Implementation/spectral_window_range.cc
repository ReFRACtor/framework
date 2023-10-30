#include "spectral_window_range.h"
#include "sample_grid.h"
#include "fp_exception.h"
#include "fp_serialize_support.h"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SpectralWindowRange::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpectralWindow)
    & FP_NVP_(range) & FP_NVP_(bad_sample_mask) & FP_NVP_(disp);
}

FP_IMPLEMENT(SpectralWindowRange);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
typedef void (SpectralWindowRange::*cfunc)(const std::vector<boost::shared_ptr<SampleGrid> >&);
typedef const ArrayWithUnit<double, 3>& (SpectralWindowRange::*a1)(void) const;
REGISTER_LUA_DERIVED_CLASS(SpectralWindowRange, SpectralWindow)
.def(luabind::constructor<const ArrayWithUnit<double, 3>&>())
.def(luabind::constructor<const ArrayWithUnit<double, 3>&, const Array<bool, 2>&>())
.def("range_array", ((a1) &SpectralWindowRange::range_array))
.def("dispersion", ((cfunc) &SpectralWindowRange::dispersion))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Construct a new window from the supplied information. We take
/// an array, which is number_spectrometer x number_microwindow x
/// 2. For a particular microwindow, we have 2 values, a lower and
/// upper range.
///
/// Note that we require the number of microwindows to be same for all
/// the spectrometers. This is just for convenience, it makes a
/// simpler interface. It is perfectly ok to have microwindow
/// with ranges of 0 to 0 - so you can simple set the number of
/// microwindows to whatever the maximum number is and use null range
/// microwindows for microwindows that aren't needed.
//-----------------------------------------------------------------------

SpectralWindowRange::SpectralWindowRange
(const ArrayWithUnit<double, 3>& Microwindow_ranges)
    : range_(Microwindow_ranges)
{
    if(range_.value.depth() != 2) {
        Exception e;
        e << "Microwindow_ranges needs to have a depth of 2\n"
          << "(for lower and upper bound). Found depth of " << range_.value.depth();
        throw e;
    }
}

//-----------------------------------------------------------------------
/// In addition to constructing the object using the microwindow
/// ranges, adds a bad sample mask argument.
///
/// The bad sample mask is sized num_bands x num_samples. A value of
/// true in the array means the sample is marked as bad.
//-----------------------------------------------------------------------

SpectralWindowRange::SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges,
                                         const Array<bool, 2>& Bad_sample_mask)
    : range_(Microwindow_ranges)
{
    if(range_.value.depth() != 2) {
        Exception err;
        err << "Microwindow_ranges needs to have a depth of 2\n"
            << "(for lower and upper bound). Found depth of " << range_.value.depth();
        throw err;
    }

    if(Bad_sample_mask.rows() != range_.rows()) {
        Exception err;
        err << "Number of channels in bad sample mask: " << Bad_sample_mask.rows()
            << " must match the number in the spectral window ranges: " << range_.rows();
        throw err;
    }

    for(int chan_idx = 0; chan_idx < Bad_sample_mask.rows(); chan_idx++) {
        bad_sample_mask_.push_back( Bad_sample_mask(chan_idx, Range::all()) );
    }
}

//-----------------------------------------------------------------------
/// In addition to constructing the object using the microwindow
/// ranges, adds a bad sample mask argument.
///
/// The bad sample mask is a vector sized for the number of instrument
/// channels with each array having a true/false for each sample posistion.
/// A value of true in the array means the sample is marked as bad.
//-----------------------------------------------------------------------

SpectralWindowRange::SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges,
                                         const std::vector<Array<bool, 1> >& Bad_sample_mask)
    : range_(Microwindow_ranges), bad_sample_mask_(Bad_sample_mask)
{
    if(range_.value.depth() != 2) {
        Exception err;
        err << "Microwindow_ranges needs to have a depth of 2\n"
            << "(for lower and upper bound). Found depth of " << range_.value.depth();
        throw err;
    }

    if((int) Bad_sample_mask.size() != range_.rows()) {
        Exception err;
        err << "Number of channels in bad sample mask: " << Bad_sample_mask.size()
            << " must match the number in the spectral window ranges: " << range_.rows();
        throw err;
    }
}


// See base class for a description.
SpectralBound SpectralWindowRange::spectral_bound() const
{
    Range ra = Range::all();
    std::vector<DoubleWithUnit> lbound, ubound;

    for(int i = 0; i < number_spectrometer(); ++i) {
        DoubleWithUnit lv(min(range_.value(i, ra, ra)), range_.units);
        DoubleWithUnit uv(max(range_.value(i, ra, ra)), range_.units);

        if(range_.units.is_commensurate(units::sample_index) && disp_.size() > 0) {
            SpectralDomain pgrid = disp_[i]->pixel_grid();

            // Special handling for SampleIndex, so it isn't out of range
            if(lv.units.is_commensurate(units::sample_index)) {
                if(lv.value < 1) {
                    lv.value = 1;
                }

                if(lv.value > pgrid.data().rows()) {
                    lv.value = pgrid.data().rows();
                }

                if(uv.value < 1) {
                    uv.value = 1;
                }

                if(uv.value > pgrid.data().rows()) {
                    uv.value = pgrid.data().rows();
                }
            }

            lv = lv.convert_wave(pgrid.units(), pgrid);
            uv = uv.convert_wave(pgrid.units(), pgrid);

            // Might switch the direction of upper and lower
            if(lv.value > uv.value) {
                std::swap(lv.value, uv.value);
            }
        }

        lbound.push_back(lv);
        ubound.push_back(uv);
    }

    return SpectralBound(lbound, ubound);
}

// See base class for a description.

std::vector<int>
SpectralWindowRange::grid_indexes(const SpectralDomain& Grid, int sensor_index) const
{
    range_check(sensor_index, 0, number_spectrometer());
    std::vector<int> res;
    Array<double, 1> g;

    if(range_.units.is_commensurate(units::sample_index)) {
        g.resize(Grid.sample_index().shape());
        g = blitz::cast<double>(Grid.sample_index());
    }
    else if(range_.units.is_commensurate(units::inv_cm)) {
        g.reference(Grid.wavenumber(range_.units));
    }
    else {
        g.reference(Grid.wavelength(range_.units));
    }
    if(bad_sample_mask_.size() != 0 && g.rows() > bad_sample_mask(sensor_index).rows()) {
      Exception err;
      err << "Bad sample mask and input SpectralDomain have a size mismatch.\n"
	  << "  sensor index:         " << sensor_index << "\n"
	  << "  SpectralDomain size:  " << g.rows() << "\n"
	  << "  Bad sample mask size: " << bad_sample_mask(sensor_index).rows() << "\n";
      throw err;
    }

    for(int samp_idx = 0; samp_idx < g.rows(); ++samp_idx) {
        bool ok = false;

        // Check if sample is bad first, if it is not then check if
        // it falls within the spectral ranges
        if(bad_sample_mask_.size() == 0 || !bad_sample_mask_[sensor_index](samp_idx)) {
            for(int win_idx = 0; win_idx < range_.value.cols(); ++win_idx) {
                if( g(samp_idx) >= range_.value(sensor_index, win_idx, 0) &&
                    g(samp_idx) <= range_.value(sensor_index, win_idx, 1)) {

                    ok = true;
                }
            }
        }

        if(ok) {
            res.push_back(samp_idx);
        }
    }

    return res;
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void SpectralWindowRange::print(std::ostream& Os) const
{
    Os << "SpectralWindowRange:\n"
       << "   Units: " << range_.units;

    for(int i = 0; i < range_.value.rows(); ++i) {
        Os << "  Spec[" << i << "]:\n";

        for(int j = 0; j < range_.value.cols(); ++j)
            Os << "    Microwindow[" << j << "]: ("
               << range_.value(i, j, 0) << ", " << range_.value(i, j, 1) << ")\n";
    }
}
