#include "absco_coeff.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include <iomanip>
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbscoCoeff::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Absco)
    & FP_NVP(cache_nline) & FP_NVP(hfile)
    & FP_NVP(sb) & FP_NVP_(table_scale);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void AbscoCoeff::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void AbscoCoeff::load(Archive & UNUSED(ar),
		    const unsigned int UNUSED(version))
{
  load_file();
}

FP_IMPLEMENT(AbscoCoeff);
#endif

//-----------------------------------------------------------------------
/// Read the given Absco file. You can optionally set the number of
/// lines of data to cache, and a scaling to apply to the underlying
/// table. 
//-----------------------------------------------------------------------

AbscoCoeff::AbscoCoeff(const std::string& Fname, double Table_scale,
                   int Cache_nline)
{
  load_file(Fname, Table_scale, Cache_nline);
}

//-----------------------------------------------------------------------
/// Read the given Absco file. You can optionally set the number of
/// lines of data to cache, and a scaling to apply to the underlying
/// table. 
//-----------------------------------------------------------------------

AbscoCoeff::AbscoCoeff(const std::string& Fname, 
		       const SpectralBound& Spectral_bound,
		       const std::vector<double>& Table_scale,
		       int Cache_nline)
{
  load_file(Fname, Spectral_bound, Table_scale, Cache_nline);
}

//-----------------------------------------------------------------------
/// Load file, using the same table scaling as we already have.  Note
/// that you don't normally call this directly, but rather use the
/// constructor. However there are certain uses where changing the
/// absco table can be useful, such as the python code
/// l2_fp_fm_errors.py.
//-----------------------------------------------------------------------

void AbscoCoeff::load_file(const std::string& Fname)
{
  load_file(Fname, sb, table_scale_, cache_nline);
}

//-----------------------------------------------------------------------
/// Read the given Absco file. Note that you don't normally call this
/// directly, but rather use the constructor. However there are
/// certain uses where changing the absco table can be useful, such
/// as the python code l2_fp_fm_errors.py.
//-----------------------------------------------------------------------

void AbscoCoeff::load_file(const std::string& Fname, double Table_scale, 
                         int Cache_nline)
{
  SpectralBound empty;
  std::vector<double> tscale;
  tscale.push_back(Table_scale);
  load_file(Fname, empty, tscale, Cache_nline);
}

//-----------------------------------------------------------------------
/// Read the given Absco file. Note that you don't normally call this
/// directly, but rather use the constructor. However there are
/// certain uses where changing the absco table can be useful, such
/// as the python code l2_fp_fm_errors.py.
//-----------------------------------------------------------------------

void AbscoCoeff::load_file(const std::string& Fname, 
                         const SpectralBound& Spectral_bound,
                         const std::vector<double>& Table_scale,
                         int Cache_nline)
{
  sb = Spectral_bound;
  std::vector<double> tcopy(Table_scale);
  table_scale_.swap(tcopy);
  cache_nline = Cache_nline;
  hfile = boost::make_shared<HdfFile>(Fname);
  load_file();
}

//-----------------------------------------------------------------------
/// Internal reading of file, once we have set up our internal variables.
//-----------------------------------------------------------------------
 
void AbscoCoeff::load_file()
{  
  if(sb.number_spectrometer() > 0 && 
     (int) table_scale_.size() != sb.number_spectrometer()) {
    Exception e;
    e << "Table_scale size needs to match the number of spectrometers\n"
      << "  Table_scale size: " << table_scale_.size() << "\n"
      << "  Number spectrometer: " << sb.number_spectrometer();
    throw e;
  }
  if(sb.number_spectrometer() ==0 && (int) table_scale_.size() != 1) {
    Exception e;
    e << "Table_scale size needs to be exactly 1\n"
      << "  Table_scale size: " << table_scale_.size() << "\n";
    throw e;
  }
  cache_double_lbound = 0;
  cache_double_ubound = 0;
  cache_float_lbound = 0;
  cache_float_ubound = 0;

  // Reset caches
  read_cache_float.reference(blitz::Array<float, 4>());
  read_cache_double.reference(blitz::Array<double, 4>());

  // Read optional metadata
  if (hfile->has_object("Extent_Ranges")) {
      extent_range.reference(hfile->read_field<double, 2>("Extent_Ranges"));
      extent_index.reference(hfile->read_field<int, 2>("Extent_Indices"));
  }

  // Read data fields
  ArrayWithUnit<double, 2> pgrid_unit =
    hfile->read_field_with_unit<double, 2>("P_layer").convert(Unit("Pa"));
  // Note that P_layer only has data where the temperature is also
  // available. So go through each layer and find the first P_layer
  // that is finite. We then save that value.
  pgrid.resize(pgrid_unit.rows());
  for(int i = 0; i < pgrid_unit.rows(); ++i) {
    bool have_value = false;
    for(int j = 0; j < pgrid_unit.cols(); ++j) {
      if(!have_value && isfinite(pgrid_unit.value(i,j))) {
	pgrid(i) = pgrid_unit.value(i,j);
	have_value = true;
      }
    }
    if(!have_value) {
      Exception e;
      e << "No pressure values found for layer " << i;
      throw e;
    }
  }
  ArrayWithUnit<double, 2> tgrid_unit =
    hfile->read_field_with_unit<double, 2>("T_layer").convert(Unit("K"));
  tgrid.reference(tgrid_unit.value);
  ArrayWithUnit<double, 1> sg =
    hfile->read_field_with_unit<double, 1>("Spectral_Grid").convert_wave(Unit("cm^-1"));
  wngrid.reference(sg.value);
  wnfront = &wngrid(0);
  cross_sec_coeff.reference(hfile.read_field<double, 3>("Cross_Section_coeffs"));

  // This may have a "\0" in it, so we create a string twice, the
  // second one ends at the first "\0".
  field_name = "Cross_Section_repvecs";
  // Determine if data is float or double. We use this to optimize the
  // read. 
  H5::DataSet d = hfile->h5_file().openDataSet(field_name);
  if(d.getDataType().getSize() == 4)
    is_float_ = true;
  else
    is_float_ = false;
  if(hfile->has_object("H2O_VMR")) {
    bname.push_back("h2o");
    bvmr.push_back(hfile->read_field<double, 1>("H2O_VMR").copy());
  } if(hfile->has_object("O2_VMR")) {
    bname.push_back("o2");
    bvmr.push_back(hfile->read_field<double, 1>("O2_VMR").copy());
  }
}


//-----------------------------------------------------------------------
/// Find the array locations where the wavenumber is contained If no
/// extents are defined, this just returns beginning and end of whole
/// absco array. (version friendlier to swig/python)
//-----------------------------------------------------------------------

void AbscoCoeff::wn_extent(double Wn_in, double& X, double& Y) const
{
  std::pair<double*, double*> t = wn_extent(Wn_in);
  X = *t.first;
  Y = *t.second;
}

//-----------------------------------------------------------------------
/// Find the array locations where the wavenumber is contained If no
/// extents are defined, this just returns beginning and end of whole
/// absco array.
//-----------------------------------------------------------------------

const std::pair<double*, double*> AbscoCoeff::wn_extent(double Wn_in) const
{
    if (extent_range.rows() > 0) {
        // Find first valid range and return
        for (int range_idx = 0; range_idx < extent_range.rows(); range_idx++) {
            if (extent_range(range_idx, 0) <= Wn_in && Wn_in <= extent_range(range_idx, 1)) {
                int beg_idx = extent_index(range_idx, 0);
                int end_idx = extent_index(range_idx, 1);
                // Remove const to satisfy std::pair constructor
                double* beg_loc = const_cast<double*>(&wngrid(beg_idx));
                double* end_loc = const_cast<double*>(&wngrid(end_idx));
                return std::pair<double*, double*>(beg_loc, end_loc);
            }
        }

        // If nothing matches in loop above then the wave number was not
        // found within one of the ranges
        Exception err;
        err << "The ABSCO file " << hfile->file_name() 
            << " does not have an extent that contains the wavenumber: "
            << std::setprecision(8) << Wn_in;
        throw err;
    } else {
        double* wnback = const_cast<double *>(&wngrid(wngrid.rows() - 1));

        if (Wn_in >= *wnfront && Wn_in <= *wnback) {
            return std::pair<double*, double*>(wnfront, wnback);
        } else {
            Exception err;
            err << "The ASBCO file " << hfile->file_name() << "\n"
                << " does not have data for wavenumber: " << Wn_in;
            throw err;
        }
    }
}

double AbscoCoeff::table_scale(double wn) const 
{
  if(sb.number_spectrometer() > 0) {
    int index = sb.spectral_index(DoubleWithUnit(wn, "cm^-1"));
    if(index < 0)
      return 1.0;
    return table_scale_[index];
  }
  return table_scale_[0];
}

// See base class for description
bool AbscoCoeff::have_data(double wn) const
{
    try {
      // If wn_extent does not exception, then we know we have
      // a wavenumber that falls with in a range, no need to check again
      // now just check that the table scale for a wn is non zero
      auto UNUSED(extents) = wn_extent(wn);
      return table_scale(wn) > 0.0;
    } catch (const Exception& e) {
      // If there is an exception that means wn_extent did 
      // not find a valid wavenumber
      return false;
    }
}

// Calculate the wn index number of the data.
void AbscoCoeff::wn_index(double Wn_in, int& Wn_index, double & F) const
{
  // This will throw an exception if it can not find the wavenumber's locations
  auto extent = wn_extent(Wn_in);

  double *wnptr = std::lower_bound(extent.first, extent.second, Wn_in);
  if(wnptr == extent.first) {
    F = 0.0;
    Wn_index = 0;
  } else {
    F = (Wn_in - *(wnptr - 1)) / (*wnptr - *(wnptr - 1));
    --wnptr;
    Wn_index = (int) (wnptr - wnfront);
  }
}

// See base class for description
Array<double, 3> AbscoCoeff::read_double(double Wn_in) const
{
  double f;
  int wi;
  wn_index(Wn_in, wi, f);
  if(wi < cache_double_lbound ||
     wi >= cache_double_ubound ||
     ((wi+1) < wngrid.rows() &&
      (wi+1) >= cache_double_ubound))
    swap<float>(wi);
  return read_cache<double>()(Range::all(), Range::all(), Range::all(), 
                              wi - cache_double_lbound);
}

// See base class for description
Array<float, 3> AbscoCoeff::read_float(double Wn_in) const
{
  // Temp, we'll skip float. The test data from Matt is double, and
  // read_float is just read_double with a different type. We can
  // fill this in once read_double is debugged fully
  throw Exception("Not implemented yet");
  double f;
  int wi;
  wn_index(Wn_in, wi, f);
  if(wi < cache_float_lbound ||
     wi >= cache_float_ubound ||
     ((wi+1) < wngrid.rows() &&
      (wi+1) >= cache_float_ubound))
    swap<float>(wi);
  return read_cache<float>()(Range::all(), Range::all(), Range::all(), 
                             wi - cache_float_lbound);
}

//-----------------------------------------------------------------------
/// Make sure the row i is found in the read_cache, possibly reading
/// data if it isn't found.
//-----------------------------------------------------------------------

template<class T> void AbscoCoeff::swap(int i) const
{
  // First time through, set up space for cache. We add an extra point
  // at the end so we have the data needed for interpolating in the
  // frequency range.
  if(read_cache<T>().cols() == 0)
    read_cache<T>().resize(cross_sec_coeff.rows(), cache_nline + 1);
  int nl = read_cache<T>().cols();
  TinyVector<int, 2> start, size;
  start = 0, (i / cache_nline) * cache_nline;
  size = cross_sec_coeff.rows(), std::min(nl, wngrid.rows() - start(1));
  bound_set<T>((i / cache_nline) * cache_nline, size(1));
  read_cache<T>()(Range::all(), Range(0, size(1) - 1)) =
    hfile->read_field<T, 2>(field_name, start, size);
}

void AbscoCoeff::print(std::ostream& Os) const
{
  Os << "AbscoCoeff" << "\n"
     << "  File name:    " << hfile->file_name() << "\n";
  if(sb.number_spectrometer() ==0)
    Os << "  Scale factor: " << table_scale_[0] << "\n";
  else {
    Os << "  Scale factor: [";
    for(int i = 0; i < (int) table_scale_.size(); ++i) {
      Os << table_scale_[i];
      if(i < (int) table_scale_.size() - 1)
        Os << ", ";
    }
    Os << "]\n";
  }
}
