#include "absco_aer.h"
#include "fp_exception.h"
#include <iomanip>
using namespace FullPhysics;
using namespace blitz;
#ifdef HAVE_LUA
#include "register_lua.h"
typedef void (AbscoAer::*a1)(AbscoAer::InterpolationType);

namespace FullPhysics {
void register_lua_AbscoAer(lua_State *ls) { \
luabind::module(ls) [ 
luabind::class_<AbscoAer,Absco,boost::shared_ptr<GasAbsorption> >
("AbscoAer")
.def(luabind::constructor<const std::string&>())
.def(luabind::constructor<const std::string&,double>())
.def(luabind::constructor<const std::string&,const SpectralBound&, 
     const std::vector<double>&>())
.def("interpolation_type", ((a1) & AbscoAer::interpolation_type))
];}
}
#endif

//-----------------------------------------------------------------------
/// Read the given Absco file. You can optionally set the number of
/// lines of data to cache, and a scaling to apply to the underlying
/// table. 
//-----------------------------------------------------------------------

AbscoAer::AbscoAer(const std::string& Fname, double Table_scale,
                   int Cache_nline, InterpolationType Itype)
  : itype_(Itype)
{
  load_file(Fname, Table_scale, Cache_nline);
}

//-----------------------------------------------------------------------
/// Read the given Absco file. You can optionally set the number of
/// lines of data to cache, and a scaling to apply to the underlying
/// table. 
//-----------------------------------------------------------------------

AbscoAer::AbscoAer(const std::string& Fname, 
                   const SpectralBound& Spectral_bound,
                   const std::vector<double>& Table_scale,
                   int Cache_nline, InterpolationType Itype)
  : itype_(Itype)
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

void AbscoAer::load_file(const std::string& Fname)
{
  load_file(Fname, sb, table_scale_, cache_nline);
}

//-----------------------------------------------------------------------
/// Read the given Absco file. Note that you don't normally call this
/// directly, but rather use the constructor. However there are
/// certain uses where changing the absco table can be useful, such
/// as the python code l2_fp_fm_errors.py.
//-----------------------------------------------------------------------

void AbscoAer::load_file(const std::string& Fname, double Table_scale, 
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

void AbscoAer::load_file(const std::string& Fname, 
                         const SpectralBound& Spectral_bound,
                         const std::vector<double>& Table_scale,
                         int Cache_nline)
{
  sb = Spectral_bound;
  std::vector<double> tcopy(Table_scale);
  table_scale_.swap(tcopy);
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
  cache_nline = Cache_nline;
  cache_double_lbound = 0;
  cache_double_ubound = 0;
  cache_float_lbound = 0;
  cache_float_ubound = 0;

  // Reset caches
  read_cache_float.resize(0,0,0,0,0);
  read_cache_double.resize(0,0,0,0,0);

  hfile.reset(new HdfFile(Fname));

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

  // This may have a "\0" in it, so we create a string twice, the
  // second one ends at the first "\0".
  field_name = "Cross_Section";
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
  // The VMR in the AbscoAer file is in ppmv, while we expect just a
  // straight vmr. So multiple 1e-6 to get to the write units.
  for(int i = 0; i < (int) bvmr.size(); ++i)
    bvmr[i] *= 1e-6;
}

//-----------------------------------------------------------------------
/// Find the array locations where the wavenumber is contained If no
/// extents are defined, this just returns beginning and end of whole
/// absco array. (version friendlier to swig/python)
//-----------------------------------------------------------------------

void AbscoAer::wn_extent(double Wn_in, double& X, double& Y) const
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

const std::pair<double*, double*> AbscoAer::wn_extent(double Wn_in) const
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

double AbscoAer::table_scale(double wn) const 
{
  if(sb.number_spectrometer() > 0) {
    int index = sb.spectral_index(DoubleWithUnit(wn, "cm^-1"));
    if(index < 0)
      return 0.0;
    return table_scale_[index];
  }
  return table_scale_[0];
}

// See base class for description
bool AbscoAer::have_data(double wn) const
{
  try {
    // If wn_extent does not exception, then we know we have
    // a wavenumber that falls with in a range, no need to check again
    // now just check that the table scale for a wn is non zero
    auto UNUSED(extents) = wn_extent(wn);
    return table_scale(wn) > 0.0;
  } catch (Exception e) {
    // If there is an exception that means wn_extent did 
    // not find a valid wavenumber
    return false;
  }
}

// Calculate the wn index number of the data.
int AbscoAer::wn_index(double Wn_in) const
{
  // This will throw an exception if it can not find the wavenumber's locations
  auto extent = wn_extent(Wn_in);

  double *wnptr = std::lower_bound(extent.first, extent.second, Wn_in);
  double f;
  if(wnptr == extent.first)
    f = 1.0;
  else
    f = (Wn_in - *(wnptr - 1)) / (*wnptr - *(wnptr - 1));
  if(itype_ == THROW_ERROR_IF_NOT_ON_WN_GRID && f > 0.1 && f < 0.9) {
    Exception e;
    e << std::setprecision(8)
      << "AbscoAER does not interpolate in wavenumber direction.\n"
      << "The interpolation doesn't work well near peaks, so we\n"
      << "just don't do it. The requested Wn needs to be within 10%\n"
      << "of the absco wn grid, return the value for the closest number.\n"
      << "Wn:        " << Wn_in << "\n"
      << "Wnptr - 1: " << *(wnptr - 1) << "\n"
      << "Wnptr:     " << *wnptr << "\n"
      << "Frac:      " << f << " (should be < 0.1 or > 0.9)\n"
      << "AbscoAer:\n"
      << *this << "\n";
    throw e;
  }
  if(f < 0.5)                        // First point is closest.
    --wnptr;
  return (int) (wnptr - wnfront);
}

// Calculate the wn index number of the data.
int AbscoAer::wn_index(double Wn_in, double& F) const
{
  // This will throw an exception if it can not find the wavenumber's locations
  auto extent = wn_extent(Wn_in);

  double *wnptr = std::lower_bound(extent.first, extent.second, Wn_in);
  if(wnptr == extent.first) {
    F = 0.0;
  } else {
    F = (Wn_in - *(wnptr - 1)) / (*wnptr - *(wnptr - 1));
    --wnptr;
  }
  return (int) (wnptr - wnfront);
}

// See base class for description
Array<double, 3> AbscoAer::read_double(double Wn_in) const
{
  int wi;
  double f;
  if(itype_ == INTERPOLATE_WN)
    wi = wn_index(Wn_in, f);
  else
    wi = wn_index(Wn_in);
  if(wi < cache_double_lbound ||
     wi >= cache_double_ubound)
    swap<double>(wi);
  // The Aer data is in temperature x pressure x broadner order. The
  // base class expects this to be pressure x temperature x broadner,
  // so we transpose this.
  if(itype_ != INTERPOLATE_WN) 
    return read_cache<double>()(wi - cache_double_lbound, Range::all(),
	Range::all(), Range::all(), 0).transpose(secondDim, firstDim, thirdDim);
  Array<double, 3> r1;
  if(wi+1 < cache_double_lbound ||
     wi+1 >= cache_double_ubound)
    r1.reference(read_cache<double>()(wi - cache_double_lbound, Range::all(),
      Range::all(), Range::all(), 0).transpose(secondDim, firstDim, thirdDim));
  else {
    r1.reference(read_cache<double>()(wi - cache_double_lbound, Range::all(),
      Range::all(), Range::all(), 0).transpose(secondDim, firstDim,
					    thirdDim).copy());
    swap<double>(wi+1);
  }
  Array<double, 3> r2 = read_cache<double>()(wi+1 - cache_double_lbound,
      Range::all(),
      Range::all(), Range::all(), 0).transpose(secondDim, firstDim,
					    thirdDim);
  Array<double, 3> res(r2.shape());
  res = r1 * (1 - f) + r2 * f;
  return res;
}

// See base class for description
Array<float, 3> AbscoAer::read_float(double Wn_in) const
{
  int wi;
  double f;
  if(itype_ == INTERPOLATE_WN)
    wi = wn_index(Wn_in, f);
  else
    wi = wn_index(Wn_in);
  if(wi < cache_float_lbound ||
     wi >= cache_float_ubound)
    swap<float>(wi);
  // The Aer data is in temperature x pressure x broadner order. The
  // base class expects this to be pressure x temperature x broadner,
  // so we transpose this.
  if(itype_ != INTERPOLATE_WN) 
    return read_cache<float>()(wi - cache_float_lbound, Range::all(),
       Range::all(), Range::all(), 0).transpose(secondDim, firstDim, thirdDim);
  Array<float, 3> r1;
  if(wi+1 < cache_float_lbound ||
     wi+1 >= cache_float_ubound)
    r1.reference(read_cache<float>()(wi - cache_float_lbound, Range::all(),
     Range::all(), Range::all(), 0).transpose(secondDim, firstDim, thirdDim));
  else {
    r1.reference(read_cache<float>()(wi - cache_float_lbound, Range::all(),
    Range::all(), Range::all(), 0).transpose(secondDim, firstDim,
					    thirdDim).copy());
    swap<float>(wi+1);
  }
  Array<float, 3> r2 = read_cache<float>()(wi+1 - cache_float_lbound,
      Range::all(),
      Range::all(), Range::all(), 0).transpose(secondDim, firstDim,
					    thirdDim);
  Array<float, 3> res(r2.shape());
  res = r1 * (1 - (float) f) + r2 * ((float) f);
  return res;
}

// See base class for description
Array<double, 4> AbscoAer::read_double_2b(double Wn_in) const
{
  int wi;
  double f;
  if(itype_ == INTERPOLATE_WN)
    wi = wn_index(Wn_in, f);
  else
    wi = wn_index(Wn_in);
  if(wi < cache_double_lbound ||
     wi >= cache_double_ubound)
    swap<double>(wi);
  // The Aer data is in temperature x pressure x broadner order. The
  // base class expects this to be pressure x temperature x broadner,
  // so we transpose this.
  if(itype_ != INTERPOLATE_WN) 
    return read_cache<double>()(wi - cache_double_lbound, Range::all(),
	Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim, thirdDim, fourthDim);
  Array<double, 4> r1;
  if(wi+1 < cache_double_lbound ||
     wi+1 >= cache_double_ubound)
    r1.reference(read_cache<double>()(wi - cache_double_lbound, Range::all(),
       Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim, thirdDim, fourthDim));
  else {
    r1.reference(read_cache<double>()(wi - cache_double_lbound, Range::all(),
	      Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim,
					  thirdDim, fourthDim).copy());
    swap<double>(wi+1);
  }
  Array<double, 4> r2 = read_cache<double>()(wi+1 - cache_double_lbound,
      Range::all(),
      Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim,
												 thirdDim, fourthDim);
  Array<double, 4> res(r2.shape());
  res = r1 * (1 - f) + r2 * f;
  return res;
}

// See base class for description
Array<float, 4> AbscoAer::read_float_2b(double Wn_in) const
{
  int wi;
  double f;
  if(itype_ == INTERPOLATE_WN)
    wi = wn_index(Wn_in, f);
  else
    wi = wn_index(Wn_in);
  if(wi < cache_double_lbound ||
     wi >= cache_double_ubound)
    swap<float>(wi);
  // The Aer data is in temperature x pressure x broadner order. The
  // base class expects this to be pressure x temperature x broadner,
  // so we transpose this.
  if(itype_ != INTERPOLATE_WN) 
    return read_cache<float>()(wi - cache_float_lbound, Range::all(),
	Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim, thirdDim, fourthDim);
  Array<float, 4> r1;
  if(wi+1 < cache_float_lbound ||
     wi+1 >= cache_float_ubound)
    r1.reference(read_cache<float>()(wi - cache_float_lbound, Range::all(),
       Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim, thirdDim, fourthDim));
  else {
    r1.reference(read_cache<float>()(wi - cache_float_lbound, Range::all(),
	      Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim,
					  thirdDim, fourthDim).copy());
    swap<float>(wi+1);
  }
  Array<float, 4> r2 = read_cache<float>()(wi+1 - cache_float_lbound,
      Range::all(),
      Range::all(), Range::all(), Range::all()).transpose(secondDim, firstDim,
												 thirdDim, fourthDim);
  Array<float, 4> res(r2.shape());
  res = r1 * (1 - ((float) f)) + r2 * ((float) f);
  return res;
}

//-----------------------------------------------------------------------
/// Make sure the row i is found in the read_cache, possibly reading
/// data if it isn't found.
//-----------------------------------------------------------------------

template<class T> void AbscoAer::swap(int i) const
{
  // First time through, set up space for cache.
 if(read_cache<T>().extent(firstDim) == 0)
   read_cache<T>().resize(cache_nline, tgrid.cols(), tgrid.rows(),
			  (number_broadener() >= 1 ? number_broadener_vmr(0) : 1),
			  (number_broadener() >= 2 ? number_broadener_vmr(1) : 1));
  int nl = read_cache<T>().extent(firstDim);
  // Either read 3d, 4d or 5d data. We tell which kind by the number
  // of broadners.
  int st0 = (i / nl) * nl;
  int sz0 = std::min(nl, wngrid.rows() - st0);
  bound_set<T>(st0, sz0);
  if(number_broadener() > 2 || number_broadener() < 0)
    throw Exception("Out of range number_broadener(). This shouldn't be possible");
  if(number_broadener() == 2) {
    TinyVector<int, 5> start, size;
    start = st0, 0, 0, 0, 0;
    size = sz0,
      read_cache<T>().shape()[1],
      read_cache<T>().shape()[2],
      read_cache<T>().shape()[3],
      read_cache<T>().shape()[4];
    read_cache<T>()(Range(0, size(0) - 1), Range::all(),
		    Range::all(), Range::all(), Range::all()) =
      hfile->read_field<T, 5>(field_name, start, size);
  }
  if(number_broadener() == 1) {
    // 4d case
    TinyVector<int, 4> start, size;
    start = st0, 0, 0, 0;
    size = sz0,
      read_cache<T>().shape()[1],
      read_cache<T>().shape()[2],
      read_cache<T>().shape()[3];
    read_cache<T>()(Range(0, size(0) - 1), Range::all(),
		    Range::all(), Range::all(), 0) =
      hfile->read_field<T, 4>(field_name, start, size);
  }
  if(number_broadener() == 0) {
    // 3d case
    TinyVector<int, 3> start, size;
    start = st0, 0, 0;
    size = sz0,
      read_cache<T>().shape()[1],
      read_cache<T>().shape()[2];
    read_cache<T>()(Range(0, size(0) - 1), Range::all(), Range::all(), 0, 0) =
      hfile->read_field<T, 3>(field_name, start, size);
  }
}

void AbscoAer::print(std::ostream& Os) const
{
  Os << "AbscoAer" << "\n"
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
