#include "absco_aer.h"
#include "fp_exception.h"
#include <iomanip>
using namespace FullPhysics;
using namespace blitz;
#ifdef HAVE_LUA
#include "register_lua.h"
namespace FullPhysics {
void register_lua_AbscoAer(lua_State *ls) { \
luabind::module(ls) [ 
luabind::class_<AbscoAer,Absco,boost::shared_ptr<GasAbsorption> >
("AbscoAer")
.def(luabind::constructor<const std::string&>())
.def(luabind::constructor<const std::string&,double>())
.def(luabind::constructor<const std::string&,const SpectralBound&, 
     const std::vector<double>&>())
];}
}
#endif

//-----------------------------------------------------------------------
/// Read the given Absco file. You can optionally set the number of
/// lines of data to cache, and a scaling to apply to the underlying
/// table. 
//-----------------------------------------------------------------------

AbscoAer::AbscoAer(const std::string& Fname, double Table_scale,
                   int Cache_nline)
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
  read_cache_float.resize(0,0,0,0);
  read_cache_double.resize(0,0,0,0);

  hfile.reset(new HdfFile(Fname));

  // Read optional metadata
  if (hfile->has_object("Extent_Range")) {
      extent_range.reference(hfile->read_field<double, 2>("Extent_Range"));
      extent_index.reference(hfile->read_field<int, 2>("Extent_Index"));
  }

  // Read data fields
  ArrayWithUnit<double, 1> pgrid_unit =
    hfile->read_field_with_unit<double, 1>("P_level").convert(Unit("Pa"));
  pgrid.reference(pgrid_unit.value);
  ArrayWithUnit<double, 2> tgrid_unit =
    hfile->read_field_with_unit<double, 2>("Temperature").convert(Unit("K"));
  tgrid.reference(tgrid_unit.value);
  // TODO - This might be freq instead of wn. Should perhaps read
  // array with units
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
  std::string bindex = "";
  if(hfile->has_object("H2O_VMR")) {
    bname = "h2o";
    bvmr.reference(hfile->read_field<double, 1>("H2O_VMR"));
  } else if(hfile->has_object("O2_VMR")) {
    bname = "o2";
    bvmr.reference(hfile->read_field<double, 1>("O2_VMR"));
  }
}

// Find the array locations where the wavenumber is contained
// If no extents are defined, this just returns beginning and end of whole absco array
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
        auto extents = wn_extent(wn);
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
  double f = (Wn_in - *(wnptr - 1)) / (*wnptr - *(wnptr - 1));
  if(f > 0.1 && f < 0.9) {
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

// See base class for description
Array<double, 3> AbscoAer::read_double(double Wn_in) const
{
  int wi = wn_index(Wn_in);
  if(wi < cache_double_lbound ||
     wi >= cache_double_ubound)
    swap<double>(wi);
  // The Aer data is in temperature x pressure x broadner order. The
  // base class expects this to be pressure x temperature x broadner,
  // so we transpose this.
  return read_cache<double>()(wi - cache_double_lbound, Range::all(),
      Range::all(), Range::all()).transpose(secondDim, firstDim, thirdDim);
}

// See base class for description
Array<float, 3> AbscoAer::read_float(double Wn_in) const
{
  int wi = wn_index(Wn_in);
  if(wi < cache_float_lbound ||
     wi >= cache_float_ubound)
    swap<float>(wi);
  // The Aer data is in temperature x pressure x broadner order. The
  // base class expects this to be pressure x temperature x broadner,
  // so we transpose this.
  return read_cache<float>()(wi - cache_float_lbound, Range::all(),
      Range::all(), Range::all()).transpose(secondDim, firstDim, thirdDim);
}

//-----------------------------------------------------------------------
/// Make sure the row i is found in the read_cache, possibly reading
/// data if it isn't found.
//-----------------------------------------------------------------------

template<class T> void AbscoAer::swap(int i) const
{
  // First time through, set up space for cache.
 if(read_cache<T>().extent(fourthDim) == 0)
   read_cache<T>().resize(cache_nline, tgrid.cols(), tgrid.rows() - 1,
			  std::max(1, number_broadener_vmr()));
  int nl = read_cache<T>().extent(firstDim);
  // Either read 3d or 4d data. We tell which kind by whether or not
  // we have a number_broadener_vmr() > 0 or not.
  TinyVector<int, 4> start, size;
  start = (i / nl) * nl, 0, 0, 0;
  size = std::min(nl, wngrid.rows() - start(3)), tgrid.cols(), tgrid.rows() - 1,
    std::max(number_broadener_vmr(), 1);
  bound_set<T>((i / nl) * nl, size(3));
  if(number_broadener_vmr() > 0) {
    // 4d case
    read_cache<T>()(Range(0, size(0) - 1), Range::all(),
		    Range::all(), Range::all()) =
      hfile->read_field<T, 4>(field_name, start, size);
  } else {
    // 3d case
    TinyVector<int, 3> start2, size2;
    start2 = 0, 0, (i / nl) * nl;
    size2 = std::min(nl, wngrid.rows() - start(3)), tgrid.cols(),
      tgrid.rows() - 1;
    read_cache<T>()(Range(0, size(3) - 1), Range::all(), Range::all(), 0) =
      hfile->read_field<T, 3>(field_name, start2, size2);
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
