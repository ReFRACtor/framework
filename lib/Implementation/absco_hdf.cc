#include "absco_hdf.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include <iomanip>
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbscoHdf::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Absco)
    & FP_NVP(cache_nline) & FP_NVP_(itype) & FP_NVP(hfile)
    & FP_NVP(sb) & FP_NVP_(table_scale);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void AbscoHdf::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void AbscoHdf::load(Archive & UNUSED(ar),
		    const unsigned int UNUSED(version))
{
  load_file();
}

FP_IMPLEMENT(AbscoHdf);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
namespace FullPhysics {
void register_lua_AbscoHdf(lua_State *ls) { \
luabind::module(ls) [ 
luabind::class_<AbscoHdf,Absco,boost::shared_ptr<GasAbsorption> >
("AbscoHdf")
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

AbscoHdf::AbscoHdf(const std::string& Fname, double Table_scale,
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

AbscoHdf::AbscoHdf(const std::string& Fname, 
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

void AbscoHdf::load_file(const std::string& Fname)
{
  load_file(Fname, sb, table_scale_, cache_nline);
}

//-----------------------------------------------------------------------
/// Read the given Absco file. Note that you don't normally call this
/// directly, but rather use the constructor. However there are
/// certain uses where changing the absco table can be useful, such
/// as the python code l2_fp_fm_errors.py.
//-----------------------------------------------------------------------

void AbscoHdf::load_file(const std::string& Fname, double Table_scale, 
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

void AbscoHdf::load_file(const std::string& Fname, 
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
 
void AbscoHdf::load_file()
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
  read_cache_float.resize(0,0,0,0);
  read_cache_double.resize(0,0,0,0);

  // Read optional metadata
  if (hfile->has_object("Extent_Range")) {
      extent_range.reference(hfile->read_field<double, 2>("Extent_Range"));
      extent_index.reference(hfile->read_field<int, 2>("Extent_Index"));
  }

  // Read data fields
  pgrid.reference(hfile->read_field<double, 1>("Pressure"));
  tgrid.reference(hfile->read_field<double, 2>("Temperature"));
  wngrid.reference(hfile->read_field<double, 1>("Wavenumber"));
  wnfront = &wngrid(0);

  // This may have a "\0" in it, so we create a string twice, the
  // second one ends at the first "\0".
  std::string t = hfile->read_field<std::string>("Gas_Index");
  field_name = std::string("Gas_") + t.c_str() + "_Absorption";
  // Determine if data is float or double. We use this to optimize the
  // read. 
  H5::DataSet d = hfile->h5_file().openDataSet(field_name);
  if(d.getDataType().getSize() == 4)
    is_float_ = true;
  else
    is_float_ = false;
  std::string bindex = "";
  try {
    bindex = hfile->read_field<std::string>("Broadener_Index");
    bindex = std::string(bindex.c_str()); // In case this ends in "\0"
  } catch(const Exception& E) {
    // It is ok if we don't have a broadener, this just means that we
    // are reading a 3D file.
  }
  if(bindex == "")
    ;                                // Nothing to do
  else if(bindex == "01")        
    // This is HITRAN index 01, which is H2O. This is documented in
    // HITRAN papers, such as "The HITRAN 2008 molecular spectroscopic
    // database", Journal of Quantitative Spectroscopy and Radiative
    // Transfer, vol. 110, pp. 533-572 (2009)
    bname.push_back("h2o");
  else {
    Exception e;
    e << "Right now we only support H2O as a broadener (HITRAN index \"01\"). "
      << "File " << hfile->file_name() << " has index of \""
      << bindex << "\"";
    throw e;
  }

  if(bname.size() > 0)
    bvmr.push_back(hfile->read_field<double, 1>("Broadener_" + bindex + "_VMR").copy());
}


//-----------------------------------------------------------------------
/// Find the array locations where the wavenumber is contained If no
/// extents are defined, this just returns beginning and end of whole
/// absco array. (version friendlier to swig/python)
//-----------------------------------------------------------------------

void AbscoHdf::wn_extent(double Wn_in, double& X, double& Y) const
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

const std::pair<double*, double*> AbscoHdf::wn_extent(double Wn_in) const
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

double AbscoHdf::table_scale(double wn) const 
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
bool AbscoHdf::have_data(double wn) const
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
int AbscoHdf::wn_index(double Wn_in) const
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
      << "AbscoHDF does not interpolate in wavenumber direction.\n"
      << "The interpolation doesn't work well near peaks, so we\n"
      << "just don't do it. The requested Wn needs to be within 10%\n"
      << "of the absco wn grid, return the value for the closest number.\n"
      << "Wn:        " << Wn_in << "\n"
      << "Wnptr - 1: " << *(wnptr - 1) << "\n"
      << "Wnptr:     " << *wnptr << "\n"
      << "Frac:      " << f << " (should be < 0.1 or > 0.9)\n"
      << "AbscoHdf:\n"
      << *this << "\n";
    throw e;
  }
  if(f < 0.5)                        // First point is closest.
    --wnptr;
  return (int) (wnptr - wnfront);
}

// See base class for description
Array<double, 3> AbscoHdf::read_double(double Wn_in) const
{
  int wi = wn_index(Wn_in);
  if(wi < cache_double_lbound ||
     wi >= cache_double_ubound)
    swap<double>(wi);
  return read_cache<double>()(Range::all(), Range::all(), Range::all(), 
                              wi - cache_double_lbound);
}

// See base class for description
Array<float, 3> AbscoHdf::read_float(double Wn_in) const
{
  int wi = wn_index(Wn_in);
  if(wi < cache_float_lbound ||
     wi >= cache_float_ubound)
    swap<float>(wi);
  return read_cache<float>()(Range::all(), Range::all(), Range::all(), 
                             wi - cache_float_lbound);
}

//-----------------------------------------------------------------------
/// Make sure the row i is found in the read_cache, possibly reading
/// data if it isn't found.
//-----------------------------------------------------------------------

template<class T> void AbscoHdf::swap(int i) const
{
  // First time through, set up space for cache.
 if(read_cache<T>().extent(fourthDim) == 0)
      read_cache<T>().resize(tgrid.rows(), tgrid.cols(),
			     (number_broadener() > 0 ? number_broadener_vmr(0) :
			      1), cache_nline);
  int nl = read_cache<T>().extent(fourthDim);
  // Either read 3d or 4d data. We tell which kind by whether or not
  // we have a number_broadener_vmr() > 0 or not.
  TinyVector<int, 4> start, size;
  start = 0, 0, 0, (i / nl) * nl;
  size = tgrid.rows(), tgrid.cols(), (number_broadener() > 0 ? number_broadener_vmr(0) : 1),
    std::min(nl, wngrid.rows() - start(3));
  bound_set<T>((i / nl) * nl, size(3));
  if(number_broadener() > 0) {
    // 4d case
    read_cache<T>()(Range::all(), Range::all(), Range::all(), Range(0, size(3) - 1)) =
      hfile->read_field<T, 4>(field_name, start, size);
  } else {
    // 3d case
    TinyVector<int, 3> start2, size2;
    start2 = 0, 0, (i / nl) * nl;
    size2 = tgrid.rows(), tgrid.cols(),
      std::min(nl, wngrid.rows() - start(3));
    read_cache<T>()(Range::all(), Range::all(), 0, Range(0, size(3) - 1)) =
      hfile->read_field<T, 3>(field_name, start2, size2);
  }
}

void AbscoHdf::print(std::ostream& Os) const
{
  Os << "AbscoHdf" << "\n"
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
