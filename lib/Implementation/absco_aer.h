#ifndef ABSCO_AER_H
#define ABSCO_AER_H
#include "absco.h"
#include "hdf_file.h"
#include "spectral_bound.h"
#include <blitz/array.h>
#include <vector>


namespace FullPhysics {
/****************************************************************//**
  This class is used to read the absco tables. This reads the HDF
  version of the files.

  Note that performance reasons we cache the data as we read it. The
  default cache is about 50 MB, which is a bit large but not too
  large. This can be adjusted if needed, either up for better
  performance or down for less memory.

  For AbscoHdf we constrain the wave numbers we can retrieve on to
  match what is in the Absco tables. This is done because
  interpolation can work very poorly near a peak where the absorption
  is rapidly changing. Also, the Absco tables are a bit of a bottle
  neck and avoiding interpolation in the wavenumber direction is
  desirable. 

  For AbscoAer, the spacing between points is not always
  constant. We'd like to explore the option of interpolation or
  simplest using nearest neighbor. Another option would be to base the
  high resolution grid on what is available in the Absco table (e.g., 
  some variation of a SpectrumSampling class).

  To allow us to explore these options, we have an "interplation_type"
  argument that controls the behavior of the class. Default behavior
  is the original "THROW_ERROR_IF_NOT_ON_WN_GRID", but we can also
  return the nearest neighbor or do a linear interpolation.
*******************************************************************/
class AbscoAer: public Absco {
public:
  enum InterpolationType {THROW_ERROR_IF_NOT_ON_WN_GRID, NEAREST_NEIGHBOR_WN,
			  INTERPOLATE_WN};
  AbscoAer(const std::string& Fname, double Table_scale = 1.0, 
	   int Cache_nline = 5000,
	   InterpolationType Itype = THROW_ERROR_IF_NOT_ON_WN_GRID);
  AbscoAer(const std::string& Fname, 
	   const SpectralBound& Spectral_bound,
	   const std::vector<double>& Table_scale,
	   int Cache_nline = 5000,
	   InterpolationType Itype = THROW_ERROR_IF_NOT_ON_WN_GRID);
  void load_file(const std::string& Fname);
  void load_file(const std::string& Fname, double Table_scale,
		 int Cache_nline = 5000);
  void load_file(const std::string& Fname, 
		 const SpectralBound& Spectral_bound,
		 const std::vector<double>& Table_scale,
		 int Cache_nline = 5000);
  virtual ~AbscoAer() {}
  virtual std::string broadener_name() const { return bname; }
  virtual blitz::Array<double, 1> broadener_vmr_grid() const
  { return bvmr; }

//-----------------------------------------------------------------------
/// Return interpolation type for how we handle wave numbers not on
/// the wave number grid.
//-----------------------------------------------------------------------
  
  InterpolationType interpolation_type() const {return itype_;}

//-----------------------------------------------------------------------
/// Set interpolation type for how we handle wave numbers not on
/// the wave number grid.
//-----------------------------------------------------------------------

  void interpolation_type(InterpolationType Itype) {itype_ = Itype;}
  
  virtual const std::pair<double*, double*> wn_extent(double Wn_in) const;
  virtual void wn_extent(double Wn_in, double& X, double& Y) const;
  virtual double table_scale(double wn) const;
  virtual blitz::Array<double, 1> pressure_grid() const {return pgrid;}
  virtual blitz::Array<double, 2> temperature_grid() const {return tgrid;}
  blitz::Array<double, 1> wavenumber_grid() const { return wngrid; }
  virtual bool have_data(double wn) const;
  virtual bool is_float() const { return is_float_;}
  virtual std::string file_name() const { return hfile->file_name(); } 
  virtual void print(std::ostream& Os) const;
protected:
  virtual blitz::Array<double, 3> read_double(double wn) const;
  virtual blitz::Array<float, 3> read_float(double wn) const;
private:
  bool is_float_;
  int cache_nline;
  InterpolationType itype_;
  blitz::Array<double, 1> bvmr;
  boost::shared_ptr<HdfFile> hfile;
  mutable int cache_double_lbound;
  mutable int cache_double_ubound;
  mutable int cache_float_lbound;
  mutable int cache_float_ubound;
  template<class T> void bound_set(int lbound, int sz) const;
  mutable blitz::Array<double, 5> read_cache_double;
  mutable blitz::Array<float, 5> read_cache_float;
  template<class T> blitz::Array<T, 5>& read_cache() const;
  template<class T> void swap(int i) const;
  int wn_index(double Wn_in) const;
  int wn_index(double Wn_in, double& F) const;
  std::string field_name;
  std::string bname;
  blitz::Array<double, 1> pgrid;
  blitz::Array<double, 2> tgrid;
  blitz::Array<double, 1> wngrid;
  SpectralBound sb;
  std::vector<double> table_scale_;
  blitz::Array<double, 2> extent_range;
  blitz::Array<int, 2> extent_index;
  // This just points to front of wn. We store this here so
  // don't need to keep calculating this.
  double *wnfront;
};

template<> inline blitz::Array<double, 5>& AbscoAer::read_cache<double>() const
{ return read_cache_double; }

template<> inline blitz::Array<float, 5>& AbscoAer::read_cache<float>() const
{ return read_cache_float; }

template<> inline void AbscoAer::bound_set<double>(int lbound, int sz) const
{
  cache_double_lbound = lbound;
  cache_double_ubound = lbound + sz;
}

template<> inline void AbscoAer::bound_set<float>(int lbound, int sz) const
{
  cache_float_lbound = lbound;
  cache_float_ubound = lbound + sz;
}

}
#endif
