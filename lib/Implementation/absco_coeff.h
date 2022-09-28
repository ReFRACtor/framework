#ifndef ABSCO_COEFF_H
#define ABSCO_COEFF_H
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

  This is similar to AbscoHdf, but instead of directly having the
  cross section tables in the file we have coefficients and
  representative vector (so analogous to doing a principle components
  analysis). We do this to compress the data files.
*******************************************************************/
class AbscoCoeff: public Absco {
public:
  AbscoCoeff(const std::string& Fname, double Table_scale = 1.0, 
	     int Cache_nline = 50000);
  AbscoCoeff(const std::string& Fname, 
	     const SpectralBound& Spectral_bound,
	     const std::vector<double>& Table_scale,
	     int Cache_nline = 50000);
  void load_file(const std::string& Fname);
  void load_file(const std::string& Fname, double Table_scale,
		 int Cache_nline = 5000);
  void load_file(const std::string& Fname, 
		 const SpectralBound& Spectral_bound,
		 const std::vector<double>& Table_scale,
		 int Cache_nline = 5000);
  virtual ~AbscoCoeff() {}
  virtual int number_broadener() const
  {
    return (int) bname.size();
  }
  virtual std::string broadener_name(int Broadener_index) const
  {
    range_check(Broadener_index, 0, number_broadener());
    return bname[Broadener_index];
  }
  virtual blitz::Array<double, 1> broadener_vmr_grid(int Broadener_index) const
  {
    range_check(Broadener_index, 0, number_broadener());
    return bvmr[Broadener_index];
  }
  virtual const std::pair<double*, double*> wn_extent(double Wn_in) const;
  virtual void wn_extent(double Wn_in, double& X, double& Y) const;
  virtual double table_scale(double wn) const;
  virtual blitz::Array<double, 1> pressure_grid() const {return pgrid;}
  virtual blitz::Array<double, 2> temperature_grid() const {return tgrid;}
  blitz::Array<double, 1> wavenumber_grid() const { return wngrid; }
  virtual bool have_data(double wn) const;
  void wn_index(double Wn_in, int& Wn_index, double& F) const;
  virtual bool is_float() const { return is_float_;}
  virtual std::string file_name() const { return hfile->file_name(); } 
  virtual void print(std::ostream& Os) const;

  virtual blitz::Array<double, 3> read_double(double wn) const;
  virtual blitz::Array<float, 3> read_float(double wn) const;
private:
  void load_file();
  bool is_float_;
  int cache_nline;
  std::vector<blitz::Array<double, 1> > bvmr;
  boost::shared_ptr<HdfFile> hfile;
  mutable int cache_double_lbound;
  mutable int cache_double_ubound;
  mutable int cache_float_lbound;
  mutable int cache_float_ubound;
  template<class T> void bound_set(int lbound, int sz) const;
  mutable blitz::Array<double, 2> read_cache_double;
  mutable blitz::Array<float, 2> read_cache_float;
  template<class T> blitz::Array<T, 2>& read_cache() const;
  template<class T> void swap(int i) const;
  std::string field_name;
  std::vector<std::string> bname;
  blitz::Array<double, 1> pgrid;
  blitz::Array<double, 2> tgrid;
  blitz::Array<double, 3> cross_sec_coeff;
  blitz::Array<double, 1> wngrid;
  SpectralBound sb;
  std::vector<double> table_scale_;
  blitz::Array<double, 2> extent_range;
  blitz::Array<int, 2> extent_index;
  // This just points to front of wn. We store this here so
  // don't need to keep calculating this.
  double *wnfront;
  AbscoCoeff() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};

template<> inline blitz::Array<double, 2>& AbscoCoeff::read_cache<double>() const
{ return read_cache_double; }

template<> inline blitz::Array<float, 2>& AbscoCoeff::read_cache<float>() const
{ return read_cache_float; }

template<> inline void AbscoCoeff::bound_set<double>(int lbound, int sz) const
{
  cache_double_lbound = lbound;
  cache_double_ubound = lbound + sz;
}

template<> inline void AbscoCoeff::bound_set<float>(int lbound, int sz) const
{
  cache_float_lbound = lbound;
  cache_float_ubound = lbound + sz;
}

}

FP_EXPORT_KEY(AbscoCoeff);
#endif
