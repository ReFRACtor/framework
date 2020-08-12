// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"

%{
#include "ils_imp_base.h"
#include "sub_state_vector_array.h"
%}

%import "double_with_unit.i"
%import "array_ad.i"
%import "spectral_domain.i"
%import "double_with_unit.i"

%base_import(ils)
%base_import(observer)
%base_import(sample_grid)

%fp_shared_ptr(FullPhysics::IlsImpBase);

%feature("director") FullPhysics::IlsImpBase;

namespace FullPhysics {
class IlsImpBase : public Ils, public Observer<SampleGrid> {
public:
  IlsImpBase(const boost::shared_ptr<SampleGrid>& Sample_grid, const DoubleWithUnit& Edge_extension);
  virtual void notify_update(const SampleGrid& D);
  virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;
  virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;
  virtual std::string band_name() const;
  virtual std::string hdf_band_name() const;
  virtual SpectralDomain pixel_grid() const;
  virtual DoubleWithUnit high_res_extension() const;
  virtual void high_res_extension(const DoubleWithUnit& extension);
  boost::shared_ptr<SampleGrid> sample_grid() const {return sample_grid_; }
    virtual boost::shared_ptr<Ils> clone() const = 0;
  %pickle_serialization();
protected:
  IlsImpBase();
};
}

// We may well move this into a central macro. But haven't worked out
// yet exactly what we want here, so we'll do all of this inline for now.
%{
#include "ils_imp_basePYTHON_wrap.h"
#include "fp_serialize_support.h"
  
namespace boost {
   namespace serialization {
     template<class Archive>
     void serialize(Archive & ar, SwigDirector_IlsImpBase & t, 
		    const unsigned int UNUSED(version))
     {
       ar &  boost::serialization::make_nvp(BOOST_PP_STRINGIZE(IlsImpBase), 
	    boost::serialization::base_object<FullPhysics::IlsImpBase>(t));
     }
     template<class Archive>
     void save_construct_data(Archive & ar, const SwigDirector_IlsImpBase* t, 
			      const unsigned int UNUSED(version))
     {
       // PyObject* obj = PyObject_CallMethodObjArgs(d->swig_get_self(),
       // 			PyString_FromString("__boost_serialize_save__"),
       // 			NULL);
       // if(PyErr_Occurred()) {
       // 	GeoCal::Exception e;
       // 	e << "Python error occurred:\n"
       // 	  << parse_python_exception();
       // 	throw e;
       // }
       //std::string python_object = PyString_AsString(obj);
       std::string python_object = "hi there";
       ar & FP_NVP(python_object);
     }
     template<class Archive>
     void load_construct_data(Archive & ar, SwigDirector_IlsImpBase* t,
			      const unsigned int UNUSED(version))
     {
       throw FullPhysics::Exception("Don't have this working yet");
     }
   }
}

BOOST_CLASS_EXPORT_KEY(SwigDirector_IlsImpBase);
BOOST_CLASS_EXPORT_IMPLEMENT(SwigDirector_IlsImpBase);

template void boost::serialization::serialize
(boost::archive::polymorphic_oarchive& ar, SwigDirector_IlsImpBase& t,
 const unsigned int version);
template void boost::serialization::serialize
 (boost::archive::polymorphic_iarchive& ar, SwigDirector_IlsImpBase& t,
  const unsigned int version);

template void boost::serialization::save_construct_data
(boost::archive::polymorphic_oarchive & ar, const SwigDirector_IlsImpBase* t,
 const unsigned int version);

template void boost::serialization::load_construct_data
(boost::archive::polymorphic_iarchive & ar, SwigDirector_IlsImpBase* t,
 const unsigned int version);
 
%}
  
