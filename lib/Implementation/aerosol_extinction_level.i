// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "aerosol_extinction_level.h"
%}

%base_import(aerosol_extinction_imp_base)
%base_import(mapping)
%base_import(mapping_linear)

%import "pressure.i"
%fp_shared_ptr(FullPhysics::AerosolExtinctionLevel)

// Force to be not abstract, SWIG had troubles seeing that the clone methods ARE implemented below
%feature("notabstract") AerosolExtinctionLevel;

namespace FullPhysics {
class AerosolExtinctionLevel : public AerosolExtinctionImpBase {
public:
  AerosolExtinctionLevel(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<bool, 1>& Flag, 
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name,
			  boost::shared_ptr<Mapping> in_map = boost::make_shared<MappingLinear>());
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
  virtual std::string state_vector_name_i(int i) const;
protected:
  virtual void calc_aerosol_extinction() const;
};
}

