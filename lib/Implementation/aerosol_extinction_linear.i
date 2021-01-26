// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "aerosol_extinction_linear.h"
%}

%base_import(aerosol_extinction_level)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AerosolExtinctionLinear)

// Force to be not abstract, SWIG had troubles seeing that the clone methods ARE implemented below
%feature("notabstract") AerosolExtinctionLinear;

namespace FullPhysics {
class AerosolExtinctionLinear : public AerosolExtinctionLevel {
public:
  AerosolExtinctionLinear(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name);

  virtual ~AerosolExtinctionLinear() = default;
    
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
  %pickle_serialization();
};
}

