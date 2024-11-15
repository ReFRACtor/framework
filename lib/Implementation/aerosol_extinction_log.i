// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "aerosol_extinction_log.h"
  %}
%base_import(aerosol_extinction_level)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AerosolExtinctionLog)

// Force to be not abstract, SWIG had troubles seeing that the clone methods ARE implemented below
%feature("notabstract") AerosolExtinctionLog;

namespace FullPhysics {
class AerosolExtinctionLog : public AerosolExtinctionLevel {
public:
  AerosolExtinctionLog(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name);

  virtual ~AerosolExtinctionLog() = default;
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
  %pickle_serialization();
};
}

