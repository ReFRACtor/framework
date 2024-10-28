// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "ground.h"
#include "pressure.h"
%}

%base_import(observer)
%base_import(state_vector)
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::Ground)
namespace FullPhysics {
  class Ground;
}
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Ground>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Ground>)

namespace FullPhysics {
%template(ObservableGround) FullPhysics::Observable<FullPhysics::Ground>;
%template(ObserverGround) FullPhysics::Observer<FullPhysics::Ground>;

enum SpurrBrdfType {
  LAMBERTIAN  = 1,
  ROSSTHIN    = 2,
  ROSSTHICK   = 3,
  LISPARSE    = 4,
  LIDENSE     = 5,
  HAPKE       = 6,
  ROUJEAN     = 7,
  RAHMAN      = 8,
  COXMUNK     = 9,
  BREONSOIL   = 10,
  BREONVEG    = 11,
  BPDFNDVI    = 12,
  NEWCMGLINT  = 13,
  RTKHOTSPOT  = 14,
  MODFRESNEL  = 15,
  SNOWBRDF    = 16,
  EMISSIVITY  = 100, // Not a real index from LIDORT family, but
                     // included to meet interface
  UNKNOWN = -999     // Type that will trigger an error, if we
                     // have a ground type that doesn't map to
                     // the spurr BRDF
};
  
class Ground : public Observable<Ground> {
public:
  virtual ~Ground();
  virtual void add_observer(Observer<Ground>& Obs); 
  virtual void remove_observer(Observer<Ground>& Obs);
  std::string print_to_string() const;
  std::string print_parent() const;
  %python_attribute(spurr_brdf_type, SpurrBrdfType);
  virtual ArrayAd<double, 1> surface_parameter
    (const double wn, const int spec_index) const = 0;
  virtual boost::shared_ptr<Ground> clone() const;
  virtual void print(std::ostream& Os) const;
  %pickle_serialization();
};
}

