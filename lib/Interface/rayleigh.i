% include "fp_common.i"

% {
#include "rayleigh.h"
    %
}

% fp_shared_ptr(FullPhysics::Rayleigh);

namespace FullPhysics {

    class Rayleigh: public Printable<Rayleigh> {
    public:

        virtual ArrayAd<double, 1> optical_depth_each_layer(double wn, int spec_index) const;

        virtual void print(std::ostream& Os) const;

    };
}
