#ifndef RAYLEIGH_H
#define RAYLEIGH_H

#include "printable.h"
#include "array_ad.h"

namespace FullPhysics {

/****************************************************************//**
  This class represents the interface for the calculation of
  the Rayleigh portion of the optical depth
*******************************************************************/

class Rayleigh: public Printable<Rayleigh> {
public:

    virtual ArrayAd<double, 1> optical_depth_each_layer(double wn, int spec_index) const = 0;

    virtual void print(std::ostream& Os) const
    {
        Os << "Rayleigh";
    }

};
}
#endif
