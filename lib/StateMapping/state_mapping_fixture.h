#ifndef STATE_MAPPING_FIXTURE_H
#define STATE_MAPPING_FIXTURE_H

#include "pressure_sigma.h"

namespace FullPhysics {

class StateMappingInterpFixture: public GlobalFixture {
public:
    StateMappingInterpFixture()
    {
        blitz::Array<double, 1> a_from(5), b_from(5);
        a_from = 0;
        b_from = 0.1  , 0.325, 0.55 , 0.775, 1.0;
        press_from.reset(new PressureSigma(a_from, b_from, 1.0));

        blitz::Array<double, 1> a_to(10), b_to(10);
        a_to = 0;
        b_to = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0;
        press_to.reset(new PressureSigma(a_to, b_to, 1.0));

        blitz::Array<double, 1> vals_from(5);
        vals_from = 1.  ,  3.25,  5.5 ,  7.75, 10.0;
        vals_from_ad = ArrayAd<double, 1>(vals_from);
    }

    boost::shared_ptr<PressureSigma> press_from;
    boost::shared_ptr<PressureSigma> press_to;
    ArrayAd<double, 1> vals_from_ad;
};

}

#endif
