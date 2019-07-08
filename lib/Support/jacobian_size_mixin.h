#ifndef JACOBIAN_SIZE_MIXIN
#define JACOBIAN_SIZE_MIXIN

#include <boost/optional.hpp>

#include "fp_exception.h"

namespace FullPhysics {
/****************************************************************//**
 This mixin marks a class as needing to be supplied with the 
 number of jacobian variables used in the construction of ArrayAd
 or AutoDerivative output values. This class contains the machinery
 for handling the tracking of that number. It will throw an error
 if the value has not yet been assigned and there is an attempt
 to access it.
*******************************************************************/

class JacobianSizeMixin {
public:
    virtual ~JacobianSizeMixin() = default;

    /// Set the number of jacobian variables the inheriting class should use
    virtual void set_jacobian_size(int number_variable) { number_variable_ = number_variable; }

    virtual bool is_jacobian_size_set() const { return bool(number_variable_); }

    /// Return the number of jacobian variables and thrown an exception when
    /// it has not yet been defined
    virtual const int jacobian_size() const {
        if (!number_variable_) {
            throw Exception("Number of jacobian variables has not yet been set");
        }

        return *number_variable_;
    }
private:
    boost::optional<int> number_variable_;
};

}

#endif
