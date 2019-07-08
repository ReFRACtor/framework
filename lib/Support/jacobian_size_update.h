#ifndef JACOBIAN_SIZE_UPDATE
#define JACOBIAN_SIZE_UPDATE

#include "state_vector_observer.h"
#include "jacobian_size_mixin.h"

namespace FullPhysics {

/****************************************************************//**
 This class updates the jacobian size for users of the 
 JacobianSizeMixin interface whenever the state vector changes.
*******************************************************************/

class JacobianSizeUpdate : virtual public StateVectorObserver {
public:
    /// Add a JacobianSizeMixin class so its jacobian size gets updated when the StateVector changes
    void add_update_target(const boost::shared_ptr<JacobianSizeMixin>& target) { update_targets.push_back(target); }
    virtual void notify_update(const StateVector& sv) { 
        int nvar = sv.state().rows();
        for (int i = 0; i < update_targets.size(); i ++) {
            update_targets[i]->set_jacobian_size(nvar);
        }
    };

private:
    std::vector<boost::shared_ptr<JacobianSizeMixin> > update_targets;
};

}

#endif
