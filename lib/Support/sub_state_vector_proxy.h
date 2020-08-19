#ifndef SUB_STATE_VECTOR_PROXY_H
#define SUB_STATE_VECTOR_PROXY_H

#include "sub_state_vector_observer.h"

namespace FullPhysics {

/****************************************************************//**
  Proxies for multiple SubStateVectorObserver classes that are
  combined for a result, but handle their state vector data
  seperately. Can only be used as a base for another class.

  Proxied classes should be pushed into proxied_observers in
  consturctor.
*******************************************************************/
class SubStateVectorProxy : virtual public SubStateVectorObserver {
public:
  virtual ~SubStateVectorProxy() {}
  virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub,
				const blitz::Array<double, 2>& Cov_sub);
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name)
    const;
  virtual void print(std::ostream& Os) const;
protected:
  // Ensure that class can only be inherited
  SubStateVectorProxy() {}
  void initialize
  (const std::vector<boost::shared_ptr<SubStateVectorObserver> >& Proxied);
private:
  std::vector<boost::shared_ptr<SubStateVectorObserver> > proxied_observers;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(SubStateVectorProxy);
#endif
