#ifndef OBSERVATION_ID_H
#define OBSERVATION_ID_H

namespace FullPhysics {
/****************************************************************//**
  For use with L1B readers, it is useful to have a base class
  for Observation IDs so the abstract base L1B reader class can
  specify it as the generic return type for some of its functions.
*******************************************************************/
class ObservationId {
public:
  // Have a virtual member function, which forces RTTI information to
  // be available.
  virtual ~ObservationId() {}
};

}
#endif
