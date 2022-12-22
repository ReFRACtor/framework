%include "fp_common.i"

%{
#include "state_mapping_at_indexes.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingAtIndexes);

namespace FullPhysics {

%pythonbegin %{
import numpy as np
%}

%feature("pythonprepend") StateMappingAtIndexes::StateMappingAtIndexes %{
  # Make sure argument is a numpy array to avoid unintented consequences
  if not isinstance(args[0], np.ndarray):
    raise ValueError(f"StateMappingAtIndexes argument 0 must be a numpy array iwth dtype bool or integer, not value with type {type(args[0])}")
  elif hasattr(args[0], "dtype") and not np.issubdtype(args[0].dtype, np.integer) and not np.issubdtype(args[0].dtype, bool):
    raise ValueError(f"StateMappingAtIndexes argument 0 must be a numpy array with dtype bool or integer, not a numpy array with dtype {args[0].dtype}")

  # Convert any integer arrays to bool arrays due to operator overloading issues
  # on the C++ side that will cast integer arrays into bool arrays due to 
  # type precedence, the int array method will never be called
  if hasattr(args[0], "dtype") and np.issubdtype(args[0].dtype, np.integer):
    # Create a flag array up to the size of maximum index
    max_index = np.max(args[0])
    flags = np.zeros(max_index+1, dtype=bool)

    # Set flag true at indexes
    flags[args[0]] = True

    # Assign flag as new argument
    args = (flags,)
%}

class StateMappingAtIndexes : public StateMapping  {
public:
  StateMappingAtIndexes(const blitz::Array<int, 1>& indexes);
  StateMappingAtIndexes(const blitz::Array<bool, 1>& flags);

  %python_attribute(retrieval_indexes, const blitz::Array<int, 1>);
  virtual boost::shared_ptr<StateMapping> clone();
  %pickle_serialization();
};
}
