import re

import numpy as np

from refractor import framework as rf

class ParamError(BaseException):
    pass

class ConfigParam(object):
    "Base class for configuration parameters"

    def __init__(self, required=True, default=None):
        self.accessor_func = None
        self.default = default

        if default is not None:
            self.required = False
        else:
            self.required = required

    def evaluate(self, in_value, **kwargs):

        # If value to be evaluated is callable, call and pass
        # along the creator requesting the parameter
        if callable(in_value):
            out_value = in_value(**kwargs)
        else:
            out_value = in_value

        self.check_type(out_value)

        return out_value

    def check_type(self, value):
        raise NotImplementedError("check_type must be implemented in inheriting class")

class AnyValue(ConfigParam):
    "Bypasses type checking for the parameter"

    def check_type(self, value):
        pass

class NoneValue(ConfigParam):
    "Configuration parameter that is equal to None, only useful when combined with the Choice parameter type"

    def check_type(self, value):
        if value is not None:
            raise ParamError("Parameter was expected to be None")

class Choice(ConfigParam):
    "Allows the choice of multiple parameter types"

    def __init__(self, *vargs, **kwargs):
        super().__init__(**kwargs)

        self.choices = []
        for choice in vargs:
            if isinstance(choice, ConfigParam):
                self.choices.append(choice)
            else:
                raise ParamError("Arguments to Choice must be instance of ConfigParam types.")

    def check_type(self, value):

        valid = False
        for choice in self.choices:
            try:
                choice.check_type(value)
            except ParamError:
                pass
            else:
                valid = True
                break

        if not valid:
            raise ParamError("Parameter value does not match any of the parameter choices")

class Scalar(ConfigParam):
    "Configuration parameter that resolves to a single scalar value"

    def __init__(self, dtype=None, **kwargs):
        super().__init__(**kwargs)

        self.dtype = dtype

    def check_type(self, value):
        if not np.isscalar(value):
            raise ParamError("Value is not a scalar: %s" % value)

        if self.dtype is not None and not isinstance(value, self.dtype):
            raise ParamError("Type of parameter %s does not match expected %s for value: %s" % (type(value), self.dtype, value))

class Array(ConfigParam):
    "Configuration parameter that resolves to a numpy array with optional checking on dimensionality and units"
    
    def __init__(self, dims=None, dtype=None, **kwargs):
        super().__init__(**kwargs)

        self.dims = dims
        self.dtype = dtype

    def check_type(self, value):

        if not isinstance(value, np.ndarray):
            raise ParamError("Value is not a numpy array: %s" % value)

        if self.dims is not None and len(value.shape) != self.dims:
            raise ParamError("Number of dimensions %d do not match the expected number %s for value: %s" % (len(value.shape), self.dims, value))

        if self.dtype is not None and not isinstance(value.dtype, self.dtype):
            raise ParamError("Type of parameter %s does not match expected %s for value: %s" % (value.dtype, self.dtype, value))

class Iterable(ConfigParam):
    "Configuration parameter that resolve to an iterable object like a tuple or list, but that is not required to be an array"

    def __init__(self, val_type=None, **kwargs):
        super().__init__(**kwargs)
        self.val_type = val_type

        if val_type is not None and not isinstance(val_type, type):
            raise ParamError(f"Type supplied must be a type object or None, not {val_type}")

    def check_type(self, value):

        if not hasattr(value, "__iter__"):
            raise ParamError(f"Expected an iterable for value: {value}")

        if self.val_type is not None:
            for idx, iter_val in enumerate(value):
                if not isinstance(iter_val, self.val_type):
                    raise ParamError(f"Expected an instance of {self.val_type} for value {iter_val} at index {idx} of iterable")

class InstanceOf(ConfigParam):
    "Configuration parameter that must be an instance of a specific type of class"

    def __init__(self, val_type, **kwargs):
        super().__init__(**kwargs)
        self.val_type = val_type

    def check_type(self, value):

        if not isinstance(value, self.val_type):
            raise ParamError("Expected an instance of %s for value: %s" % (self.val_type, value))

class Dict(InstanceOf):
    "Configuration parameter that resolves to a dict"

    def __init__(self, **kwargs):
        super().__init__(dict, **kwargs)

class ArrayWithUnit(ConfigParam):

    awu_type_pattern = r"ArrayWithUnit_\w+_\d+"

    def __init__(self, dims=None, dtype=None, **kwargs):
        super().__init__(**kwargs)

        self.dims = dims
        self.dtype = dtype

    def check_type(self, value):

        # Check that type string matches the pattern, can't just use ininstance
        # on the SWIG type
        type_str = str(type(value))

        if not re.search(self.awu_type_pattern, type_str):
            raise ParamError("Value is not an ArrayWithUnit type: %s" % type(value))

        # Check dims and dtype against the value within the awu
        if self.dims is not None and len(value.value.shape) != self.dims:
            raise ParamError("Number of dimensions %d do not match the expected number %s for ArrayWithUnit value: %s" % (len(value.value.shape), self.dims, value))

        if self.dtype is not None and not isinstance(value.value.dtype, self.dtype):
            raise ParamError("Type of parameter %s does not match expected %s for ArrayWithUnit value: %s" % (value.value.dtype, self.dtype, value))

class DoubleWithUnit(InstanceOf):
    "Short cut for the DoubleWithUnit type to parallel the ArrayWithUnit type parameter"

    def __init__(self, **kwargs):
        super().__init__(rf.DoubleWithUnit, **kwargs)

class ObjectVector(ConfigParam):
    "Checks that value is a C++ vector of a certain type"

    def __init__(self, vec_type=None, **kwargs):
        super().__init__(**kwargs)

        self.vec_type = vec_type

    def check_type(self, value):
        
        # Check that the type string has vector_ in it and vec_type if that is supplied
        if self.vec_type is not None:
            check_str = r"\.vector_%s" % self.vec_type
        else:
            check_str = r"\.vector_.*"

        type_str = str(type(value))

        if not re.search(check_str, type_str):
            raise ParamError("Value with type string %s is not a C++ vector with type vector_%s" % (type_str, self.vec_type and self.vec_type or ""))
