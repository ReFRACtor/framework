import numpy as np

class ParamError(BaseException):
    pass

class ConfigParam(object):
    "Base class for configuration parameters"

    def __init__(self, required=True):
        self.required = required

    def evaluate(self, in_value, creator):

        # If value to be evaluated is callable, call and pass
        # along the creator requesting the parameter
        if callable(in_value):
            out_value = in_value(creator)
        elif isinstance(in_value, dict) and "creator" in in_value:
            out_value = in_value["creator"](in_value, common_store=creator.common_store).create()
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

class Scalar(ConfigParam):
    "Configuration parameter that resolves to a single scalar value"

    def __init__(self, dtype=None, **kwargs):
        super().__init__(**kwargs)

        self.dtype = dtype

    def check_type(self, value):
        if not np.isscalar(value):
            raise ParamError("Value is not a scalar: %s" % value)

        if self.dtype is not None and not isinstance(type(value), self.dtype):
            raise ParamError("Type of parameter %s does not match expected %s for value: %s" % (type(value), self.dtype, value))

class Array(ConfigParam):
    "Configuration parameter that resolves to a numpy array with optional checking on dimensionality and units"
    
    def __init__(self, dims=None, dtype=None, has_units=None, **kwargs):
        super().__init__(**kwargs)

        self.dims = dims
        self.dtype = dtype
        self.has_units = has_units

    def check_type(self, value):

        if not isinstance(value, np.ndarray):
            raise ParamError("Value us not a numpy array: %s" % value)

        if self.dims is not None and len(value.shape) != self.dims:
            raise ParamError("Number of dimensions %d do not match the expected number %s for value: %s" % (len(value.shape), self.dims, value))

        if self.dtype is not None and not isinstance(value.dtype, self.dtype):
            raise ParamError("Type of parameter %s does not match expected %s for value: %s" % (value.dtype, self.dtype, value))

class Iterable(ConfigParam):
    "Configuration parameter that resolve to an iterable object like a tuple or list, but that is not required to be an array"

    def check_type(self, value):

        if not hasattr(value, "__iter__"):
            raise ParamError("Expected an iterable for value: %s" % value)

class InstanceOf(ConfigParam):
    "Configuration parameter that must be an instance of a specific type of class"

    def __init__(self, cls_type, **kwargs):
        super().__init__(**kwargs)
        self.cls_type = cls_type

    def check_type(self, value):

        if not isinstance(value, self.cls_type):
            raise ParamError("Expected an instance of %s for value: %s" % (self.cls_type, value))
