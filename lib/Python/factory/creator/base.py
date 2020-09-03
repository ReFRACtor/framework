import logging
import inspect
from collections import OrderedDict

# Use AttrDict as the dictionary class used to store config instantiation
# to make accessing config values easier
from attrdict import AttrDict as ConfigDict

from ..param import ConfigParam, ParamError, AnyValue, Iterable, InstanceOf, Scalar

logger = logging.getLogger('factory.creator.base')

# Other creators subscribed to the output of creators by type
# Maintain this once in the module so the reference to it contained in every single creator
subscribed_observers = OrderedDict()

# Name that an not be used for Creator parameters or else they conflict with class internals
RESERVED_PARAM_NAMES = [
    "parameters",
    "common_store",
    "create",
]

class CreatorError(Exception):
    pass

class ParameterAccessor(object):
    "Proxy class to stand in place of ConfigParam defitions in Creator instances, to allow accessing as if a method of the creator class"

    def __init__(self, param_name, param_def, creator):
        self.param_name = param_name
        self.param_def = param_def
        self.config_def = creator.config_def
        self.common_store = creator.common_store
        self.creator_name = creator.__class__.__name__

    def value(self, **kwargs):
        "Retrieve a parameter from the configuration definition"

        # Get parameter from configuration definition first, then trying common store and
        # if not required return None without error
        if self.param_name in self.config_def:
            param_val = self.config_def[self.param_name]
        elif self.param_name in self.common_store:
            param_val = self.common_store[self.param_name]
        elif self.param_def.required:
            raise KeyError("The parameter name %s requested from config definition by %s was not found and is required" % (self.param_name, self.creator_name))
        else:
            return self.param_def.default

        try:
            return self.param_def.evaluate(param_val, **kwargs)
        except ParamError as exc:
            raise ParamError("The parameter named %s requested by creator %s fails with error: %s" % (self.param_name, self.creator_name, exc))

    def __call__(self, **kwargs):
        return self.value(**kwargs)

class Creator(object):
    "Base creator object that handles ensuring the contract of creators parameters is kept and type checking parameters requested"

    def __init__(self, config_def, common_store=None):

        # Check that no parameters have been defined that will conflict with class internals
        # Do this check before we shadow any of these definitions below
        self._check_parameter_defs()

        # Create a new common store if none are passed
        if common_store is not None:
            if not isinstance(common_store, dict):
                raise TypeError("The common store passed must be a instance of a dict")
            self.common_store = common_store
        else:
            self.common_store = {}

        # Ensure the config definition is of the appropriate type
        if isinstance(config_def, dict):
            self.config_def = self._process_config(config_def)
        else:
            raise TypeError("The config definition (config_def) argument passed to %s should be of type dict" % self.__class__.__name__)

        # Load parameters defined in class definition
        self.parameters = {}
        self._load_parameters()

        logger.debug("Initialized creator %s with parameters: %s" % (self.__class__.__name__, self.parameters.keys()))

    def _process_config(self, in_config_def):
        "Process config definition, create nested Creators"

        out_config_def = ConfigDict()

        for config_name, config_val in in_config_def.items():
            if isinstance(config_val, dict) and "creator" in config_val:
                # Creator in a dictionary block
                logger.debug("Initializing nested creator %s for %s" % (config_name, self.__class__.__name__))
                out_config_def[config_name] = config_val["creator"](config_val, common_store=self.common_store)
            elif inspect.isclass(config_val) and issubclass(config_val, Creator):
                # Creator by itself with no arguments from a dictionary block, all arguments from common store
                logger.debug("Initializing nested creator %s for %s" % (config_name, self.__class__.__name__))
                out_config_def[config_name] = config_val({}, common_store=self.common_store)
            else:
                out_config_def[config_name] = config_val

        return out_config_def

    def register_parameter(self, param_name, param_def):

        if param_name in RESERVED_PARAM_NAMES:
            raise ParamError(f"Can not use '{attr_name}' as a parameter name in the {self.__class__.__name__} creator class as it is a reserved name used internally")

        param_proxy = ParameterAccessor(param_name, param_def, self)
        self.parameters[param_name] = param_proxy
        setattr(self, param_name, param_proxy)

    def _check_parameter_defs(self):
        
         for attr_name in dir(self):
            attr_val = getattr(self, attr_name)
            if isinstance(attr_val, ConfigParam) and attr_name in RESERVED_PARAM_NAMES:
                raise ParamError(f"Can not use '{attr_name}' as a parameter name in the {self.__class__.__name__} creator class as it is a reserved name used internally")
        
    def _load_parameters(self):
        "Gather class parameter types and ensure config_def has necessary items"

        # Look through class attributes to find the ones that have values that
        # inherit from ConfigParam to determine which are the parameters defined.
        # Replace in the object the parameter definition with an ParamAccessor proxy
        for attr_name in dir(self):
            attr_val = getattr(self, attr_name)
            if isinstance(attr_val, ConfigParam):
                self.register_parameter(attr_name, attr_val)

    def param(self, param_name, **kwargs):
        "Retrieve a parameter from the configuration definition"

        try:
            param_proxy = self.parameters[param_name]
        except KeyError:
            # Creator requested a parameter that was not defined by a ConfigParam attribute
            raise KeyError("Unregistered parameter %s requested by %s" % (param_name, self.__class__.__name__))
        else:
            return param_proxy.value(**kwargs)

    def register_to_receive(self, RfType):
        dispatch_list = subscribed_observers.get(RfType, [])
        dispatch_list.append(self)
        subscribed_observers[RfType] = dispatch_list

    def deregister_to_receive(self, RfType=None):

        if RfType is None:
            # Deregister this creator from all observations
            for obj_type, dispatch_list in subscribed_observers.items():
                if self in dispatch_list:
                    dispatch_list.remove(self)
        else:
            dispatch_list = subscribed_observers.get(RfType, [])
            if self in dispatch_list:
                dispatch_list.remove(self)

        # Clean out empty lists of creators 
        types_to_del = []
        for obj_type in subscribed_observers.keys():
            if len(subscribed_observers[obj_type]) == 0:
                types_to_del.append(obj_type)

        for obj_type in types_to_del:
            del subscribed_observers[obj_type]

    def _dispatch(self, rf_obj):
        "Called when the creator's object has been created to be sent to other creators who have registered to recieve objects of that type"

        for RfType, dispatch_list in subscribed_observers.items():
            if isinstance(rf_obj, RfType):
                for observer in dispatch_list:
                    observer.receive(rf_obj)

    def receive(self, rf_obj):
        "Recieves a dispatched object, needs to be implemented in inheriting class"
        pass

    def create(self, **kwargs):
        "Implements the action to create the object referred to by the Creator class. Should not be called directly."
        raise NotImplementedError("Create must be defined in inheriting Creator classes")

    def __call__(self, **kwargs):
        """Turns creators into callables so that they can be evaluated by ConfigParam as any other callable without it
        needing to know any details of this class."""

        result = self.create(**kwargs)
    
        if isinstance(result, list):
            for result_item in result:
                self._dispatch(result_item)
        elif isinstance(result, dict):
            for result_item in result.items():
                self._dispatch(result_item)

        # For iterable types (list, dict), also dispatch the list and dict itself
        # in case these are subclasses that indicate a grouping of data that needs
        # to be captured
        self._dispatch(result)

        # Remove this class from the subscribed observers once its create routine has run
        self.deregister_to_receive()

        return result


class ParamIterateCreator(Creator):
    "Base class for creators that iterate over their parameters"

    order = Iterable(required=False)

    def __init__(self, config_def, param_order=None, **kwargs):
        super().__init__(config_def, **kwargs)

        if param_order is None:
            if "order" in self.config_def:
                param_order = self.param("order")
            else:
                param_order = self.config_def.keys()

        # Filter out parameter names that should not be processed,
        # such as the creator param
        self.param_names = []
        for param_name in param_order:
            if param_name != "creator":
                self.param_names.append(param_name)

                # Param type not defined so allow any value
                if param_name not in self.parameters:
                    self.register_parameter(param_name, AnyValue())

class ParamPassThru(ParamIterateCreator):
    "Evaluates and passes configurations parameter through as the creator result"

    def create(self, **kwargs):

        result = ConfigDict()
        for param_name in self.param_names:
            logger.debug("Passing through parameter %s" % param_name)
            ret_obj = self.param(param_name, **kwargs)
            self._dispatch(ret_obj)
            result[param_name] = ret_obj

        return result

class SaveToCommon(ParamPassThru):
    "Evaluates parameters and saves them into the common store, creator has no return value"

    def create(self, **kwargs):

        result = ConfigDict()
        for param_name in self.param_names:
            logger.debug("Saving to the common store parameter %s" % param_name)
            ret_obj = self.param(param_name, **kwargs)
            self._dispatch(ret_obj)
            result[param_name] = ret_obj
            self.common_store[param_name] = ret_obj

        return result

class PickChild(Creator):
    "Simply picks one of many nested child elements and returns that, can be used as a place holder before implementing some sort of choice logic Creator"

    child = Scalar(str)

    def create(self, **kwargs):
        self.register_parameter(self.child(), AnyValue())
        ret_obj = self.param(self.child(), **kwargs)
        self._dispatch(ret_obj)
        return ret_obj
