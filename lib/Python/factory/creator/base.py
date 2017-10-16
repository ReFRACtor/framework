import logging

from ..param import ConfigParam, ParamError, AnyValue, Iterable, InstanceOf

logger = logging.getLogger('factory.creator.base')

class ParameterAccessor(object):
    "Proxy class to stand in place of ConfigParam defitions in Creator instances, to allow accessing as if a method of the creator class"

    def __init__(self, param_name, param_def, config_def, common_store, creator):
        self.param_name = param_name
        self.param_def = param_def
        self.config_def = config_def
        self.common_store = common_store
        self.creator = creator

    def value(self):
        "Retrieve a parameter from the configuration definition"

        # Get parameter from configuration definition first, then trying common store and
        # if not required return None without error
        if self.param_name in self.config_def:
            param_val = self.config_def[self.param_name]
        elif self.param_name in self.common_store:
            param_val = self.common_store[self.param_name]
        elif self.param_def.required:
            raise KeyError("The parameter name %s requested from config definition by %s was not found and is required" % (self.param_name, self.creator.__class__.__name__))
        else:
            return None

        try:
            return self.param_def.evaluate(param_val, self.creator)
        except ParamError as exc:
            raise ParamError("The parameter named %s requested by creator %s fails with error: %s" % (self.param_name, self.creator.__class__.__name__, exc))

    def __call__(self):
        return self.value()

class Creator(object):
    "Base creator object that handles ensuring the contract of creators parameters is kept and type checking parameters requested"

    def __init__(self, config_def, common_store=None):

        # Create a new common store if non are passed
        if common_store is not None:
            if not isinstance(common_store, dict):
                raise TypeError("The common store passed must be a instance of a dict")
            self.common_store = common_store
        else:
            self.common_store = {}

        # Ensure the config definition is of the appropriate type
        if type(config_def) is dict:
            self.config_def = self._process_config(config_def)
        else:
            raise TypeError("The config definition (config_def) argument passed to %s should be of type dict" % self.__class__.__name__)

        self.parameters = self._load_parameters()

        logger.debug("Initialized creator %s with parameters: %s" % (self.__class__.__name__, self.parameters.keys()))

    def _process_config(self, in_config_def):
        "Process config definition, create nested Creators"

        out_config_def = {}

        for config_name, config_val in in_config_def.items():
            if isinstance(config_val, dict) and "creator" in config_val:
                logger.debug("Initializing nested creator %s for %s" % (config_name, self.__class__.__name__))
                out_config_def[config_name] = config_val["creator"](config_val, common_store=self.common_store)
            else:
                out_config_def[config_name] = config_val

        return out_config_def

    def _load_parameters(self):
        "Gather class parameter types and ensure config_def has necessary items"

        # Look through class attributes to find the ones that have values that
        # inherit from ConfigParam to determine which are the parameters defined
        # in subclasses.
        # Replace in the object the parameter definition with an ParamAccessor proxy
        parameters = {}
        for attr_name in dir(self):
            attr_val = getattr(self, attr_name)
            if isinstance(attr_val, ConfigParam):
                param_proxy = ParameterAccessor(attr_name, attr_val, self.config_def, self.common_store, self)
                parameters[attr_name] = param_proxy
                setattr(self, attr_name, param_proxy)

        return parameters

    def param(self, param_name):
        "Retrieve a parameter from the configuration definition"

        try:
            param_proxy = self.parameters[param_name]
        except KeyError:
            # Creator requested a parameter that was not defined by a ConfigParam attribute
            raise KeyError("Unregistered parameter %s requested by %s" % (param_name, self.__class__.__name__))
        else:
            return param_proxy.value()

    def create(self):
        raise NotImplementedError("Create must be defined in inheriting Creator classes")

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
                    param_proxy = ParameterAccessor(param_name, AnyValue(), self.config_def, self.common_store, self)
                    self.parameters[param_name] = param_proxy
                    setattr(self, param_name, param_proxy)

class ParamPassThru(ParamIterateCreator):
    "Evaluates and passes configurations parameter through as the creator result"

    def create(self):

        result = {} 
        for param_name in self.param_names:
            result[param_name] = self.param(param_name)

        return result

class SaveToCommon(ParamPassThru):
    "Evalualtes parameters and saves them into the common store, creator has no return value"

    def create(self):
        param_vals = super().create()

        for param_name in self.param_names:
            self.common_store[param_name] = param_vals[param_name]
        
        return param_vals
