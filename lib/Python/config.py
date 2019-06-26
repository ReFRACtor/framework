import os
import sys
import string
import random

import importlib.util

CONFIG_MARKER = "__refractor_config"
STRATEGY_MARKER = "__refractor_strategy"

def refractor_config(func):
    "Decorator to mark a function as a returning a ReFRACtor config"

    setattr(func, CONFIG_MARKER, True)

    return func

def strategy_list(func):
    "Decorator to mark a function as a returning a strategy list"

    setattr(func, STRATEGY_MARKER, True)

    return func

def load_config_module(filename):
    "Loads in a Python module by filename in a manner expected of ReFRACtor config function containing modules"

    # Add some random characters to the module name to make each read instance unique
    uniq_module_name = "config_" + ''.join([ random.choice(string.ascii_lowercase) for i in range(5) ])

    # Add the path where the file is found to the system path list so imports for
    # other modules in the same directory can be found
    # There is probably a better way to do this that does not pollute the systemwide
    # path variable as a side effect
    file_path = os.path.realpath(os.path.dirname(filename))

    if not file_path in sys.path:
        sys.path.append(file_path)

    # Import the file
    spec = importlib.util.spec_from_file_location(uniq_module_name, filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module

def find_marked_uniq_function(module, marker):
    "Finds a function marked with a certain attribute in a named module and returns the only unique one. An error occurs if more than one marked function is present in a module."

    # Search for a function tagged with __refractor_config by the refractor_config decorator above
    # Throw an error if multiple are found
    config_func = None
    config_found = False
    for mod_item_name in dir(module):
        mod_item = getattr(module, mod_item_name)
        if getattr(mod_item, marker, False) and getattr(mod_item, "__module__", None) == module.__name__:
            if not config_found:
                config_func = mod_item
                config_found = True
            else:
                raise Exception("Only one ReFRACtor configuration function should be present in file: {}".format(filename))
    
    if not config_found:
        raise Exception("No ReFRACtor configuration found in file: {}".format(filename))

    return config_func

def find_config_function(module):
    "Finds the item in a ReFRACtor configuration tagged as a configuration function"
    
    return find_marked_uniq_function(module, CONFIG_MARKER)

def find_strategy_function(module):
    "Finds the item in a ReFRACtor configuration tagged as a strategy function"
    
    return find_marked_uniq_function(module, STRATEGY_MARKER)
