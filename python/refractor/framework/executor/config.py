import os
import sys
import string
import random

import importlib.util

from ..config import CONFIG_MARKER, STRATEGY_MARKER

def load_config_module(filename):
    '''Loads in a Python module by filename in a manner expected of 
    ReFRACtor config function containing modules.

    If filename is a directory, we assume we are actually loading __init__.py
    found in that directory'''

    # Fully expand filename to avoid problems introduced by relative paths
    filename = os.path.realpath(filename)

    # Add some random characters to the module name to make each read instance unique
    uniq_module_name = "config_" + ''.join([ random.choice(string.ascii_lowercase) for i in range(5) ])

    # This may need some cleanup. What we are trying to do here is to allow
    # a config file to load other files in the same directory as it is found,
    # without conflicting with other configurations (e.g., both omps_nm and
    # cris have a "base_config.py" file). So we treat the directory as its
    # own module, to keep these separate. We load the directory as a module
    # if there is a __init__.py file. The config file then does relative
    # loading of other files, each "from .base_config import blah"
    #
    # Not really sure about the logic here, I've always found python 
    # modules/search paths really arcane
    #
    if os.path.isdir(filename):
        # If a directory, import the module
        spec = importlib.util.spec_from_file_location(uniq_module_name, filename + "/__init__.py")
        module = importlib.util.module_from_spec(spec)
    else:
        spec = importlib.util.spec_from_file_location(uniq_module_name, filename)
        # If we find a __init__.py, treat the directory as a module that
        # is the scope of filename
        init_filename = os.path.dirname(filename) + "/__init__.py"
        if(os.path.isfile(init_filename)):
            spec2 = importlib.util.spec_from_file_location(uniq_module_name, os.path.dirname(filename) + "/__init__.py")
            module = importlib.util.module_from_spec(spec2)
        else:
            module = importlib.util.module_from_spec(spec)

    sys.modules[spec.name] = module    
    spec.loader.exec_module(module)

    return module

class NoConfigFoundError(Exception):
    pass

class MultipleConfigFoundError(Exception):
    pass

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
                raise MultipleConfigFoundError("Only one ReFRACtor configuration function should be present in file: {}".format(module.__file__))
    
    if not config_found:
        raise NoConfigFoundError("No ReFRACtor configuration found in file: {}".format(module.__file__))

    return config_func

def find_config_function(module):
    "Finds the item in a ReFRACtor configuration tagged as a configuration function"
    
    return find_marked_uniq_function(module, CONFIG_MARKER)

def find_strategy_function(module):
    "Finds the item in a ReFRACtor configuration tagged as a strategy function"
    
    return find_marked_uniq_function(module, STRATEGY_MARKER)
