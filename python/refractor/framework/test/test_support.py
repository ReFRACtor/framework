import os
import logging
from functools import wraps

from refractor.framework import load_config_module, find_config_function, NoConfigFoundError

logger = logging.getLogger(__name__)

def dummy_env_var(var_name, dummy_val=""):
    """Set up dummy paths for environment variables only for the duration of a function"""

    def decorator(func):

        @wraps(func)
        def wrapper(*vargs, **kwargs):

            # Only define a dummy value if environment variable is not already present
            var_defined = True
            if not var_name in os.environ:
                os.environ[var_name] = dummy_val
                var_defined = False

            func(*vargs, **kwargs)

            # Delete if it was not already defined to put state back to what it was
            if not var_defined:
                del os.environ[var_name]
        return wrapper
    return decorator

def check_config_loading(config_filenames, config_function_args=None):
    """Loads a list of ReFRACtor configuration files iteratively to check for errors. This is
    intended to be used in unit test to look for import errors in config files.

    config_function_args is a dictionary that can define values to be supplied to configuration
    functions that require arguments.
    """

    if config_function_args is None:
        config_function_args = {}

    for config_fn in config_filenames:
        logger.debug(f"Trying to load {config_fn}")
        config_module = load_config_module(config_fn)
        try:
            config_func = find_config_function(config_module)
            args = config_function_args.get(config_func.__name__, ())
            logger.debug(f"Calling config function {config_func.__name__} with {len(args)} arguments")
            config_def = config_func(*args)
        except NoConfigFoundError as exc:
            logger.debug(f"No config found in {config_fn}, skipping")
