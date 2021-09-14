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
