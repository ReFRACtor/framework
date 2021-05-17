#!/usr/bin/env python3

import os
import re
import sys
import logging

import pygraphviz as pgv

from refractor.executor import StrategyExecutor
from refractor.executor.config import find_config_function
from refractor.factory.creator.base import Creator

logger = logging.getLogger(__file__)

CREATOR_COMMON_ATTRIBUTES = {
}

def config_label(cr):
    pass

def create_config_edge(graph, name, cr):
    if not hasattr(cr, "config_def"):
        return []
    
    for sub_name, sub_cr in cr.config_def.items():
        if isinstance(sub_cr, Creator):
            # Connect current creator instance to parent creator instance
            graph.add_edge(name, sub_name)

            # Connect all children Creators
            create_config_edge(graph, sub_name, sub_cr)

            # Create a label that shows the Creator type and
            # includes all parameters and their types
            label = f"{{\\N\ ::\ {sub_cr.__class__.__name__}|"

            for opt_name, opt_obj in sub_cr.parameters.items():
                label += f" | {opt_name}\ ::\ {opt_obj.param_def}"

            label += "}"

            graph.add_node(sub_name, shape="record", label=label, **CREATOR_COMMON_ATTRIBUTES)

def write_config_diagram(exc, output_filename, step_index):

    logger.debug("Loading configuration definition")

    step_keywords = exc.strategy_list[step_index]
    config_func = find_config_function(exc.config_module)
    config_def = config_func(**step_keywords)

    ConfigCreator = config_def['creator']

    c = ConfigCreator(config_def)

    logger.debug("Creating graph")

    graph = pgv.AGraph(directed=True,strict=True)

    create_config_edge(graph, "config", c)

    logger.debug(f"Writing to: {output_filename}")

    graph.layout("dot")
    graph.draw(output_filename)

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Write optical properties set up by configuration")
    
    parser.add_argument("config_filename", metavar="FILE",
        help="File containing a configuration method returning a configuration instance")

    parser.add_argument("--strategy_filename", "-s", metavar="FILE", default=None,
        help="File containing the strategy defining steps to take for the retrival, if not defined a single step will be taken")

    parser.add_argument("--step_index", metavar="INT", type=int, default=0,
        help="Strategy step index to use when a strategy table is defined.")

    parser.add_argument("--output_filename", "-o", metavar="FILE", default=None,
        help="File to write results of configuration execution")

    parser.add_argument("--verbose", "-v", action="store_true", default=False,
        help="Increase level of verbosity")

    args = parser.parse_args()
 
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stdout)

    exc = StrategyExecutor(args.config_filename, 
                           strategy_filename=args.strategy_filename)

    write_config_diagram(exc, args.output_filename, args.step_index)

if __name__ == "__main__":
    main()
