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

def creator_node_attrs_simple(creator):
    return {"shape": "box"}

def creator_node_attrs_record(creator):
    "Node attributes to represent a Creator with the Graphviz record shape"

    # Create a label that shows the Creator type and
    # includes all parameters and their types
    label = f"{{\\N\ ::\ {creator.__class__.__name__}|"

    for param_name, param_obj in creator.parameters.items():
        label += f" | {param_name}\ ::\ {param_obj.param_def}"

    label += "}"
    
    return {"shape": "record", "label": label}

def creator_node_attrs_table(creator):
    "Node attributes to represent a Creator with HTML like tables"

    title = f"\\N :: {creator.__class__.__name__}"

    rows = []
    for param_name, param_obj in creator.parameters.items():
        rows.append(f"<TR><TD ALIGN=\"LEFT\">{param_name}</TD><TD ALIGN=\"LEFT\">{param_obj.param_def}</TD></TR>")

    label = """<<TABLE BORDER="1" CELLBORDER="0" CELLPADDING="0">
<TR><TD BORDER="1" COLSPAN="2"><B>{title}</B></TD></TR>
{rows}
</TABLE>>""".format(title=title, rows="\n".join(rows))

    return {"shape": "plaintext", "label": label}


def add_creator_nodes(graph, name, cr, node_attrs_func):
    "Recursively connects edges between creators"

    if not hasattr(cr, "config_def"):
        return
    
    for sub_name, sub_cr in cr.config_def.items():
        if isinstance(sub_cr, Creator):
            # Connect current creator instance to parent creator instance
            graph.add_edge(name, sub_name)

            # Connect all children Creators
            add_creator_nodes(graph, sub_name, sub_cr, node_attrs_func)

            # Add attributes to the node to control it's shape
            graph.add_node(sub_name, **node_attrs_func(sub_cr))

def add_parameter_edges(graph, name, cr):
    "Add edges to any existing nodes from Creators"

    if not hasattr(cr, "config_def"):
        return

    for param_name, param_obj in cr.parameters.items():

        if param_name in graph.nodes() and not param_name in cr.config_def.keys():
            graph.add_edge(name, param_name)
    
    for sub_name, sub_cr in cr.config_def.items():
        if isinstance(sub_cr, Creator):
            add_parameter_edges(graph, sub_name, sub_cr)

def write_config_diagram(exc, output_filename, step_index, node_attrs_func=creator_node_attrs_simple, param_edges=False):

    logger.debug("Loading configuration definition")

    # Load configuration definition
    step_keywords = exc.strategy_list[step_index]
    config_func = find_config_function(exc.config_module)
    config_def = config_func(**step_keywords)

    # Instead of running process_config we use the same creator type defined at the 
    # root level so contain the configuration
    ConfigCreator = config_def['creator']
    root_creator = ConfigCreator(config_def)

    logger.debug("Creating graph")

    # Recursively descend the configuration gathering edges between Creators
    graph = pgv.AGraph(directed=True,strict=True)

    logger.debug("Adding creator nodes and edges")
    add_creator_nodes(graph, "config", root_creator, node_attrs_func)

    if param_edges:
        logger.debug("Adding parameter edges")
        add_parameter_edges(graph, "config", root_creator)

    logger.debug(f"Writing to: {output_filename}")

    # The output filename will determine the file format
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

    parser.add_argument("--shape", default="simple", choices=["simple", "record", "table"],
        help="Controls how the graph nodes are displayed. Default = simple")

    parser.add_argument("--param_edges", default=False, action="store_true",
        help="Add edges from creators for parameters that are named but not present in config definition")

    parser.add_argument("--verbose", "-v", action="store_true", default=False,
        help="Increase level of verbosity")

    args = parser.parse_args()
 
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    if args.shape == "record":
        logger.debug("Using record method for node shapes")
        node_attrs_func = creator_node_attrs_record
    elif args.shape == "table":
        logger.debug("Using table method for node shapes")
        node_attrs_func = creator_node_attrs_table
    else:
        logger.debug("Using simple method for node shapes")
        node_attrs_func = creator_node_attrs_simple

    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stdout)

    exc = StrategyExecutor(args.config_filename, 
                           strategy_filename=args.strategy_filename)

    write_config_diagram(exc, args.output_filename, args.step_index, node_attrs_func=node_attrs_func, param_edges=args.param_edges)

if __name__ == "__main__":
    main()
