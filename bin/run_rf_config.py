#!/usr/bin/env python3

import sys
import logging

from refractor.executor import StrategyExecutor

def main():

    from argparse import ArgumentParser

    parser = ArgumentParser(description="ReFRACtor retrieval execution")
    
    parser.add_argument("config_filename", metavar="FILE",
        help="File containing a configuration method returning a configuration instance")

    parser.add_argument("config_args", metavar="ARG", nargs="*",
                        help="Additional arguments passed to configuration file")
    parser.add_argument("--strategy_filename", "-s", metavar="FILE", default=None,
        help="File containing the strategy defining steps to take for the retrival, if not defined a single step will be taken")

    parser.add_argument("--output_filename", "-o", metavar="FILE", default=None,
        help="File to write results of configuration execution")

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--retrieval', action='store_true', default=True,
        help="Run a retrieval (default)")
    group.add_argument('--simulation', action='store_true', default=False,
        help="Run a foward model only smulation")

    parser.add_argument("--step_index", metavar="INT", type=int, action="append", default=None,
        help="Strategy step indexes to us instead of all. Specify multiple times for multiple")

    parser.add_argument("--verbose", "-v", action="store_true", default=False,
        help="Increase level of verbosity")

    args = parser.parse_args()

    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stdout)

    exc = StrategyExecutor(args.config_filename, config_args=args.config_args,
                           output_filename=args.output_filename,
                           strategy_filename=args.strategy_filename)

    if args.simulation:
        exc.execute_simulation(args.step_index)
    else:
        exc.execute_retrieval(args.step_index)

if __name__ == "__main__":
    main()
