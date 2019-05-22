#!/usr/bin/env python3

import sys
import logging

from refractor.executor import RetrievalStrategyExecutor

def main():

    from argparse import ArgumentParser

    parser = ArgumentParser(description="ReFRACtor retrieval execution")
    
    parser.add_argument("--config_filename", "-c", metavar="FILE", required=True,
        help="File containing a configuration method returning a configuration instance")

    parser.add_argument("--strategy_filename", "-s", metavar="FILE", default=None,
        help="File containing the strategy defining steps to take for the retrival, if not defined a single step will be taken")

    parser.add_argument("--verbose", "-v", action="store_true", default=False,
        help="Increase level of verbosity")

    args = parser.parse_args()

    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stdout)

    ret = RetrievalStrategyExecutor(args.config_filename, args.strategy_filename)

    ret.execute_retrieval()

if __name__ == "__main__":
    main()
