#!/usr/bin/env python3

import os
import sys
import logging

from refractor import framework as rf
from refractor.executor.config import load_config_module, find_config_function
from refractor.factory import process_config
from refractor import write_shelve

logging.basicConfig(level=logging.DEBUG)

config_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../config"))
input_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../data/in"))

config_fn = os.path.join(config_dir, "uv_atmosphere.py")
serialized_fn = os.path.join(input_dir, "uv_atmosphere/uv_atmosphere.bin.gz")

logging.info(f"Loading configuration: {config_fn}")

config_module = load_config_module(config_fn)
config_func = find_config_function(config_module)

config_def = config_func()

config_inst = process_config(config_def)

logging.info(f"Writing serialized representation: {serialized_fn}")
write_shelve(serialized_fn, config_inst.atmosphere)
