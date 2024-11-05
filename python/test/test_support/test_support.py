# This contains support routines for unit tests.
import numpy as np
from unittest import SkipTest
from refractor import framework as rf
import os.path
import os
import sys
import subprocess
import filelock
import contextlib
import pytest
    
# Location of test data that is part of source
unit_test_data = os.path.abspath(os.path.dirname(__file__) + "/../../../test/unit/data/") + "/"

# Location of absco test data
absco_data_dir = os.environ.get("ABSCODIR", None)

# Define this environment variable so that it can be used within deserialized data to set paths that get expanded out during test execution
os.environ['abs_top_srcdir'] = os.path.realpath(os.path.join(os.path.dirname(__file__), "../../.."))

# Marker that skips a test if we have a build without boost serialization
# support
require_serialize = pytest.mark.skipif(not rf.have_serialize_supported(),
    reason="need a framework build with boost serialization support to run")

# Short hand for marking as unconditional skipping. Good for tests we
# don't normally run, but might want to comment out for a specific debugging
# reason.
skip = pytest.mark.skip
skipif = pytest.mark.skipif

@pytest.fixture(scope="function")
def isolated_dir(tmpdir):
    '''This is a fixture that creates a temporary directory, and uses this
    while running a unit tests. Useful for tests that write out a test file
    and then try to read it.

    This fixture changes into the temporary directory, and at the end of
    the test it changes back to the current directory.

    Note that this uses the fixture tmpdir, which keeps around the last few
    temporary directories (cleaning up after a fixed number are generated).
    So if a test fails, you can look at the output at the location of tmpdir, 
    e.g. /tmp/pytest-of-smyth
    '''
    curdir = os.getcwd()
    try:
        tmpdir.chdir()
        yield curdir
    finally:
        os.chdir(curdir)


@pytest.fixture(scope='session')
def lock(tmp_path_factory):
    base_temp = tmp_path_factory.getbasetemp()
    lock_file = base_temp.parent / 'serial.lock'
    yield filelock.FileLock(lock_file=str(lock_file))
    with contextlib.suppress(OSError):
        os.remove(path=lock_file)


# Tests with serial all run in serial, needed for our lifetime tests
# because we have a class level variable to keep track of the number
# of objects deleted
@pytest.fixture()
def serial(lock):
    with lock.acquire(poll_interval=0.1):
        yield
        
