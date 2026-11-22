import subprocess
import sys

from test_support import *


def run(expression):
    '''Evaluate an expression in a separate interpreter: the failure this tests for is a
    crash of the interpreter, which would end the test run.'''
    code = ("import refractor.framework as rf\n"
            "u = rf.Unit('m')\n"
            "a = rf.AutoDerivativeDouble(1.0)\n"
            f"print({expression})")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, f"{expression} ended with status {r.returncode}: {r.stderr[-300:]}"
    return r.stdout.strip()


@pytest.mark.parametrize("expression, expected", [
    ("u == None", "False"),
    ("u != None", "True"),
    ("None == u", "False"),
    ("None != u", "True"),
    ("a == None", "False"),
    ("a != None", "True"),
    ("None == a", "False"),
    ("u == 'm'", "False"),
    ("u == rf.Unit('m')", "True"),
    ("u != rf.Unit('cm')", "True"),
    ("a == 1.0", "True"),
    ("a == 2.0", "False"),
    ("a == rf.AutoDerivativeDouble(1.0)", "True"),
    ("a != rf.AutoDerivativeDouble(2.0)", "True"),
])
def test_comparison_does_not_crash(expression, expected):
    '''Comparing a Unit or an AutoDerivative with None used to crash the interpreter.'''
    assert run(expression) == expected
