from test_support import *
from refractor.factory.creator.scenario import ScenarioFromL1b
from refractor.framework import *

def test_scenerio_from_l1b():
    l1b = ExampleLevel1b(HdfFile(unit_test_data +
                                 "in/common/l1b_example_data.h5"),
                         "2014090915251774")
    t = ScenarioFromL1b({})
    t.common_store["l1b"] = l1b
    s = t.create()
    expected = ArrayWithUnit_double_1([-25.049287796020508,
                                       -25.048643112182617,
                                       -25.047758102416992], "Degrees")
    assert s["longitude"] == approx(expected)

