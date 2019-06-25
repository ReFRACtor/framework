from .swig_shelve import *
import time
from .test_support import *

def test_sqlite_shelf(isolated_dir):
    d = SQLiteShelf("sqlite_shelf.db")
    d["value1"] = 1
    time.sleep(0.1)
    d["value2"] = [1, 2, 3, "blah"]
    d["value3"] = 3
    d["value3"] = 4
    d["value4"] = 5
    del d["value4"]
    d.close()
    d = SQLiteShelf("sqlite_shelf.db", "r")
    assert d["value1"] == 1
    assert d["value2"] == [1, 2, 3, "blah"]
    assert d["value3"] == 4
    with pytest.raises(KeyError) as e_info:
        d["value4"]
    with pytest.raises(RuntimeError) as e_info:
        d["value4"] = 1
    # This should be close to 0.1 seconds. But don't depend on that for
    # a test, instead just say we are ordered
    #print d.update_time_unix("value2") - d.update_time_unix("value1")
    assert d.update_time_julian("value2") > d.update_time_julian("value1")
    assert shelve_time_after("sqlite_shelf.db:value2",
                             "sqlite_shelf.db:value1")
    time.sleep(0.1)
    d.close()
    d = SQLiteShelf("sqlite_shelf.db", "r+")
    d.touch("value1")
    d.close()
    d = SQLiteShelf("sqlite_shelf.db", "r")
    assert d.update_time_julian("value1") > d.update_time_julian("value2")

def test_read_write_shelf(isolated_dir):
    write_shelve("sqlite_shelf.db:value1", 1)
    write_shelve("sqlite_shelf.db:value1", 3)
    write_shelve("sqlite_shelf.db:value2", [1, 2, 3, "blah"])
    assert read_shelve("sqlite_shelf.db:value1") == 3
    assert read_shelve("sqlite_shelf.db:value2") == [1, 2, 3, "blah"]

@require_serialize    
def test_read_write_xml(isolated_dir):
    psigma = rf.PressureSigma([0,0,0],[0.3,0.6,1.0], 10, True)
    write_shelve("sqlite_shelf_test.xml", psigma)
    psigma2 = read_shelve("sqlite_shelf_test.xml")
    assert_array_almost_equal(psigma.a, psigma2.a)
    assert_array_almost_equal(psigma.b, psigma2.b)

