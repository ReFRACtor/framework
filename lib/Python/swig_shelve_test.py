from .swig_shelve import *
import time
from .test_support import *
import pickle
from refractor_swig import *

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
    assert psigma.a == approx(psigma2.a)
    assert psigma.b == approx(psigma2.b)

@require_serialize
def test_all_formats(isolated_dir):
    psigma = rf.PressureSigma([0,0,0],[0.3,0.6,1.0], 10, True)
    write_shelve("psigma.xml", psigma)
    # Need to fix this
    # write_shelve("psigma.bin", psigma)
    write_shelve("psigma.json", psigma)
    write_shelve("data.db:psigma", psigma)
    pickle.dump(psigma, open( "psigma.pkl", "wb" ))

    psigma1 = read_shelve("psigma.xml")
    # Need to fix this
    # psigma2 = read_shelve("psigma.bin")
    psigma3 = read_shelve("psigma.json")
    psigma4 = read_shelve("data.db:psigma")
    psigma5 = pickle.load(open("psigma.pkl", "rb"))
    assert psigma.a == approx(psigma1.a)
    assert psigma.b == approx(psigma1.b)
    assert psigma.a == approx(psigma3.a)
    assert psigma.b == approx(psigma3.b)
    assert psigma.a == approx(psigma4.a)
    assert psigma.b == approx(psigma4.b)
    assert psigma.a == approx(psigma5.a)
    assert psigma.b == approx(psigma5.b)

@require_serialize
def test_generic_object_map(isolated_dir, sample_aerosol):
    a = sample_aerosol
    p = a.pressure
    write_shelve("a.xml", a)
    write_shelve("p.xml", p)
    a2 = read_shelve("a.xml")
    p2 = read_shelve("p.xml")
    assert a2.pressure.surface_pressure == approx(p2.surface_pressure)
    # Updating p2 won't change pressure in a2
    p2.surface_pressure = AutoDerivativeDouble(100)
    assert p2.surface_pressure.value == approx(100)
    assert a2.pressure.surface_pressure != approx(p2.surface_pressure)
    # Repeat, using GenericObjectMap
    m = GenericObjectMap()
    m["a"] = a
    m["p"] = p
    write_shelve("m.xml", m)
    m3 = read_shelve("m.xml")
    a3 = m3.a
    p3 = m3.p
    assert a3.pressure.surface_pressure == approx(p3.surface_pressure)
    # Updating p3 will change pressure in a3
    p3.surface_pressure = AutoDerivativeDouble(100)
    assert p3.surface_pressure.value == approx(100)
    assert a3.pressure.surface_pressure == approx(p3.surface_pressure)
    
