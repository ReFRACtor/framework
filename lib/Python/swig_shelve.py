import sys
# The Shelve package is a very useful way to easily implement persistence.
# But it has the disadvantage that depending on the system we are on different
# databases will be available (e.g., not every system has berkely db).
# Also the Shelve package is apparently deprecated in future version of python. 
# I spent some time googling, and several places recommended using sqlite3
# to implement this. In addition to being available (we check this at as 
# part of the configure script) it is also faster.
#
# This class gives a Shelve like interface to sqlite.
#
# In addition, we give and xml and json interface also

# This is python 3 only, but our code already requires that so no problem

from collections.abc import MutableMapping as DictMixin

import pickle
import sqlite3
import os.path
import refractor_swig

# Allow jsonpickle to be used, if available.
try:
    import jsonpickle
    have_jsonpickle = True
    jsonpickle.set_encoder_options('json', sort_keys=True, 
                                   indent=4, separators=(',', ': '))
except ImportError:
    have_jsonpickle = False

def to_db_type(value):
    """If value's type is supported natively in SQLite, return value.
    Otherwise, return a pickled representation. """
    if value is None or isinstance(value, (int, int, float, bytes, str)):
        return value
    else:
        # Note the types really are different for python 2 vs. 3. This is 
        # because bytes isn't really different than str in python 2, but is 
        # in 3
        if sys.version_info > (3,):
            return bytes(pickle.dumps(value))
        else:
            return buffer(pickle.dumps(value))

def from_db_type(value):
    """ Converts a value from the database to a Python object. """
    if sys.version_info > (3,):
        if isinstance(value, bytes):
            return pickle.loads(value)
    else:
        if isinstance(value, buffer):
            return pickle.loads(value)
    return value

def read_shelve(f):
    '''This handles reading a value from a shelve file. The string f should
    be of the form file_name:key, file_name.xml, file_name.xml.gz,
    file_name.bin, file_name.bin.gz or file_name.json. We 
    open the given file, and read the value for the given key.

    A problem with the python shelve/pickle files is that it can't
    communicate directly with a C++ program, and also the files aren't
    human readable or portable.  So we also support xml files, we key
    off of the file name and if it is something like "foo.xml" we read
    that rather than a shelve file. "foo.bin" is the binary version of
    of the serialized data.

    We also support files with the extension "foo.json" as json pickle
    files. This doesn't have any additional functionality over the
    sqlite python, but it does make for more human readable files that
    are sometimes preferable. We use jsonpickle for reading these
    files.

    Note that it can be useful to execute python code before using
    a shelve file, e.g., we are using python modules not already 
    included in AFIDS. If the special key "_extra_python_init" is 
    found, we execute the code found there. This can do things like 
    import modules. For xml, bin and json files, we look for the file 
    "extra_python_init.py" found in the same directory.

    Because we often use relative names for files, we first chdir to 
    the same directory as the database file (if different than the current
    one). We change back to the current directory when done.

    '''
    fname = f.split(':')[0]
    dirn, fb = os.path.split(fname)
    curdir = os.getcwd()
    try:
        if(dirn):
            os.chdir(dirn)
        f2, ext = os.path.splitext(f)
        if(ext not in (".gz")):
            f2 = f
        if(os.path.splitext(f2)[1] == ".xml"):
            if(os.path.exists("extra_python_init.py")):
                exec(open("extra_python_init.py").read())
            return refractor_swig.serialize_read_generic(fb)
        if(os.path.splitext(f2)[1] == ".bin"):
            if(os.path.exists("extra_python_init.py")):
                exec(open("extra_python_init.py").read())
            return refractor_swig.serialize_read_binary_generic(fb)
        if(os.path.splitext(f)[1] == ".json"):
            if(os.path.exists("extra_python_init.py")):
                exec(open("extra_python_init.py").read())
            if(have_jsonpickle):
                return jsonpickle.decode(open(fb).read())
            else:
                raise RuntimeError("Use of .json file requires jsonpickle package to be installed")
        t = SQLiteShelf(fb, "r")
        if("_extra_python_init" in list(t.keys())):
            exec(t["_extra_python_init"])
        key = f.split(':')[1]
        return t[key]
    finally:
        os.chdir(curdir)

def shelve_time_after(f1, f2):
    '''Compare the update time on 2 shelve objects, return if f1 update time
    >= f2 update time. Note that either f1 or f2 can be files, in which case
    we use the file modify time instead.
    It is ok if f1 doesn't exist, in that case always return False. '''
    if(':' in f1):
        fname, key = f1.split(':')
        t = SQLiteShelf(fname, "r")
        if key in t:
            f1time = t.update_time_unix(key)
        else:
            return False
    else:
        if(os.path.exists(f1)):
            f1time = os.path.getmtime(f1)
        else:
            return False
    if(':' in f2):
        fname, key = f2.split(':')
        t = SQLiteShelf(fname, "r")
        f2time = t.update_time_unix(key)
    else:
        f2time = os.path.getmtime(f2)
    return f1time >= f2time
        

def write_shelve(f, val):
    '''This handles writing a value to a shelve file, possibly creating the
    file is it doesn't exist. The string f should be of the form
    file_name:key. We open/create the given file and write the value for
    the given key.

    A problem with the shelve files is that it can't communicate directly 
    with a C++ program, and also the files aren't human readable or portable.
    So we also support xml and bin files, we key off of the file name and if it 
    is something like "foo.xml" we write that rather than a shelve file. You
    can add a ".gz" to use gzip compression on the file.
    '''
    # Strip off any compression part, when determining file type
    f2, ext = os.path.splitext(f)
    if(ext not in (".gz")):
        f2 = f
    if(os.path.splitext(f2)[1] == ".xml"):
        refractor_swig.serialize_write(f, val)
        return
    if(os.path.splitext(f2)[1] == ".bin"):
        refractor_swig.serialize_write_binary(f, val)
        return
    if(os.path.splitext(f)[1] == ".json"):
        if(have_jsonpickle):
            with open(f, "w") as fh:
                fh.write(jsonpickle.encode(val))
            return 
        else:
            raise RuntimeError("Use of .json file requires jsonpickle package to be installed")
    fname, key = f.split(':')
    d = SQLiteShelf(fname)
    d[key] = val
    d.close()

class SQLiteShelf(DictMixin):
    """Shelf implementation using an SQLite3 database. """
    def __init__(self, filename, mode = "r+"):
        '''Open an existing file, or create a new one if it doesn't exist.
        The mode can be "r+" for read/write or "r" for read only.
        '''
        self._database = sqlite3.connect(filename)
        self._database.execute("CREATE TABLE IF NOT EXISTS Shelf "
                               "(Key TEXT PRIMARY KEY NOT NULL, Value BLOB, Updated datetime)")
        self._open = True
        self._read_only = (mode == "r")

    def __del__(self):
        self.close()

    def __getitem__(self, key):
        row = self._database.execute("SELECT Value FROM Shelf WHERE Key=?",
                                     [key]).fetchone()
        if row:
            return from_db_type(row[0])
        else:
            raise KeyError(key)

    def update_time(self, key):
        '''Return updated time as a string.'''
        row = self._database.execute("SELECT Updated FROM Shelf WHERE Key=?",
                                     [key]).fetchone()
        if row:
            return from_db_type(row[0])
        else:
            raise KeyError(key)

    def update_time_julian(self, key):
        '''Return updated time as Julian day, including fraction'''
        row = self._database.execute("SELECT julianday(Updated) FROM Shelf WHERE Key=?",
                                     [key]).fetchone()
        if row:
            return from_db_type(row[0])
        else:
            raise KeyError(key)

    def update_time_unix(self, key):
        '''Return updated time as unix time, including fraction'''
        # Unix epoch in Julian days is 2440587.5
        return (self.update_time_julian(key) - 2440587.5) * 86400.0

    def touch(self, key):
        '''Change the updated time to now (like the unix command "touch")'''
        if(self._read_only):
            raise RuntimeError("Attempt to write to read only shelve.")
        self._database.execute("UPDATE Shelf SET Updated=strftime('%Y-%m-%d %H:%M:%f', 'now') WHERE Key=?", [key])

    def __setitem__(self, key, value):
        if(self._read_only):
            raise RuntimeError("Attempt to write to read only shelve.")
        self._database.execute("INSERT OR REPLACE INTO Shelf VALUES (?, ?, strftime('%Y-%m-%d %H:%M:%f', 'now'))",
                               [key, to_db_type(value)])

    def __delitem__(self, key):
        if(self._read_only):
            raise RuntimeError("Attempt to delete from read only shelve.")
        self._database.execute("DELETE FROM Shelf WHERE Key=?", [key])

    def keys(self):
        """Return a list of keys in the shelf."""
        return [row[0] for row in
                self._database.execute("SELECT Key FROM Shelf")]

    def close(self):
        """Commit changes and close the file."""
        if self._database is not None:
            self._database.commit()
            self._database.close()
            self._database = None 

    # These are needed by python 3, but not python 2
    def __len__(self):
        raise RuntimeError("Not implemented yet")

    def __iter__(self):
        raise RuntimeError("Not implemented yet")

__all__ = ["read_shelve", "shelve_time_after", "write_shelve", "SQLiteShelf",
           "have_jsonpickle"]
           
