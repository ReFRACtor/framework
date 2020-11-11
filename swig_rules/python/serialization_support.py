import pickle

def _pickle_no_cxx_load(cls, d):
    inst = cls.__new__(cls)
    inst.__dict__ = d
    return inst
    
class SwigDirectorSerialization:
    '''Helper class for adding serialization to a python director class.
    The issue is that we need to separate out the C++ stuff that gets 
    serialized by boost serialization and the python that gets serialized
    by pickling.'''
    def pickle_no_cxx_save(self):
        d = self.__dict__.copy()
        # TODO Should check for other objects in the dictionary that might
        # be C++ (not sure how). Then move these to a GenericObjectMap in
        # the base class. For now, we just assume that "this" is the only
        # thing we need to remove
        del d['this']
        return pickle.dumps((_pickle_no_cxx_load, (self.__class__, d)))

    def add_this(self, this):
        # Copy syntax from swig. I think this is to support inheriting from
        # multiple directors. We don't have an example of this, but set up
        # the code the same way.
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        

__all__ = ["SwigDirectorSerialization",]
