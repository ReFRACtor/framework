# This contains a few useful function to extend GenericObject from swig

from refractor.framework_swig.generic_object import GenericObject

def _subobject(self, recursive=True):
    '''Return all the (possibly empty) subobjects found for this object.
    If recursive is True, then walk through each subobject recursively 
    returning all of its subobjects.'''
    if(hasattr(self, "subobject_list")):
        for obj in self.subobject_list:
            obj_sp = GenericObject.convert_to_most_specific_class(obj)
            if(obj_sp is not None):
                yield obj_sp
                if(recursive):
                    yield from obj_sp.subobject(recursive=True)

GenericObject.subobject = _subobject

def _find_all_subobject_of_type(self, typ):
    '''Find all subobjects that are an instance of the given type.
    We remove duplicates'''
    found = list()
    for obj in self.subobject():
        if(isinstance(obj, typ) and
           not next((x for x in found if x.is_same_ptr(obj)), False)):
            found.append(obj)
            yield obj

GenericObject.find_all_subobject_of_type = _find_all_subobject_of_type

def _find_subobject_of_type(self, typ):
    '''Return the first subobject we find of the given type, or None if
    we don't find any'''
    return next(self.find_all_subobject_of_type(typ), None)

GenericObject.find_subobject_of_type = _find_subobject_of_type

__all__ = []
