from string import Formatter

from .base import Creator

import refractor.framework as rf

def as_vector_string(string_vals):
    "Convert a list of strings into a C++ vector of strings"

    vec_str = rf.vector_string()
    for str_val in string_vals:

        if isinstance(str_val, bytes):
            str_val = str_val.decode("UTF-8")
        elif not isinstance(str_val, str):
            raise TypeError("Cannot add incompatible type to vector_string: %s" % str_val)

        vec_str.push_back(str_val)
    return vec_str

class ExtendedFormatter(Formatter):
    """An extended format string formatter

    Formatter with extended conversion symbol
    """
    def convert_field(self, value, conversion):
        """ Extend conversion symbol
        Following additional symbol has been added
        * l: convert to string and low case
        * u: convert to string and up case

        default are:
        * s: convert with str()
        * r: convert with repr()
        * a: convert with ascii()
        """

        if conversion == "u":
            return str(value).upper()
        elif conversion == "l":
            return str(value).lower()
        # Do the default conversion or raise error if no matching conversion found
        super(ExtendedFormatter, self).convert_field(value, conversion)

        # return for None case
        return value

def ObjectCapture(capture_class):
    "Generate a Creator that watches for an object to be emitted elsewhere then stores it internally to use as its return object"

    class ObjectCaptureCreator(Creator):
        def __init__(self, *vargs, **kwargs):
            super().__init__(*vargs, **kwargs)

            self.register_to_receive(capture_class)
            self.captured_object = None

        def receive(self, rec_obj):
            self.captured_object = rec_obj

        def create(self, **kwargs): 
            return self.captured_object

    return ObjectCaptureCreator
