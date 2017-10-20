from refractor import framework as rf

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
