from .base import ParamIterateCreator
from .. import param

singleton_instances = {}
singleton_captured_objects = {}

def singleton(SingletonCls, instances=singleton_instances, captured_objects=singleton_captured_objects, key=None):
    """Turns any creator into a singleton that only creates an object once, storage can be overridden
    from the default with the storage keyword argument. Instances are stored by the original class's
    name by default unless a different key is supplied.
    """

    storage_key = key or SingletonCls.__name__
    
    # This must be defined inside the decorator function so a fresh version is always returned
    class SingletonCreator(ParamIterateCreator):

        def receive(self, obj):

            obj_list = captured_objects.get(storage_key, [])
            obj_list.append(obj)
            captured_objects[storage_key] = obj_list

        def create(self, **kwargs):

            singleton_value = instances.get(storage_key, None)

            # Defer parameter access until the wrapped creator calls for the item
            # There may be some internal order needed for accessing the parameters
            # that is encoded internally to the wrapped creator
            def param_wrapper(param_name):
                def get_param():
                    return self.param(param_name)
                return get_param

            if singleton_value is None:
                # Collect parameter values and pass to proxied Creator
                param_values = {}
                for param_name in self.param_names:
                    param_values[param_name] = param_wrapper(param_name)

                # Record objects created from the proxied creator for replaying during caching
                self.register_to_receive(object)

                creator = SingletonCls(param_values, common_store=self.common_store)
                singleton_value = creator.create(**kwargs)
                instances[storage_key] = singleton_value

                # Unregister to keep from capturing objects outside the scope of the proxied creator
                self.deregister_to_receive(object)
            else:
                # Emit captured objects from first storage, replaying what happened
                # when the proxied creator was actually created
                obj_list = captured_objects.get(storage_key, [])
                for obj in obj_list:
                    self._dispatch(obj)

            return singleton_value

    # Copy parameter specification from original Creator to our proxied one so they
    # will be registered when it gets 
    for attr_name in dir(SingletonCls):
        attr_val = getattr(SingletonCls, attr_name)
        if isinstance(attr_val, param.ConfigParam):
            setattr(SingletonCreator, attr_name, attr_val)

    return SingletonCreator
