%include "fp_common.i"
%{
#include "first_order_driver.h"
%}

%import "first_order_interface.i"

%base_import(spurr_driver)

%fp_shared_ptr(FullPhysics::FirstOrderDriver);

%include "first_order_driver.h"
