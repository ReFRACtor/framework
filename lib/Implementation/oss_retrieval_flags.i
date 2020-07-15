%include "fp_common.i"
%{
#include "oss_retrieval_flags.h"
%}

%base_import(generic_object)

%fp_shared_ptr(FullPhysics::OssRetrievalFlags);

%include "oss_retrieval_flags.h"