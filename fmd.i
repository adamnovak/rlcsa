%module fmd
%include "std_string.i"
%include "std_vector.i"

%{
#include "fmd.h"
%}

%template(MappingVector) std::vector<Mapping>; 

%include "fmd.h"
