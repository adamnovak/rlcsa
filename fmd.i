// Name the module something Java-y that isn't going to clash with the FMD
// class.
%module FMDUtil

// Set up STL compatibility
%include "std_string.i"
%include "std_vector.i"

%{
#include "fmd.h"
using namespace CSA;
%}

// Java needs to work with vectors of mappings coming back from the map method.
%template(MappingVector) std::vector<Mapping>; 

// Java can't handle this operator name.
%rename(leftShift) operator<<(std::ostream& o, FMDPosition const& position);

%import "rlcsa.h"
%include "fmd.h"

