// Name the module something Java-y that isn't going to clash with the RLCSA
// class.
%module RLCSAUtil

// Set up STL compatibility
%include "std_string.i"
%include "std_pair.i"

%{
#include "rlcsa.h"
using namespace CSA;
%}

// Java needs to work with pair_types that are count results.
%template(PairType) std::pair<usint, usint>; 

%include "rlcsa.h"

