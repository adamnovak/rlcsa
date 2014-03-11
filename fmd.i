// Name the module something Java-y that isn't going to clash with the FMD
// class.
%module FMDUtil

// Set up STL compatibility
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"

%{
#include "fmd.h"
using namespace CSA;
%}

#ifdef MASSIVE_DATA_RLCSA

typedef unsigned long usint;
typedef signed long   sint;

#else

typedef unsigned int  usint;
typedef signed int    sint;

#endif

// Java needs to work with vectors of mappings coming back from the map method.
%template(MappingVector) std::vector<CSA::Mapping>; 

// Java needs to work with pair_types that are locate in text results.
typedef std::pair<usint, usint> pair_type;
%template(pair_type) std::pair<usint, usint>;

// Java can't handle this operator name.
%rename(leftShift) operator<<(std::ostream& o, FMDPosition const& position);

// Whenever any of the JNI classes loads, load the native library.
%pragma(java) jniclasscode=%{
  static {
    RLCSANativeLoader.load();
  }
%}

%import "rlcsa.h"
%include "fmd.h"

