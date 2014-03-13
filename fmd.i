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

// Make sure to rename these things that we swap out with #defines to a
// consistent name. Also note that these "vectors" are bit vectors, not
// std::vectors.
#ifdef USE_NIBBLE_VECTORS
  // Use Nibble Vectors to encode our range endpoint bitmaps
  %rename(RangeVector) NibbleVector;
  typedef NibbleVector RangeVector;
  %rename(RangeEncoder) NibbleEncoder;
  typedef NibbleEncoder RangeEncoder;
#else
  // Use RLEVectors to encode our range endpoint bitmaps
  %rename(RangeVector) RLEVector;
  typedef RLEVector RangeVector;
  %rename(RangeEncoder) RLEEncoder;
  typedef RLEEncoder RangeEncoder;
#endif

// Java needs to work with vectors of mappings coming back from the map method.
%template(MappingVector) std::vector<CSA::Mapping>; 

// Java also needs to work with vectors of sints coming back from the map method
// when working on ranges.
%template(SintVector) std::vector<sint>; 

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

// Skip inner classes we can't wrap anyway (so no way to read the bit vectors.
// Sorry.) We'd %ignore them, but Swig ignores that.
#pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS

%include "bits/bitvector.h"
#ifdef USE_NIBBLE_VECTORS
  %include "bits/nibblevector.h"
#else
  %include "bits/rlevector.h"
#endif

%import "rlcsa.h"
%include "fmd.h"

