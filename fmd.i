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
  %{
    typedef NibbleVector::Iterator RangeVectorIterator;
  %}
#else
  // Use RLEVectors to encode our range endpoint bitmaps
  %rename(RangeVector) RLEVector;
  typedef RLEVector RangeVector;
  %rename(RangeEncoder) RLEEncoder;
  typedef RLEEncoder RangeEncoder;
  %{
    typedef RLEVector::Iterator RangeVectorIterator;
  %}
#endif

// We need to use the inner vector iterator classes to look at vectors. Give a
// partial definition under a new name that works with the conditional typedefs
// above.
class RangeVectorIterator
{
public:
  explicit RangeVectorIterator(const CSA::RangeVector& par);
  ~RangeVectorIterator();

  usint rank(usint value, bool at_least = false);

  usint select(usint index);
  usint selectNext();

  pair_type valueBefore(usint value);
  pair_type valueAfter(usint value);
  pair_type nextValue();

  pair_type selectRun(usint index, usint max_length);
  pair_type selectNextRun(usint max_length);

  bool isSet(usint value);

  usint countRuns();
};

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

// We already worked around the inner classes thing.
#pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS

%include "bits/bitvector.h"
#ifdef USE_NIBBLE_VECTORS
  %include "bits/nibblevector.h"
#else
  %include "bits/rlevector.h"
#endif

%import "rlcsa.h"
%include "fmd.h"

