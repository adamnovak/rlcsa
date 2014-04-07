// Name the module something Java-y that isn't going to clash with the RLCSA
// class.
%module RLCSAUtil

// Set up STL compatibility
%include "std_string.i"
%include "std_pair.i"

// And pointer arrays
%include "carrays.i"

// And pointers
%include "cpointer.i"

// And integers of various sizes
%include <stdint.i>


%{
#include "rlcsa.h"
using namespace CSA;
%}

#ifdef MASSIVE_DATA_RLCSA

typedef unsigned long usint;
typedef signed long   sint;

#else

typedef unsigned int  usint;
typedef signed int    sint;

#endif

// uchars are unsigned chars.
typedef unsigned char uchar;

// Java needs to work with pair_types that are count results.
typedef std::pair<usint, usint> pair_type;
%template(pair_type) std::pair<usint, usint>;

// We need to be able to work with arrays of usint
%array_functions(usint,USIntArray);

// And with pointers to usint
%pointer_functions(usint,USIntPointer);

// And with arrays of characters (instead of null-terminated strings) for the
// display output.
%array_functions(uchar,UCharArray)

// Whenever any of the JNI classes loads, load the native library.
%pragma(java) jniclasscode=%{
  static {
    RLCSANativeLoader.load();
  }
%}

// Mark some of the display methods as making new objects we need to manage the
// memory of. TODO: get them all.
%newobject RLCSA::display(usint sequence, bool include_end_marker = false);
%newobject RLCSA::display(usint sequence, pair_type range) const;

uchar* RLCSA::display(usint sequence, pair_type range) const;

// Those display methods don't null-terminate their strings. Make sure they are
// arrays by not saying those uchar*s should be strings.

// The other uchar*s are actually char*s (null-terminated strings)
%apply char* { uchar* pattern, uchar* base_name };


%include "rlcsa.h"

// We also take the important inlines in misc/definitions.h for working with
// pair_type ranges.
%include "misc/definitions.h"
