#ifndef FMD_H
#define FMD_H

/**
 * Define a macro for easily compiling in/out detailed debugging information
 * about FMD search.
 */
#define DEBUG(op) op
//#define DEBUG(op)

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "bits/deltavector.h"
#include "bits/rlevector.h"
#include "bits/nibblevector.h"
#include "bits/succinctvector.h"

#include "sasamples.h"
#include "alphabet.h"
#include "lcpsamples.h"
#include "misc/parameters.h"
#include "sampler.h"
#include "suffixarray.h"
#include "rlcsa.h"
#include "misc/definitions.h"

namespace CSA
{

const static usint NUM_BASES = 5;

// This holds the bases in alphabetcal order by reverse complement. The only
// time the order of the bases matters is when doing the iterative scoping out
// of the reverse complement intervals in the extension procedure, and there we
// need to go through them in this order.
const static std::string BASES = "TGCNA";


/**
 * Return true if a character is a valid DNA base, and false otherwise. Only
 * capital letters are allowed, and N counts.
 */
inline bool isBase(usint input)
{
  for(std::string::const_iterator i = BASES.begin(); i != BASES.end(); ++i)
  {
    if((char)input == *i)
    {
      return true;
    }
  }
  return false;
}


/**
 * Return the "reverse" complement of a single character. Only capital letters
 * are allowed, and N is its own reverse complement.
 */
inline usint reverse_complement(usint input) {
  switch((char)input)
  {
    case 'A':
      return (usint)'T';
    case 'C':
      return (usint)'G';
    case 'G':
      return (usint)'C';
    case 'T':
      return (usint)'A';
    case 'N':
      return (usint)'N';
    default:
      throw "Invalid character to reverse complement";
  }
}

/**
 * Represents the state (or result) of an FMD-index search, which is two ranges
 * (one for the forward sequence, and one for the reverse complement) of equal
 * length. The ranges are stored as two start indices and a length. They can be
 * in either SA space (not counting the text start symbols at the beginning of
 * the BWT) or in BWT space.
 *
 * Range semantics are inclusive, so a length = 0 range holds 1 thing and its
 * reverse complement.
 */
struct FMDPosition
{
  usint forward_start;
  usint reverse_start;
  // Length 0 = only the entry at start/end
  sint length;
  FMDPosition();
  FMDPosition(usint forward_start, usint reverse_start, usint length);
  /** 
   * Flip the FMDPosition around so the reverse complement interval is the
   * forward interval and visa versa.
   */
  FMDPosition flip() const;
};

/**
 * Provide pretty-printing for FMDPositions. See
 * <http://www.parashift.com/c++-faq/output-operator.html>
 */
std::ostream& operator<< (std::ostream& o, FMDPosition const& position);

const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, -1);

/**
 * Is an FMDPosition empty?
 */
inline bool isEmpty(const FMDPosition& position)
{
  return position.length < 0;
}

/**
 * Return the actual number of matches represented by an FMDPosition.
 */
inline usint length(const FMDPosition& position)
{
  return position.length + 1;
}

/**
 * Represents a mapping between a base in a string and a (text, index) position
 * in the FMD-index. Contains the text and offset to which a character maps, and
 * a flag to say if it represents a real mapping or a result of "unmapped".
 */
struct Mapping
{
  // Holds (text, position)
  pair_type location;
  bool is_mapped;
  Mapping();
  Mapping(pair_type location, bool is_mapped=true);
};

/**
 * A triple to hold the return values from FMD::mapPosition(). Holds a flag for
 * whether the mapping succeeded or not, an FMDPosition corresponding either to
 * where the character mapped or the longest search starting at the character
 * that did actually return results, and the number of characters in the
 * FMDPosition's search pattern.
 */
struct MapAttemptResult
{
  bool is_mapped;
  FMDPosition position;
  uint characters;
};

/**
 * Defines an RLCSA index derivative that represents an FMD-index: an index of
 * DNA sequences (over the alphabet {A, C, G, T, N}) where all texts are present
 * with their reverse complements.
 *
 * In such an index, an ongoing search can be extended or retracted at either
 * end in O(1) time.
 *
 * See the paper "Exploring single-sample SNP and INDEL calling with whole-
 * genome de novo assembly" (2012), by Heng Li, which defines the FMD-index.
 */
class FMD : public RLCSA
{
  public:
    // We can only be constructed on a previously generated RLCSA index that
    // just happens to meet our requirements.
    explicit FMD(const std::string& base_name, bool print = false);
    
    /**
     * Extend a search by a character, either backward or forward. Ranges are in
     * BWT coordinates.
     */
    FMDPosition extend(FMDPosition range, usint c, bool backward) const;
    
    /**
     * Retract a search by a character, either backward or forward. Reverses a
     * call to extend, but can also retract one way when the extend call was
     * made going the other way. Ranges are in BWT coordinates.
     * 
     * TODO: Doesn't actually work yet.
     */
    FMDPosition retract(FMDPosition range, usint c, bool backward) const;
    
    /**
     * Count occurrences of a pattern using the FMD search algorithm, iterating
     * through the pattern either forward or backward.
     */
    FMDPosition fmdCount(const std::string& pattern, bool backward = true)
      const;
      
    /**
     * Try left-mapping the given index in the given string, starting from
     * scratch. Start a backwards search at that index in the string and extend
     * left until we map to exactly one or zero places. Returns true or false
     * depending on whether we map, an FMDPosition (in BWT coordinates) that, if
     * nonempty, can be extended right to try and map the next base to the
     * right, and the number of characters in the pattern used to make that
     * FMDPosition.
     *
     * If the mapping succeeded, the FMDPosition returned has one thing in it,
     * which is the mapping upstream context.
     *
     * Index must be a valid character position in the string.
     */
    MapAttemptResult mapPosition(const std::string& pattern,
      usint index) const;
      
    /**
     * Attempt to map each base in the query string to a (text, position) pair.
     * The vector returned will have one entry for each character in the
     * selected range.
     * 
     * Optionally a start and length for the region to map can be specified. The
     * whole string will be used as context, but only that region will actually
     * be mapped. A length of -1 means to use the entire string after the start,
     * and is the default.
     */
    std::vector<Mapping> map(const std::string& query, usint start = 0,
      sint length = -1) const;
      
  private:
    /**
     * Get an FMDPosition covering the whole SA.
     */
    FMDPosition getSAPosition() const;
    
    /**
     * Get an FMDPosition for the part of the BWT for things starting with the
     * given character.
     */
    FMDPosition getCharPosition(usint c) const;
    
    /**
     * Convert an FMDPosition in BWT coordinates to one in SA coordinates, in
     * place.
     */
    void convertToSAPosition(FMDPosition& bwt_position) const;
  
    // These are not allowed.
    FMD();
    FMD(const FMD&);
    FMD& operator = (const FMD&);
};

}

#endif
