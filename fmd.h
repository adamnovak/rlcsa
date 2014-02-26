#ifndef FMD_H
#define FMD_H

#include <fstream>
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

const static std::string BASES = "ACGTN";

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
 */
struct FMDPosition {
  usint forward_start;
  usint reverse_start;
  usint length;
  FMDPosition();
  FMDPosition(usint forard_start, usint reverse_start, usint length);
};

const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, 0);

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
     * Count occurrences of a pattern using the FMD search algorithm.
     */
    FMDPosition fmdCount(const std::string& pattern) const;
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
