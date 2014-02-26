#include "fmd.h"

namespace CSA
{

FMDPosition::FMDPosition(usint forard_start, usint reverse_start, usint length):
  forward_start(forward_start), reverse_start(reverse_start), length(length)
{
}

FMDPosition::FMDPosition(): forward_start(0), reverse_start(0), length(0)
{
}

FMD::FMD(const std::string& base_name, bool print): 
  RLCSA(base_name, print)
{
}

FMDPosition
FMD::extend(FMDPosition range, usint c, bool backward) const
{

  // More or less directly implemented off of algorithms 2 and 3 in "Exploring
  // single-sample SNP and INDEL calling with whole-genome de novo assembly"
  // (Li, 2012). However, our character indices are one less, since we don't
  // allow search patterns to include the end-of-text symbol.
  
  if(backward)
  {

    // Only allow characters in the index
    if(c >= CHARS || this->array[c] == 0) { return EMPTY_FMD_POSITION; }
    // Only allow DNA bases
    if(!isBase(c)) { return EMPTY_FMD_POSITION; }
    
    // We have an array of FMDPositions, one per base, that we will fill in by a
    // tiny dynamic programming.
    FMDPosition answers[NUM_BASES];
    
    // Since we don't keep an FMDPosition for the non-base end-of-text
    // character, we need to track its length separately in order for the DP
    // algorithm given in the paper to be implementable. We calculate
    // occurrences of the text end character analytically, since we know there
    // are this->number_of_sequences of them, at the start of the BWT space.
    usint endOfTextLength = std::min((range.forward_start + range.length - 1), 
      this->number_of_sequences) - std::min(this->number_of_sequences, 
      (range.forward_start - 1));

    for(usint base = 0; base < NUM_BASES; base++)
    {
      
      // Count up the number of characters < this base, including sequence stop
      // characters.
      usint start = this->alphabet->cumulative((usint)BASES[base]) + 
        this->number_of_sequences - 1;
        
      // Get an iterator for the bit vector for this character, for calculating
      // ranks/occurrences.
      PsiVector::Iterator iter(*(this->array[(usint)BASES[base]]));
    
      // Fill in the forward-strand start positions and range lengths for each
      // base's answer.
      // TODO: do we want at_least set or not? What does it do?
      answers[base].forward_start = start + iter.rank(range.forward_start, true);
      answers[base].length = iter.rank(range.forward_start + range.length, 
        true) - iter.rank(range.forward_start, true);
    }
    
    // Set up the dynamic programming for the reverse start, to figure out where
    // the corresponding reverse intervals are.
    answers[3].reverse_start = range.reverse_start + endOfTextLength;
    
    for(usint base = 3; base >= 1; base--)
    {
      answers[base].reverse_start = answers[base + 1].reverse_start +
        answers[base + 1].length;
    }

    // N comes after (before?) everything for some reason. TODO: why? Is it
    // something to do with the encoding Heng uses?
    answers[4].reverse_start = answers[0].reverse_start + answers[0].length;

    // Return the correct answer for the base we actually want
    for(usint base = 0; base < NUM_BASES; base++)
    {
      if(BASES[base] == (char)c)
      {
        // We found the right index for this character. Return that answer
        return answers[base];
      }
    }
    
    // If we ger here, they gave us something not in BASES somehow.
    throw "Unrecognized base";
  }
  else
  {
  
    // Flip the input interval. TODO: why does GCC think these things may have
    // uninitialized members?
    FMDPosition reverse(range.reverse_start, range.forward_start, range.length);
    
    // Do backwards search with the reverse complement of the base
    FMDPosition extended = this->extend(reverse, reverse_complement(c),
      true);
    
    // Reverse the interval again
    FMDPosition forward(extended.reverse_start, extended.reverse_start,
      extended.length);
    return forward;
  
  }
}

FMDPosition
FMD::getSAPosition() const
{
  return FMDPosition(0, 0, this->data_size - 1);
}

FMDPosition
FMD::fmdCount(const std::string& pattern) const
{
  if(pattern.length() == 0) { return this->getSAPosition(); }

  std::string::const_reverse_iterator iter = pattern.rbegin();
  FMDPosition index_position = this->getCharPosition((uchar)*iter);
  if(index_position.length <= 0) { return index_position; }

  for(++iter; iter != pattern.rend(); ++iter)
  {
    // Backwards extend with subsequent characters.
    index_position = this->extend(index_position, (uchar)*iter, true);
    if(index_position.length <= 0) { return EMPTY_FMD_POSITION; }
  }
  this->convertToSAPosition(index_position);

  return index_position;
}

FMDPosition
FMD::getCharPosition(usint c) const
{
  if(c >= CHARS) { return EMPTY_FMD_POSITION; }
  if(!isBase(c)) { return EMPTY_FMD_POSITION; }
  
  pair_type forward_range = this->alphabet->getRange(c);
  this->convertToBWTRange(forward_range);
  
  pair_type reverse_range = this->alphabet->getRange(reverse_complement(c));
  this->convertToBWTRange(reverse_range);
  
  // Make the FMDPosition rolling together both ranges.
  // TODO: Make sure both ranges are the same length, as they should be.
  FMDPosition position(forward_range.first, reverse_range.first,
    forward_range.second - forward_range.first);
  
  return position;
}

void
FMD::convertToSAPosition(FMDPosition& bwt_position) const
{
  bwt_position.forward_start -= this->number_of_sequences;
  bwt_position.reverse_start -= this->number_of_sequences;
}

}
