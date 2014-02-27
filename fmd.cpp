#include "fmd.h"

namespace CSA
{

FMDPosition::FMDPosition(usint forward_start, usint reverse_start, usint length):
  forward_start(forward_start), reverse_start(reverse_start), length(length)
{
}

FMDPosition::FMDPosition(): forward_start(0), reverse_start(0), length(-1)
{
}

std::ostream& operator<< (std::ostream& o, FMDPosition const& position)
{
  // Report both the ranges that we represent.
  return o << position.forward_start << "-" << 
    (position.forward_start + position.length) << "|" << 
    position.reverse_start << "-" << (position.reverse_start + position.length);
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
  // allow search patterns to include the end-of-text symbol. We also use
  // alphabetical ordering instead of the paper's N-last ordering in the FM-
  // index, and consequently need to assign reverse ranges in alphabetical order
  // by reverse complement.
  
  if(backward)
  {

    // Only allow characters in the index
    if(c >= CHARS || this->array[c] == 0) { return EMPTY_FMD_POSITION; }
    // Only allow DNA bases
    if(!isBase(c)) { return EMPTY_FMD_POSITION; }
    
    DEBUG(std::cout << "Extending " << range << " backwards with " << (char)c <<
      std::endl;)
    
    // We have an array of FMDPositions, one per base, that we will fill in by a
    // tiny dynamic programming.
    FMDPosition answers[NUM_BASES];
    
    for(usint base = 0; base < NUM_BASES; base++)
    {
      // Go through the bases in arbitrary order.
      
      DEBUG(std::cout << "\tThinking about base " << base << "(" << 
        BASES[base] << ")" << std::endl;)
      
      // Count up the number of characters < this base, including sequence stop
      // characters.
      usint start = this->alphabet->cumulative((usint)BASES[base]) + 
        this->number_of_sequences - 1;
        
      DEBUG(std::cout << "\t\tstart = " << start << std::endl;)
      
      // Get a pointer to the bit vector for this letter, which might be NULL if
      // this base never appeared.
      PsiVector* vector = this->array[(usint)BASES[base]];
      
      if(vector == NULL)
      {
        DEBUG(std::cout << "\t\tCharacter never appeared!" << std::endl;)
        
        // Fill in forward_start and length with the knowledge that this
        // character doesn't exist. forward_start should never get used, but
        // length will get used and probably needs to be -1 for empty.
        answers[base].length = -1;
        
      }
      else
      {
        DEBUG(std::cout << "\t\tCharacter appeared." << std::endl;)
        
        // Get an iterator for the bit vector for this character, for
        // calculating ranks/occurrences.
        PsiVector::Iterator iter(*vector);
        
        DEBUG(std::cout << "\t\tGot iterator" << std::endl;)
      
        // Fill in the forward-strand start positions and range lengths for each
        // base's answer. TODO: do we want at_least set or not? What does it do?
        answers[base].forward_start = start + iter.rank(range.forward_start,
          true);
        answers[base].length = iter.rank(range.forward_start + range.length, 
          false) - iter.rank(range.forward_start, true);
          
        // Make sure rank and select work reasonably.
          
        for(int i = -2; i < range.length + 1; i++) {
          DEBUG(std::cout << "\t\trank(" << range.forward_start + i << ", true)=" <<
            iter.rank(range.forward_start + i, true) << std::endl;)
        }
          
        usint rank = iter.rank(range.forward_start, true);
        
        for(int i = rank; i >= 0; i--) {
        
          DEBUG(std::cout << "\t\tselect(" << i << ")=" << 
            iter.select(i) << std::endl;)
        }
        
      }
        
      
        
      DEBUG(std::cout << "\t\tWould go to: " << answers[base].forward_start <<
        "-" << (sint)answers[base].forward_start + answers[base].length << 
        " length " << length(answers[base]) << std::endl;)
    }
    
    // Since we don't keep an FMDPosition for the non-base end-of-text
    // character, we need to track its length separately in order for the DP
    // algorithm given in the paper to be implementable. We calculate
    // occurrences of the text end character (i.e. how much of the current range
    // is devoted to things where an endOfText comes next) implicitly: it's
    // whatever part of the length of the range is unaccounted-for by the other
    // characters. We need to use the length accessor because ranges with one
    // thing have the .length set to 0.
    usint endOfTextLength = length(range);
    
    for(usint base = 0; base < NUM_BASES; base++)
    {
      // Go through the bases in arbitrary order and account for their lengths.
      endOfTextLength -= length(answers[base]);
    }
    
      
    DEBUG(std::cout << "\tendOfTextLength = " << endOfTextLength << std::endl;)
    
    // The endOfText character is the very first character we need to account
    // for when subdividing the reverse range and picking which subdivision to
    // take.
    DEBUG(std::cout << "\tendOfText reverse_start would be " << 
      range.reverse_start << std::endl;)
    
    // Next, allocate the range for the base that comes first in alphabetical
    // order by reverse complement.
    answers[0].reverse_start = range.reverse_start + endOfTextLength;
    DEBUG(std::cout << "\t" << BASES[0] << " reverse_start is " << 
      answers[0].reverse_start << std::endl;)
    
    for(usint base = 1; base < NUM_BASES; base++)
    {
      // For each subsequent base in alphabetical order by reverse complement
      // (as stored in BASES), allocate it the next part of the reverse range.
      
      answers[base].reverse_start = answers[base - 1].reverse_start + 
        length(answers[base - 1]);
      DEBUG(std::cout << "\t" << BASES[base] << " reverse_start is " << 
        answers[base].reverse_start << std::endl;)
      
    }
    
    // Now all the per-base answers are filled in.
    
    for(usint base = 0; base < NUM_BASES; base++)
    {
      // For each base in arbitrary order
      if(BASES[base] == (char)c)
      {
        DEBUG(std::cout << "Moving " << range << " to " << answers[base] << 
          " on " << BASES[base] << std::endl;)
        // This is the base we're actually supposed to be extending with. Return
        // its answer.
        return answers[base];
      }
    }
    
    // If we get here, they gave us something not in BASES somehow.
    throw "Unrecognized base";
  }
  else
  {
  
    // Flip the input interval.
    FMDPosition reverse(range.reverse_start, range.forward_start, range.length);
    
    // Do backwards search with the reverse complement of the base
    FMDPosition extended = this->extend(reverse, reverse_complement(c),
      true);
    
    // Reverse the interval again
    FMDPosition forward(extended.reverse_start, extended.forward_start,
      extended.length);
    return forward;
  
  }
}

FMDPosition
FMD::retract(FMDPosition range, usint c, bool backward) const
{

  
  if(backward)
  {
    DEBUG(std::cout << "Going back from " << range << " on " << (char)c <<
      std::endl;)
  
    // Keep the original FMDPosition to build up. We call it "original" because
    // we logically think about undoing an extension, but in reality it may
    // never have existed.
    FMDPosition original;
  
    // Get a pointer to the bit vector for this letter, which might be NULL if
    // this base never appeared.
    PsiVector* vector = this->array[c];
    
    if(vector == NULL) { throw "Character never appeared!"; }
      
    // Get an iterator for the bit vector for this character, for selecting by
    // rank.
    PsiVector::Iterator iter(*vector);
  
    // Count up the number of characters < this base, including sequence stop
    // characters. The same as the "start" variable in extend.
    usint start = this->alphabet->cumulative(c) + this->number_of_sequences - 1;
    
    DEBUG(std::cout << "\tOriginal start was " << start << std::endl;)
      
    // Back-derive the original forward_start, which can be done with "start"
    // and the inverse of rank(i, true).
    original.forward_start = iter.select(range.forward_start - start - 1);
    
    DEBUG(std::cout << "\tOriginal forward_start was " <<
      original.forward_start << std::endl;)
      
    return original;
    
  }
  else
  {
  
    // Flip the input interval.
    FMDPosition reverse(range.reverse_start, range.forward_start, range.length);
    
    // Do backwards retraction with the reverse complement of the base
    FMDPosition retracted = this->retract(reverse, reverse_complement(c),
      true);
    
    // Reverse the interval again
    FMDPosition forward(retracted.reverse_start, retracted.forward_start,
      retracted.length);
    return forward;
  
  }
}

FMDPosition
FMD::getSAPosition() const
{
  return FMDPosition(0, 0, this->data_size - 1);
}

FMDPosition
FMD::fmdCount(const std::string& pattern, bool backward) const
{
  DEBUG(std::cout << "Counting " << pattern << std::endl;)

  if(pattern.length() == 0) { return this->getSAPosition(); }

  // Keep an FMDPosition to store our intermediate result in.
  FMDPosition index_position;

  if(backward)
  {
    // Start at the end of the pattern and work towards the front
    
    std::string::const_reverse_iterator iter = pattern.rbegin();
    index_position = this->getCharPosition((uchar)*iter);
    if(isEmpty(index_position)) { return index_position; }

    DEBUG(std::cout << "Starting with " << index_position << std::endl;)

    for(++iter; iter != pattern.rend(); ++iter)
    {
      // Backwards extend with subsequent characters.
      index_position = this->extend(index_position, *iter, true);
      DEBUG(std::cout << "Now at " << index_position << " after " << *iter <<
        std::endl;)
      // Test out retracting
      this->retract(index_position, *iter, true);
      if(isEmpty(index_position)) { return EMPTY_FMD_POSITION; }
    }
  }
  else 
  {
    // Start at the front of the pattern and work towards the end.
    
    std::string::const_iterator iter = pattern.begin();
    index_position = this->getCharPosition((uchar)*iter);
    if(isEmpty(index_position)) { return index_position; }

    DEBUG(std::cout << "Starting with " << index_position << std::endl;)

    for(++iter; iter != pattern.end(); ++iter)
    {
      // Forwards extend with subsequent characters.
      index_position = this->extend(index_position, *iter, false);
      DEBUG(std::cout << "Now at " << index_position << " after " << *iter << 
        std::endl;)
      // Test out retracting
      this->retract(index_position, *iter, false);
      if(isEmpty(index_position)) { return EMPTY_FMD_POSITION; }
    }
    
  }
  
  this->convertToSAPosition(index_position);

  return index_position;
}

FMDPosition
FMD::getCharPosition(usint c) const
{
  if(c >= CHARS || this->array[c] == 0) { return EMPTY_FMD_POSITION; }
  if(!isBase(c)) { return EMPTY_FMD_POSITION; }
  
  pair_type forward_range = this->alphabet->getRange(c);
  this->convertToBWTRange(forward_range);
  DEBUG(std::cout << (char)c << " range: " << forward_range.first << "-" << 
    forward_range.second << std::endl;)
  
  pair_type reverse_range = this->alphabet->getRange(reverse_complement(c));
  this->convertToBWTRange(reverse_range);
  DEBUG(std::cout << (char)reverse_complement(c) << " range: " << 
    reverse_range.first << "-" << reverse_range.second << std::endl;)
  
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
