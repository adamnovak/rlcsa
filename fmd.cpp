#include "fmd.h"
#include <algorithm>
#include <stdexcept>

namespace CSA
{

const std::string BASES = "TGCNA";
const std::string ALPHABETICAL_BASES = "ACGNT";

FMDPosition::FMDPosition(usint forward_start, usint reverse_start,
  usint end_offset): forward_start(forward_start), reverse_start(reverse_start),
  end_offset(end_offset)
{
}

FMDPosition::FMDPosition(): forward_start(0), reverse_start(0), end_offset(-1)
{
}

FMDPosition
FMDPosition::flip() const
{
  // Swap the two intervals of the bi-interval
  return FMDPosition(reverse_start, forward_start, end_offset);
}

bool FMDPosition::operator==(const FMDPosition& other) const
{
  // Compare all the fields.
  return
    forward_start == other.forward_start &&
    reverse_start == other.reverse_start &&
    end_offset == other.end_offset;
}

sint
FMDPosition::range(const RangeVector& ranges) const
{
  // Get an iterator for making rank queries.
  // TODO: Is it efficient to do this a lot?
  RangeVector::Iterator iter(ranges);

  // Look up the range that the forward starting position is in
  sint start_range = iter.rank(forward_start);
  
  // And the range the forward end is in
  sint end_range = iter.rank(forward_start + end_offset);
  
  if(start_range == end_range)
  {
    // Both ends of the interval are in the same range.
    return start_range;
  }
  else
  {
    // There is no range spanning our forward-strand interval.
    return -1;
  }
}

sint
FMDPosition::ranges(const RangeVector& ranges) const
{
  // Get an iterator for making rank queries.
  // TODO: Is it efficient to do this a lot?
  RangeVector::Iterator iter(ranges);

  // Look up the range that the starting position is in
  sint start_range = iter.rank(forward_start);
  
  // And the range the end is in
  sint end_range = iter.rank(forward_start + end_offset);
  
  // Return the number of ranges we intersect (1s hit plus 1)
  return(end_range - start_range + 1);
}

std::ostream& operator<< (std::ostream& o, FMDPosition const& position)
{
  // Report both the ranges that we represent.
  return o << position.forward_start << "-" << 
    (position.forward_start + position.end_offset) << "|" << 
    position.reverse_start << "-" << (position.reverse_start + 
    position.end_offset);
}

Mapping::Mapping(): location(0, 0), is_mapped(false)
{
}

Mapping::Mapping(pair_type location, bool is_mapped): location(location), 
  is_mapped(is_mapped)
{
}

bool
Mapping::operator==(const Mapping& other) const
{
  return location == other.location && is_mapped == other.is_mapped;
}

std::ostream&
operator<< (std::ostream& o, Mapping const& mapping)
{
  if(mapping.is_mapped)
  {
    o << "Text " << mapping.location.first << " offset " <<
      mapping.location.second;
  }
  else
  {
    o << "-----------------";
  }
  return o;
}

// Stuff for FMDIterators that traverse the suffix tree.

FMDIterator::FMDIterator(const FMD& parent, usint depth, bool beEnd,
    bool reportDeadEnds): 
  parent(parent), depth(depth), reportDeadEnds(reportDeadEnds), stack(),
  pattern()
{
  // By default we start out with empty everything, which is what we should have
  // at the end.
  
  if(!beEnd) 
  {
    // Start at the beginning.
    search();
  }
}

FMDIterator::FMDIterator(const FMDIterator& toCopy): parent(toCopy.parent), 
  depth(toCopy.depth), reportDeadEnds(toCopy.reportDeadEnds),
  stack(toCopy.stack), pattern(toCopy.pattern)
{
  // Already made a duplicate stack. Nothing to do.
}

FMDIterator& FMDIterator::operator++()
{
  // Advance
  search();
  
  // Return ourselves.
  return *this;
}

FMDIterator FMDIterator::operator++(int)
{
  // Copy ourselves
  FMDIterator copy(*this);
  
  // Advance
  search();
  
  // Return the copy
  return copy;
}

std::pair<std::string, FMDPosition> FMDIterator::operator*() const 
{ 
  // Just grab the value we stored.
  return toYield;
}

void FMDIterator::yield(std::pair<std::string, FMDPosition> value) {
  // Save the value for returning when dereferenced.
  toYield = value;
}

bool FMDIterator::operator==(const FMDIterator& other) const
{
  return 
    // We have the same parent addresses
    &parent == &(other.parent) && 
    // And go to the same depth
    depth == other.depth && 
    // And both report dead ends or not
    reportDeadEnds == other.reportDeadEnds &&
    // And are at the same depth
    stack.size() == other.stack.size() && 
    // And followed the same path to get there
    std::equal(stack.begin(), stack.end(), other.stack.begin()) && 
    // And we have the same string pattern (which the above should imply)
    pattern == other.pattern;
}

bool FMDIterator::operator!=(const FMDIterator& other) const
{
  // Just use the equality check.
  return !(*this == other);
}

void FMDIterator::search()
{
  // TODO: Unify with tryRecurseToDepth by making its topDepth a parameter.

  if(stack.size() == 0) 
  {
    // Try to recurse down to depth, starting with base 0 at the root branch.
    tryRecurseToDepth(0);
    
    // Now we are either at the correct depth at the leftmost suffix tree node,
    // or in the state we started in (which is equal to end). So we are done.
  } 
  else if(stack.size() == depth) 
  {
    // We were at the right depth, meaning we don't need to recurse further.
    
    do
    {
      // Pop the bottom frame (already explored).
      std::pair<FMDPosition, usint> lastFrame = pop();
      
      // Try recursing to the next thing at that same level and then down to the
      // required depth.
      if(tryRecurseToDepth(lastFrame.second + 1))
      {
        // We got something at the required depth, or ran into a shorter end of
        // text suffix to report.
        return;
      }
      
      // Otherwise, we didn't get something to yield in this subtree. That means
      // this subtree has been exhausted and we should pop and try the next one
      // over, which we will do on the next loop iteration.
    }
    while(stack.size() > 0);
    
    // If we get here, we have gone all the way back up and popped and tried
    // replacing everything with no more results. That means we have finished
    // the search, and should be equal to end.
    
  } else if(reportDeadEnds) {
    // We weren't at the right depth; we were at a node that happened to show up
    // followed by a text end somewhere. Continue recursing from here.
    
    if(tryRecurseToDepth(0)) {
      // We can just go down from here.
      return;
    }
    
    // Otherwise we need to pop up and continue the main recursion loop from
    // above. TODO: Unify.
    
    do
    {
      // Pop the bottom frame (already explored).
      std::pair<FMDPosition, usint> lastFrame = pop();
      
      // Try recursing to the next thing at that same level and then down to the
      // required depth.
      if(tryRecurseToDepth(lastFrame.second + 1))
      {
        // We got something at the required depth, or another shorter suffix
        // followed by end of text.
        return;
      }
      
      // Otherwise, we didn't get something at the required depth in this
      // subtree. That means this subtree has been exhausted and we should pop
      // and try the next one over, which we will do on the next loop iteration.
    }
    while(stack.size() > 0);    
    
  } else {
    // We broke something.
    throw std::runtime_error("Iterator was at wrong depth");
  }
}

bool FMDIterator::recurse(usint baseNumber)
{
  if(baseNumber >= NUM_BASES)
  {
    // Don't even try out-of-range bases.
    return false;
  }
  
  // What will we extend to?
  FMDPosition extension;
  
  if(stack.size() == 0)
  {
    // Our "extension" is just starting with this base.
    extension = parent.getCharPosition(ALPHABETICAL_BASES[baseNumber]);
  }
  else
  {
    // Work out what we would select if we extended forwards with this letter
    // (i.e. appended it to the suffix).
    extension = parent.extend(stack.back().first, 
      ALPHABETICAL_BASES[baseNumber], false);
  }
  
  
  
  if(extension.isEmpty())
  {
    // This would be a suffix that doesn't appear.
    return false;
  }
  
  // This would not be an empty place. Go there.
  // Add a stack frame
  stack.push_back(std::make_pair(extension, baseNumber));
  // And record the change to the pattern.
  pattern.push_back(ALPHABETICAL_BASES[baseNumber]);
  
  return true;
}

bool FMDIterator::tryRecurse(usint baseNumber)
{
  // Try every base number after the given starting one until we either recurse
  // or run out.
  for(; baseNumber < NUM_BASES && !recurse(baseNumber); baseNumber++);
  
  // If we didn't run out, we must have successfully recursed.
  return baseNumber < NUM_BASES;
}

bool FMDIterator::tryRecurseToDepth(usint baseNumber)
{
  // What depth should we not go above?
  size_t topDepth = stack.size();
  
  while(stack.size() < depth)
  {
    // Until we get to the right depth...
    
    // We never go into an empty thing.
    // So we must be at a shallower depth. Try to go deeper.
    
    if(tryRecurse(baseNumber)) {
      // We made it down. Reset baseNumber
      baseNumber = 0;
      
      if(reportDeadEnds && stack.size() < depth) {
        // If we lose some range here (i.e. some positions are followed by end
        // of text and don't show up in an extension with any base), we have to
        // stop so that we yield them. Those positions will be the first ones in
        // our range if they exist.
        
        // See what we would get if we extended. TODO: Can we show that this
        // will always work? It seems to work when there are none of the first
        // base, and should always work when there is some of the first base,
        // but I'm not sure it won't break.
        FMDPosition extension = parent.extend(stack.back().first, 
          ALPHABETICAL_BASES[0], false);
          
        if(extension.forward_start != stack.back().first.forward_start) {
          // There are some suffixes that come into this node and leave before
          // the first real base. They must end the text.
          
          INFO(
            std::cout << "End of text: " << pattern << "$" << std::endl;
            std::cout << extension << " vs. " << stack.back().first <<
            std::endl;
          )
          
          // We got to a place we want to yield.
          
          // Grab the FMDPosition we want to return
          FMDPosition toConvert = stack.back().first; 
          
          // But before converting it to SA coordinates, fix it up to indicate
          // only the part we aren't covering. Forward start is going to be the
          // same, but it is going to run until the start of extension. Subtract
          // 1 to keep this as an offset where 0 = a 1-base itnerval.
          toConvert.end_offset = extension.forward_start - 
            toConvert.forward_start - 1;
            
          // TODO: The reverse start can't be moved sanely and is thus going to
          // be invalid (we can't search anchored to text start).
          
          // We need to convert it from BWT coordinates (which we use internally
          // for extension) to SA coordinates (which are natural for locate).
          parent.convertToSAIndex(toConvert.forward_start);
          parent.convertToSAIndex(toConvert.reverse_start);
          
          // Grab the pattern and the FMDPosition from the top of the stack to
          // go with it. Yield the pair.
          yield(std::make_pair(pattern, toConvert));
          
          return true;
        }
        
      }
      
    } else {
      // We can't go to anything nonempty down. So we should try going up.
      
      if(stack.size() == topDepth) {
        // We don't want to go any higher than this. Give up on finding a
        // nonempty thing at the right depth in this subtree.
        return false;
      }
      
      // Otherwise, go up and continue on the next branch.
      baseNumber = pop().second + 1;
    }
  }
  
  // We got to the right depth, and we never go to empty things.
  
  // Grab the FMDPosition we want to return
  FMDPosition toConvert = stack.back().first; 
  
  // We need to convert it from BWT coordinates (which we use internally for
  // extension) to SA coordinates (which are natural for locate).
  parent.convertToSAIndex(toConvert.forward_start);
  parent.convertToSAIndex(toConvert.reverse_start);
  
  // Grab the pattern and the FMDPosition from the top of the stack to go with
  // it. Yiueld the pair.
  yield(std::make_pair(pattern, toConvert));
  
  return true;
}

std::pair<FMDPosition, usint> FMDIterator::pop()
{
  // Drop the character that this added to the pattern.
  pattern.erase(pattern.size() - 1, 1);
  
  // Grab the stack frame to return
  std::pair<FMDPosition, usint> toReturn(stack.back());
  
  // Drop it from the stack.
  stack.pop_back();
  
  // Return it
  return toReturn;
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
    if(c >= CHARS || this->array[c] == 0) {
        DEBUG(std::cout << "Character " << c << " not in index." << std::endl;)
        return EMPTY_FMD_POSITION;
    }
    // Only allow DNA bases
    if(!isBase(c)) {
        DEBUG(std::cout << "Character " << c << " is not a DNA base." << 
            std::endl;)
        return EMPTY_FMD_POSITION;
    }
    
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
        // end_offset will get used and probably needs to be -1 for empty.
        answers[base].end_offset = -1;
        
      }
      else
      {
        DEBUG(std::cout << "\t\tCharacter appeared." << std::endl;)
        
        // Get an iterator for the bit vector for this character, for
        // calculating ranks/occurrences.
        PsiVector::Iterator iter(*vector);
        
        DEBUG(std::cout << "\t\tGot iterator" << std::endl;)
      
        // Fill in the forward-strand start positions and range end_offsets for
        // each base's answer. TODO: do we want at_least set or not? What does
        // it do?
        
        // First cache the forward_start rank we re-use
        usint forward_start_rank = iter.rank(range.forward_start, true);
        
        answers[base].forward_start = start + forward_start_rank;
        answers[base].end_offset = iter.rank(range.forward_start + 
          range.end_offset, false) - forward_start_rank;
          
      }
        
      
        
      DEBUG(std::cout << "\t\tWould go to: " << answers[base].forward_start <<
        "-" << (sint)answers[base].forward_start + answers[base].end_offset << 
        " length " << answers[base].getLength() << std::endl;)
    }
    
    // Since we don't keep an FMDPosition for the non-base end-of-text
    // character, we need to track its length separately in order for the DP
    // algorithm given in the paper to be implementable. We calculate
    // occurrences of the text end character (i.e. how much of the current range
    // is devoted to things where an endOfText comes next) implicitly: it's
    // whatever part of the length of the range is unaccounted-for by the other
    // characters. We need to use the length accessor because ranges with one
    // thing have the .end_offset set to 0.
    usint endOfTextLength = range.getLength();
    
    for(usint base = 0; base < NUM_BASES; base++)
    {
      // Go through the bases in order and account for their lengths.
      endOfTextLength -= answers[base].getLength();
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
        answers[base - 1].getLength();
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
    // Flip the interval, do backwards search with the reverse complement of the
    // base, and then flip back.
    return this->extend(range.flip(), reverse_complement(c), true).flip();
  
  }
}

FMDPosition
FMD::retract(FMDPosition range, usint c, bool backward) const
{

  // TODO: Does not work at all.
  
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
    
    DEBUG(std::cout << "\tOriginal forward range contains " <<
      original.forward_start << std::endl;)
      
    // Work out something from reverse range
      
    return original;
    
  }
  else
  {
  
    // Flip the interval, do backwards retract with the reverse complement of
    // the base, and then flip back.
    return this->retract(range.flip(), reverse_complement(c), true).flip();
  
  }
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
    if(index_position.isEmpty()) { return index_position; }

    DEBUG(std::cout << "Starting with " << index_position << std::endl;)

    for(++iter; iter != pattern.rend(); ++iter)
    {
      // Backwards extend with subsequent characters.
      index_position = this->extend(index_position, *iter, true);
      DEBUG(std::cout << "Now at " << index_position << " after " << *iter <<
        std::endl;)
      // Test out retracting
      //DEBUG(this->retract(index_position, *iter, true);)
      if(index_position.isEmpty()) { return EMPTY_FMD_POSITION; }
    }
  }
  else 
  {
    // Start at the front of the pattern and work towards the end.
    
    std::string::const_iterator iter = pattern.begin();
    index_position = this->getCharPosition((uchar)*iter);
    if(index_position.isEmpty()) { return index_position; }

    DEBUG(std::cout << "Starting with " << index_position << std::endl;)

    for(++iter; iter != pattern.end(); ++iter)
    {
      // Forwards extend with subsequent characters.
      index_position = this->extend(index_position, *iter, false);
      DEBUG(std::cout << "Now at " << index_position << " after " << *iter << 
        std::endl;)
      // Test out retracting
      //DEBUG(this->retract(index_position, *iter, false);)
      if(index_position.isEmpty()) { return EMPTY_FMD_POSITION; }
    }
    
  }
  
  this->convertToSAPosition(index_position);

  return index_position;
}

std::pair<pair_type, usint>
FMD::countUntilUnique(const std::string& pattern, usint index) const
{
  // Mostly copied from the RLCSA count method.
  
  if(pattern.length() == 0 || index == 0) {
    return std::make_pair(this->getSARange(), 0);
  }

  // Start at the index we want to map.
  sint i = index;
  // And with the range from just that character.
  pair_type index_range = this->getCharRange((uchar)pattern[i]);
  
  // Make sure we aren't empty already.
  if(isEmpty(index_range)) {
    return std::make_pair(index_range, index - i + 1);
  }
  
  if(index_range.first == index_range.second) { 
    // We found a unique place. Return it.
    this->convertToSARange(index_range);
    return std::make_pair(index_range, index - i + 1);
  }

  for(--i; i >= 0; --i)
  {
    // For each base going left...
  
    // Apply the LF mapping to shrink the range.
    index_range = this->LF(index_range, (uchar)pattern[i]);
    
    // Stop if we're empty.
    if(isEmpty(index_range)) {
      return std::make_pair(EMPTY_PAIR, index - i + 1);
    }
    
    if(index_range.first == index_range.second) { 
      // We found a unique place. Return it.
      this->convertToSARange(index_range);
      return std::make_pair(index_range, index - i + 1);
    }
    
  }
  
  // If we get here, we hit the start and still aren't unique. Report our too-
  // big range.
  this->convertToSARange(index_range);
  return std::make_pair(index_range, index - i);
}


MapAttemptResult
FMD::mapPosition(const std::string& pattern, usint index) const
{
  DEBUG(std::cout << "Mapping " << index << " in " << pattern << std::endl;)
  
  // Initialize the struct we will use to return our somewhat complex result.
  // Contains the FMDPosition (which we work in), an is_mapped flag, and a
  // variable counting the number of extensions made to the FMDPosition.
  MapAttemptResult result;
  
  // Do a backward search.
  // Start at the given index, and get the starting range for that character.
  result.is_mapped = false;
  result.position = this->getCharPosition(pattern[index]);
  result.characters = 1;
  if(result.position.isEmpty())
  {
    // This character isn't even in it. Just return the result with an empty
    // FMDPosition; the next character we want to map is going to have to deal
    // with having some never-before-seen character right upstream of it.
    return result;
  }
  else if(result.position.getLength() == 1)
  {
    // We've already mapped.
    result.is_mapped = true;
    return result;
  }

  if(index == 0) {
    // The rest of the function deals with characters to the left of the one we
    // start at. If we start at position 0 there can be none.
    return result;
  }
  
  DEBUG(std::cout << "Starting with " << result.position << std::endl;)

  do
  {
    // Now consider the next character to the left.
    index--;
    
    // Grab the character to extend with.
    usint character = pattern[index];
    
    DEBUG(std::cout << "Index " << index << " in " << pattern << " is " << 
        (char) character << "(" << character << ")" << std::endl;)
    
    // Backwards extend with subsequent characters.
    FMDPosition next_position = this->extend(result.position, character,
      true);
    extends++;
      
    DEBUG(std::cout << "Now at " << next_position << " after " << 
      pattern[index] << std::endl;)
    if(next_position.isEmpty())
    {
      // The next place we would go is empty, so return the result holding the
      // last position.
      return result;
    }
    else if(next_position.getLength() == 1)
    {
      // We have successfully mapped to exactly one place. Update our result to
      // reflect the additional extension and our success, and return it.
      result.position = next_position;
      result.characters++;
      result.is_mapped = true;
      return result;      
    }
    
    // Otherwise, we still map to a plurality of places. Record the extension
    // and loop again.
    result.position = next_position;
    result.characters++;
    
  }
  while(index > 0);
  // Continue considering characters to the left until we hit the start of the
  // string.
  
  // If we get here, we ran out of upstream context and still map to multiple
  // places. Just give our multi-mapping FMDPosition and unmapped result.
  return result;

}

MapAttemptResult
FMD::mapPosition(const RangeVector& ranges, const std::string& pattern, 
  usint index) const
{
  // We're going to right-map so ranges match up with the things we can map to
  // (downstream contexts)

  // Initialize the struct we will use to return our somewhat complex result.
  // Contains the FMDPosition (which we work in), an is_mapped flag, and a
  // variable counting the number of extensions made to the FMDPosition.
  MapAttemptResult result;
  
  // Do a forward search.
  // Start at the given index, and get the starting range for that character.
  result.is_mapped = false;
  result.position = this->getCharPosition(pattern[index]);
  result.characters = 1;
  if(result.position.isEmpty())
  {
    // This character isn't even in it. Just return the result with an empty
    // FMDPosition; the next character we want to map is going to have to deal
    // with having some never-before-seen character right upstream of it.
    return result;
  }
  else if (result.position.range(ranges) != -1)
  {
    // We've already mapped.
    result.is_mapped = true;
    return result;
  }

  DEBUG(std::cout << "Starting with " << result.position << std::endl;)

  for(index++; index < pattern.size(); index++)
  {
    // Forwards extend with subsequent characters.
    FMDPosition next_position = this->extend(result.position, pattern[index],
      false);
    extends++;
      
    DEBUG(std::cout << "Now at " << next_position << " after " << 
      pattern[index] << std::endl;)
    if(next_position.isEmpty())
    {
      // The next place we would go is empty, so return the result holding the
      // last position.
      return result;
    }
    
    if(next_position.range(ranges) != -1)
    {
      // We have successfully mapped to exactly one range. Update our result to
      // reflect the additional extension and our success, and return it.
      result.position = next_position;
      result.characters++;
      result.is_mapped = true;
      return result;      
    }
    
    // Otherwise, we still map to a plurality of ranges. Record the extension
    // and loop again.
    result.position = next_position;
    result.characters++;
  }
  
  // If we get here, we ran out of downstream context and still map to multiple
  // ranges. Just give our multi-mapping FMDPosition and unmapped result.
  return result;

}

std::vector<Mapping>
FMD::map(const std::string& query, usint start, sint length) const
{

  if(length == -1) {
    // Fix up the length parameter if it is -1: that means the whole rest of the
    // string.
    length = query.length() - start;
  }
  
  // We need a vector to return.
  std::vector<Mapping> mappings;
  
  // Keep around the result that we get from the single-character mapping
  // function. We use it as our working state to track our FMDPosition and how
  // many characters we've extended by. We use the is_mapped flag to indicate
  // whether the current iteration is an extension or a restart.
  MapAttemptResult location;
  // Make sure the scratch position is empty so we re-start on the first base.
  // Other fields get overwritten.
  location.position = EMPTY_FMD_POSITION;
  
  for(sint i = start; i < (sint)(start + length); i++)
  {
    if(location.position.isEmpty())
    {
      INFO(std::cout << "Starting over by mapping position " << i <<
        std::endl;)
      // We do not currently have a non-empty FMDPosition to extend. Start over
      // by mapping this character by itself.
      location = this->mapPosition(query, i);
      restarts++;
    }
    else
    {
      INFO(std::cout << "Extending with position " << i << std::endl;)
      // The last base either mapped successfully or failed due to multi-
      // mapping. Try to extend the FMDPosition we have to the right (not
      // backwards) with the next base.
      location.position = this->extend(location.position, query[i], false);
      extends++;
      location.characters++;
    }
    
    if(location.is_mapped && location.position.getLength() == 1)
    {
      // It mapped. We didn't do a re-start and fail, and there's exactly one
      // thing in our interval.
        
      // Take the first (only) thing in the bi-interval's forward strand side,
      // and convert to SA coordinates.
      usint converted_start = location.position.forward_start;
      convertToSAIndex(converted_start);
      
      // Locate it, and then report position as a (text, offset) pair. This will
      // give us the position of the first base in the pattern, which lets us
      // infer the position of the last base in the pattern.
      pair_type text_location = getRelativePosition(locate(converted_start));
        
      INFO(std::cout << "Mapped " << location.characters << 
        " context to text " << text_location.first << " position " << 
        text_location.second << std::endl;)
        
      // Correct to the position of the last base in the pattern, by offsetting
      // by the length of the pattern that was used. A 2-character pattern means
      // we need to go 1 further right in the string it maps to to find where
      // its rightmost character maps.
      text_location.second += (location.characters - 1);
      
      // Add a Mapping for this mapped base.
      mappings.push_back(Mapping(text_location));
      
      // We definitely have a non-empty FMDPosition to continue from
      
    }
    else
    {
    
      INFO(std::cout << "Failed (" << location.position.getLength() << 
        " options for " << location.characters << " context)." << std::endl;)
        
      if(location.is_mapped && location.position.isEmpty())
      {
        // We extended right until we got no results. We need to try this base
        // again, in case we tried with a too-long left context.
        
        INFO(std::cout << "Restarting from here..." << std::endl;)
        
        // Move the loop index back
        i--;
        
        // Since the FMDPosition is empty, on the next iteration we will retry
        // this base.
        
      }
      else
      {
        // It didn't map for some other reason:
        // - It was an initial mapping with too little left context to be unique
        // - It was an initial mapping with a nonexistent left context
        // - It was an extension that was multimapped and still is
        
        // In none of these cases will re-starting from this base help at all.
        // If we just restarted here, we don't want to do it again. If it was
        // multimapped before, it had as much left context as it could take
        // without running out of string or getting no results.
      
        // It didn't map. Add an empty/unmapped Mapping.
        mappings.push_back(Mapping());
        
        // Mark that the next iteration will be an extension (if we had any
        // results this iteration; if not it will just restart)
        location.is_mapped = true;
        
      }
    }
  
  }
  
  // We've gone through and attempted the whole string. Give back our answers.
  return mappings;
  
}

std::vector<Mapping>
FMD::mapFM(const std::string& query, usint start, sint length) const
{

  // The same as the above function, but without using extend.

  if(length == -1) {
    // Fix up the length parameter if it is -1: that means the whole rest of the
    // string.
    length = query.length() - start;
  }
  
  // We need a vector to return.
  std::vector<Mapping> mappings;
  
  for(sint i = start; i < (sint)(start + length); i++)
  {
    // For each base to map...

    // Count left from there until we don't need to any more.
    std::pair<pair_type, usint> countResult = countUntilUnique(query, i);
    pair_type range = countResult.first;
    usint characters = countResult.second;
    
    if(range.first == range.second) {
      // We successfully mapped to just one place.
      
      // Locate it, and then report position as a (text, offset) pair. This will
      // give us the position of the first base in the pattern, which lets us
      // infer the position of the last base in the pattern.
      pair_type text_location = getRelativePosition(locate(range.first));
        
      INFO(std::cout << "Mapped to text " << text_location.first << 
        " position " << text_location.second << std::endl;)
        
      // Correct to the position of the last base in the pattern, by offsetting
      // by the length of the pattern that was used. A 2-character pattern means
      // we need to go 1 further right in the string it maps to to find where
      // its rightmost character maps.
      text_location.second += (characters - 1);
      
      // Add a Mapping for this mapped base.
      mappings.push_back(Mapping(text_location));
      
    } else {
      // We mapped to multiple or no places: no results
      mappings.push_back(Mapping());
      
    }
  }
  
  // We've gone through and attempted the whole string. Give back our answers.
  return mappings;
  
}


std::vector<sint> 
FMD::map(const RangeVector& ranges, const std::string& query, usint start, 
  sint length) const {
  
  // RIGHT-map to a range.
    
  if(length == -1) {
    // Fix up the length parameter if it is -1: that means the whole rest of the
    // string.
    length = query.length() - start;
  }
  
  // We need a vector to return.
  std::vector<sint> mappings;
  
  // Keep around the result that we get from the single-character mapping
  // function. We use it as our working state to trackour FMDPosition and how
  // many characters we've extended by. We use the is_mapped flag to indicate
  // whether the current iteration is an extension or a restart.
  MapAttemptResult location;
  // Make sure the scratch position is empty so we re-start on the first base
  location.position = EMPTY_FMD_POSITION;
  
  for(sint i = start + length - 1; i >= (sint) start; i--)
  {
    // Go from the end of our selected region to the beginning.
    
    if(location.position.isEmpty())
    {
      INFO(std::cout << "Starting over by mapping position " << i <<
        std::endl;)
      // We do not currently have a non-empty FMDPosition to extend. Start over
      // by mapping this character by itself.
      location = this->mapPosition(ranges, query, i);
      restarts++;
    }
    else
    {
      INFO(std::cout << "Extending with position " << i << std::endl;)
      // The last base either mapped successfully or failed due to multi-
      // mapping. Try to extend the FMDPosition we have to the left (backwards)
      // with the next base.
      location.position = this->extend(location.position, query[i], true);
      extends++;
      location.characters++;
    }
    
    // What range index does our current left-side position (the one we just
    // moved) correspond to, if any?
    sint range = location.position.range(ranges);
    
    if(location.is_mapped && !location.position.isEmpty() && range != -1)
    {
      // It mapped. We didn't do a re-start and fail, and our interval is
      // nonempty and subsumed by a range.
      
      INFO(std::cout << "Mapped " << location.characters << 
        " context to range #" << range << " in range vector." << std::endl;)
      
      // Remember that this base mapped to this range
      mappings.push_back(range);
      
      // We definitely have a non-empty FMDPosition to continue from
      
    }
    else
    {
    
      INFO(std::cout << "Failed (" << location.position.ranges(ranges) << 
        " options for " << location.characters << " context)." << std::endl;)
        
      if(location.is_mapped && location.position.isEmpty())
      {
        // We extended right until we got no results. We need to try this base
        // again, in case we tried with a too-long left context.
        
        INFO(std::cout << "Restarting from here..." << std::endl;)
        
        // Move the loop index towards the end we started from (right)
        i++;
        
        // Since the FMDPosition is empty, on the next iteration we will retry
        // this base.
        
      }
      else
      {
        // It didn't map for some other reason:
        // - It was an initial mapping with too little left context to be unique
        //   to a range.
        // - It was an initial mapping with a nonexistent left context
        // - It was an extension that was multimapped and still is
        
        // In none of these cases will re-starting from this base help at all.
        // If we just restarted here, we don't want to do it again. If it was
        // multimapped before, it had as much left context as it could take
        // without running out of string or getting no results.
      
        // It didn't map. Say it corresponds to no range.
        mappings.push_back(-1);
        
        // Mark that the next iteration will be an extension (if we had any
        // results this iteration; if not it will just restart)
        location.is_mapped = true;
        
      }
    }
  
  }
  
  // We've gone through and attempted the whole string. Put our results in the
  // same order as the string, instead of the backwards order we got them in.
  // See <http://www.cplusplus.com/reference/algorithm/reverse/>
  std::reverse(mappings.begin(), mappings.end());
  
  // Give back our answers.
  return mappings;
    
}

FMD::iterator FMD::begin(usint depth, bool reportDeadEnds) const
{
  // Make a new suffix tree iterator that automatically searches out the first
  // suffix of the right length.
  return FMD::iterator(*this, depth, false, reportDeadEnds);
}
     
FMD::iterator FMD::end(usint depth, bool reportDeadEnds) const
{
  // Make a new suffix tree iterator that is just a 1-past-the-end sentinel.
  return FMD::iterator(*this, depth, true, reportDeadEnds);
}

usint FMD::extends = 0;
usint FMD::restarts = 0;

pair_type
FMD::getStats() {
  // Pair up the stats.
  pair_type toReturn = std::make_pair(extends, restarts);
  
  // Clear the static counters
  extends = restarts = 0;
  
  return toReturn;
}

FMDPosition
FMD::getSAPosition() const
{
  return FMDPosition(0, 0, this->data_size - 1);
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
