#ifndef _RLCSA_SUCCINCTVECTOR_H
#define _RLCSA_SUCCINCTVECTOR_H

#include <fstream>

#include "bitvector.h"


namespace CSA
{


/*
  This class is used to construct a SuccinctVector.
  Unlike other bit vectors, we can set any bit in the current superblock, even if we have
  already set bits after it.
*/

class SuccinctEncoder : public VectorEncoder
{
  public:
    SuccinctEncoder(usint block_bytes, usint superblock_size = VectorEncoder::SUPERBLOCK_SIZE);
    ~SuccinctEncoder();

    void setBit(usint value);
    void setRun(usint start, usint len);

    void addBit(usint value);
    void addRun(usint start, usint len);
    void flush();

  protected:

    // These are not allowed.
    SuccinctEncoder();
    SuccinctEncoder(const SuccinctEncoder&);
    SuccinctEncoder& operator = (const SuccinctEncoder&);
};


/*
  This is a succinct bit vector.
*/

class SuccinctVector : public BitVector
{
  public:
    const static usint SHORT_RANGE = 16;  // Should be at least 2.

    typedef SuccinctEncoder Encoder;

    explicit SuccinctVector(std::ifstream& file);
    explicit SuccinctVector(FILE* file);
    SuccinctVector(Encoder& encoder, usint universe_size);
    explicit SuccinctVector(Encoder& encoder); // Use the array directly.
    explicit SuccinctVector(WriteBuffer& vector);
    ~SuccinctVector();

//--------------------------------------------------------------------------

    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;
    usint reportSize() const;

    // The default implementation removes rank_index that is required
    // for select() in SuccinctVector.
    void strip() {}

//--------------------------------------------------------------------------

    class Iterator
    {
      public:
        explicit Iterator(const SuccinctVector& par);
        ~Iterator();

        usint rank(usint value, bool at_least = false);

        usint select(usint index);
        usint selectNext();
        inline bool hasNext() const { return (this->cur < this->parent.items - 1); }

        pair_type valueBefore(usint value);
        pair_type valueAfter(usint value);
        pair_type nextValue();

        pair_type selectRun(usint index, usint max_length);
        pair_type selectNextRun(usint max_length);

        bool isSet(usint value);

        usint countRuns();  // Not implemented.

      protected:
        const SuccinctVector& parent;
        usint cur;

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------

  protected:

    // How many 1-bits are in the previous blocks.
    void indexForRank();

    // Which block contains the 1-bit of rank i * select_rate.
    void indexForSelect();

    // These are not allowed.
    SuccinctVector();
    SuccinctVector(const SuccinctVector&);
    SuccinctVector& operator = (const SuccinctVector&);
};


} // namespace CSA


#endif // _RLCSA_SUCCINCTVECTOR_H
