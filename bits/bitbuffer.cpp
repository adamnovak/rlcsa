#include <cstdlib>
#include <cstring>

#include "bitbuffer.h"


namespace CSA
{

//--------------------------------------------------------------------------

ReadBuffer::ReadBuffer(std::ifstream& file, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(true)
{
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  file.read((char*)buffer, this->size * sizeof(usint));
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(FILE* file, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(true)
{
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  if(file != 0)
  {
    if(!std::fread(buffer, this->size * sizeof(usint), 1, file)) { return; }
  }
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(std::ifstream& file, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(true)
{
  this->size = bitsToWords(this->items * this->item_bits);
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  file.read((char*)buffer, this->size * sizeof(usint));
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(FILE* file, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(true)
{
  this->size = bitsToWords(this->items * this->item_bits);
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  if(file != 0)
  {
    if(!std::fread(buffer, this->size * sizeof(usint), 1, file)) { return; }
  }
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(const usint* buffer, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(false)
{
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(const usint* buffer, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(false)
{
  this->size = bitsToWords(this->items * this->item_bits);
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(const ReadBuffer& original) :
  data(original.data),
  size(original.size),
  item_bits(original.item_bits),
  items(original.items),
  free_buffer(false)
{
  this->reset();
}

ReadBuffer::~ReadBuffer()
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
}

//--------------------------------------------------------------------------

void
ReadBuffer::claimData()
{
  this->free_buffer = true;
}

void
ReadBuffer::writeTo(std::ofstream& file) const
{
  file.write((const char*)this->data, this->size * sizeof(usint));
}

void
ReadBuffer::writeTo(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(this->data, this->size * sizeof(usint), 1, file);
}

void
ReadBuffer::moveBuffer(const usint* buffer)
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
  this->free_buffer = false;

  this->data = buffer;
  this->reset();
}

usint
ReadBuffer::reportSize() const
{
  usint bytes = sizeof(*this);
  if(this->free_buffer) { bytes += this->size * sizeof(usint); }
  return bytes;
}

//--------------------------------------------------------------------------

WriteBuffer::WriteBuffer(usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(true)
{
  this->data = new usint[words];
  memset(this->data, 0, this->size * sizeof(usint));
  this->reset();
}

WriteBuffer::WriteBuffer(usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(true)
{
  this->size = bitsToWords(this->items * this->item_bits);
  this->data = new usint[this->size];
  memset(this->data, 0, this->size * sizeof(usint));
  this->reset();
}

WriteBuffer::WriteBuffer(usint* buffer, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(false)
{
  this->data = buffer;
  this->reset();
}

WriteBuffer::WriteBuffer(usint* buffer, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(false)
{
  this->size = bitsToWords(this->items * this->item_bits);
  this->data = buffer;
  this->reset();
}

WriteBuffer::~WriteBuffer()
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
}

//--------------------------------------------------------------------------

ReadBuffer*
WriteBuffer::getReadBuffer()
{
  ReadBuffer* buffer;
  if(this->items > 0)
  {
    buffer = new ReadBuffer(this->data, this->items, this->item_bits);
  }
  else
  {
    buffer = new ReadBuffer(this->data, this->size);
  }

  if(this->free_buffer)
  {
    buffer->claimData();
    this->free_buffer = false;
  }

  return buffer;
}

void
WriteBuffer::writeTo(std::ofstream& file) const
{
  file.write((char*)this->data, this->size * sizeof(usint));
}

void
WriteBuffer::writeTo(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(this->data, this->size * sizeof(usint), 1, file);
}

void
WriteBuffer::moveBuffer(usint* buffer)
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
  this->free_buffer = false;

  this->data = buffer;
  this->reset();
}

usint
WriteBuffer::reportSize() const
{
  usint bytes = sizeof(*this);
  if(this->free_buffer) { bytes += this->size * sizeof(usint); }
  return bytes;
}

//--------------------------------------------------------------------------

} // namespace CSA
