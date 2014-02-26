#include <iostream>

#include "fmd.h"
#include "misc/utils.h"

using namespace CSA;


enum mode_type { COUNT, TOTAL, START, RELATIVE, DISPLAY, CONTEXT };


void printUsage()
{
  std::cout << "Usage: fmd_grep [-c|-t|-s|-r|-NUM] pattern base_name" << std::endl;
  std::cout << "  -c    print the number of matching sequences" << std::endl;
  std::cout << "  -t    print the total number of occurrences" << std::endl;
  std::cout << "  -s    print the start positions of matches" << std::endl;
  std::cout << "  -r    print the relative start positions of matches (sequence, position)" << std::endl;
  std::cout << "  -NUM  display NUM characters of leading and trailing context instead of" << std::endl;
  std::cout << "        the entire line" << std::endl;
}


int main(int argc, char** argv)
{
  int base_arg = 2, pattern_arg = 1;
  mode_type mode = DISPLAY;
  usint context = 0;

  if(argc < base_arg + 1)
  {
    printUsage();
    return 1;
  }

  if(argv[1][0] == '-')
  {
    base_arg++; pattern_arg++;
    if(std::string("-c").compare(argv[1]) == 0)
    {
      mode = COUNT;
    }
    else if(std::string("-t").compare(argv[1]) == 0)
    {
      mode = TOTAL;
    }
    else if(std::string("-s").compare(argv[1]) == 0)
    {
      mode = START;
    }
    else if(std::string("-r").compare(argv[1]) == 0)
    {
      mode = RELATIVE;
    }
    else
    {
      mode = CONTEXT;
      context = atoi(&(argv[1][1]));
    }
    if(argc < base_arg + 1)
    {
      printUsage();
      return 2;
    }
  }

  FMD fmd(argv[base_arg], false);
  if(!fmd.isOk())
  {
    return 3;
  }

  usint len = std::string(argv[pattern_arg]).length();
  pair_type result_range = fmd.count(argv[pattern_arg]);
  usint occurrences = length(result_range);
  
  // Try the alternative double-ended search
  FMDPosition fmd_result = fmd.fmdCount(argv[pattern_arg]);
  
  std::cout << "Got " << fmd_result.length << " FMD matches, " << occurrences << " RLCSA matches" << std::endl;
  
  if(mode == TOTAL)
  {
    std::cout << occurrences << std::endl;
    return 0;
  }
  
  if(occurrences == 0)
  {
    if(mode == COUNT)
    {
      std::cout << 0 << std::endl;
    }
    return 0;
  }

  usint last_row = 0;
  usint* results = fmd.locate(result_range);
  if(mode == COUNT || mode == DISPLAY)
  {
    // Make results hold text numbers instead of positions.
    fmd.getSequenceForPosition(results, occurrences);
  }
  std::sort(results, results + occurrences);
  if(mode == COUNT || mode == DISPLAY)
  {
    // Find each text only once, even if it matches multiple times.
    for(usint i = 1; i < occurrences; i++)
    {
      if(results[i] != results[last_row])
      {
        last_row++; results[last_row] = results[i];
      }
    }
  }

  if(mode == COUNT)
  {
    std::cout << (last_row + 1) << std::endl;
  }
  else if(mode == DISPLAY)
  {
    for(usint i = 0; i <= last_row; i++)
    {
      uchar* row = fmd.display(results[i]);
      std::cout.write((char*)row, length(fmd.getSequenceRange(results[i])));
      std::cout << std::endl;
      delete[] row;
    }
  }
  else if(mode == START)
  {
    for(usint i = 0; i < occurrences; i++)
    {
      std::cout << results[i] << std::endl;
    }
  }
  else if(mode == RELATIVE)
  {
    for(usint i = 0; i < occurrences; i++)
    {
      pair_type relative = fmd.getRelativePosition(results[i]);
      std::cout << relative.first << ", " << relative.second << std::endl;
    }
  }
  else if(mode == CONTEXT)
  {
    usint result_length = 0;
    for(usint i = 0; i < occurrences; i++)
    {
      uchar* text = fmd.display(results[i], len, context, result_length);
      std::cout.write((char*)text, result_length);
      std::cout << std::endl;
      delete[] text;
    }
  }

  delete[] results;
  return 0;
}
