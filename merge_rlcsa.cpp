#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif

#include "rlcsa_builder.h"
#include "misc/utils.h"


using namespace CSA;


double getRLCSA(RLCSABuilder& builder, const std::string& base_name);

const int MAX_THREADS = 64;


int
main(int argc, char** argv)
{
  std::cout << "RLCSA merger" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: merge_rlcsa original additional [threads]" << std::endl;
    return 1;
  }

  int original_parameter = 1, additional_parameter = 2, threads_parameter = 3;

  std::string additional_name = argv[additional_parameter];
  std::cout << "Index to add: " << additional_name << std::endl;

  std::string base_name = argv[original_parameter];
  std::cout << "Index to update: " << base_name << std::endl;

  usint threads = 1;
  if(argc > threads_parameter)
  {
    threads = std::min(MAX_THREADS, std::max(atoi(argv[threads_parameter]), 1));
  }
  std::cout << "Threads: " << threads << std::endl; 
  std::cout << std::endl;

  std::string parameters_name = base_name + PARAMETERS_EXTENSION;
  Parameters parameters;
  parameters.set(RLCSA_BLOCK_SIZE);
  parameters.set(SAMPLE_RATE);
  parameters.set(SUPPORT_LOCATE);
  parameters.set(SUPPORT_DISPLAY);
  parameters.set(WEIGHTED_SAMPLES);
  parameters.read(parameters_name);
  parameters.print();

  double start = readTimer();
  double megabytes = 0.0;
  
  
  std::cout << "Merging the indexes" << std::endl;
  
  double mark = readTimer();
  std::cout << "Load: " << base_name; std::cout.flush();
  RLCSA* originalIndex = new RLCSA(base_name);
  RLCSABuilder builder(parameters.get(RLCSA_BLOCK_SIZE), parameters.get(SAMPLE_RATE), 0, threads, originalIndex);
  std::cout << " (" << (readTimer() - mark) << " seconds)" << std::endl;
  
  mark = readTimer();
  std::cout << "Increment: " << additional_name; std::cout.flush();
  builder.insertFromFile(additional_name);
  std::cout << " (" << (readTimer() - mark) << " seconds)" << std::endl;
  
  std::cout << std::endl;
  megabytes = getRLCSA(builder, base_name);

  double stop = readTimer();
  std::cout << megabytes << " megabytes indexed in " << (stop - start) << " seconds (" << (megabytes / (stop - start)) << " MB/s)." << std::endl;
  std::cout << "Search time:   " << builder.getSearchTime() << " seconds" << std::endl;
  std::cout << "Sort time:     " << builder.getSortTime() << " seconds" << std::endl;
  std::cout << "Merge time:    " << builder.getMergeTime() << " seconds" << std::endl;
  std::cout << "Memory usage:  " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return 0;
}


double
getRLCSA(RLCSABuilder& builder, const std::string& base_name)
{
  double megabytes = 0.0;

  RLCSA* index = builder.getRLCSA();
  if(index != 0 && index->isOk())
  {
    index->printInfo();
    index->reportSize(true);
    index->writeTo(base_name);
    megabytes = index->getSize() / (double)MEGABYTE;
  }

  delete index;
  return megabytes;
}
