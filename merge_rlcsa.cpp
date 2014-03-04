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
    std::cout << "Usage: merge_rlcsa [-threads] original additional [additional2...]" << std::endl;
    return 1;
  }

  int original_parameter = 2, additional_parameter = 3, threads_parameter = 1;

  usint threads = 1;
  if(argc > threads_parameter && argv[threads_parameter][0] == '-')
  {
    threads = std::min(MAX_THREADS, std::max(atoi(argv[threads_parameter] + 1), 1));
  }
  else 
  {
    original_parameter--;
    additional_parameter--;
  }

  if(argc < additional_parameter)
  {
    std::cerr << "Error: specify the index to merge into and one index to add." << std::endl;
    return 1;
  }

  std::string base_name = argv[original_parameter];
  std::cout << "Index to update: " << base_name << std::endl;
  
  std::vector<std::string> additional_names;
  
  for(int i = additional_parameter; i < argc; i++)
  {
    std::string additional_name = argv[i];
    std::cout << "Index to add: " << additional_name << std::endl;
    additional_names.push_back(additional_name);
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
  
  for(std::vector<std::string>::iterator i = additional_names.begin(); i != additional_names.end(); ++i)
  {
    mark = readTimer();
    std::cout << "Increment: " << *i; std::cout.flush();
    builder.insertFromFile(*i);
    std::cout << " (" << (readTimer() - mark) << " seconds)" << std::endl;
  }
  
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
