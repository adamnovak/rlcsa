CXX = g++

# Use 64-bit integers in a 64-bit environment.
SIZE_FLAGS = -DMASSIVE_DATA_RLCSA

# Parallelism is supported by either libstdc++ Parallel Mode or MCSTL.
PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -D_GLIBCXX_PARALLEL -fopenmp
# MCSTL_ROOT = /fs-3/d/jltsiren/suds/mcstl
# PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -I$(MCSTL_ROOT)/c++ -fopenmp

# Vectors using nibble codes instead of delta codes are faster, but they also
# take up more space.
VECTOR_FLAGS = $(PSI_FLAGS) $(LCP_FLAGS) $(SA_FLAGS)
# PSI_FLAGS = -DUSE_NIBBLE_VECTORS
# LCP_FLAGS = -DSUCCINCT_LCP_VECTOR
# SA_FLAGS = -DSUCCINCT_SA_VECTOR
DEBUG_FLAGS = -g -fno-inline

# Flags to use for SWIG. Adjust for your platform
SWIG_FLAGS = -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -fno-strict-aliasing
JAVA_PACKAGE = fi.helsinki.cs.rlcsa


CXXFLAGS = -Wall -O3 -fPIC $(DEBUG_FLAGS) $(SIZE_FLAGS) $(PARALLEL_FLAGS) $(VECTOR_FLAGS)
OBJS = rlcsa.o rlcsa_builder.o fmd.o sasamples.o alphabet.o \
lcpsamples.o sampler.o suffixarray.o adaptive_samples.o docarray.o \
bits/array.o bits/bitbuffer.o bits/multiarray.o bits/bitvector.o bits/deltavector.o \
bits/rlevector.o bits/nibblevector.o bits/succinctvector.o misc/parameters.o misc/utils.o
SWIG_OBJS = rlcsa_wrap.o fmd_wrap.o

PROGRAMS = rlcsa_test lcp_test parallel_build build_rlcsa merge_rlcsa build_sa \
locate_test display_test document_graph read_bwt extract_sequence rlcsa_grep fmd_grep \
build_plcp sample_lcp sampler_test ss_test utils/extract_text utils/convert_patterns \
utils/split_text utils/sort_wikipedia utils/genpatterns

VPATH = bits:misc:utils


default: parallel_build build_rlcsa merge_rlcsa rlcsa_test sampler_test display_test \
document_graph

jar: rlcsa.jar

rlcsa.jar: rlcsa.so RLCSANativeLoader.java
	mkdir -p jar
	javac java/*.java RLCSANativeLoader.java -d jar
	# Make the directory for the Java package
	mkdir -p jar/`echo "$(JAVA_PACKAGE)" | sed s/\\\\./\\\\//g`
	cp rlcsa.so jar/`echo "$(JAVA_PACKAGE)" | sed s/\\\\./\\\\//g`/
	jar cf $@ -C jar .

# Install the jar in the Maven local repository. See
# http://maven.apache.org/guides/mini/guide-3rd-party-jars-local.html
jar-install: rlcsa.jar
	mvn install:install-file -Dfile=rlcsa.jar -DgroupId=fi.helsinki.cs \
	-DartifactId=rlcsa -Dversion=1.0.0-SNAPSHOT -Dpackaging=jar

librlcsa.a: $(OBJS)
	ar rcs librlcsa.a $(OBJS)
	
rlcsa.so: $(OBJS) $(SWIG_OBJS)
	$(CXX) $(LDFLAGS) -shared -o rlcsa.so  $(OBJS) $(SWIG_OBJS)

depend:
	g++ -MM *.cpp bits/*.cpp misc/*.cpp utils/*.cpp > dependencies.mk

rlcsa_test: rlcsa_test.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o rlcsa_test rlcsa_test.o librlcsa.a

lcp_test: lcp_test.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o lcp_test lcp_test.o librlcsa.a

parallel_build: parallel_build.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o parallel_build parallel_build.o librlcsa.a

build_rlcsa: build_rlcsa.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o build_rlcsa build_rlcsa.o librlcsa.a
	
merge_rlcsa: merge_rlcsa.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o merge_rlcsa merge_rlcsa.o librlcsa.a

build_sa: build_sa.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o build_sa build_sa.o librlcsa.a

locate_test: locate_test.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o locate_test locate_test.o librlcsa.a

display_test: display_test.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o display_test display_test.o librlcsa.a

document_graph: document_graph.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o $@ $@.o librlcsa.a

read_bwt: read_bwt.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o read_bwt read_bwt.o librlcsa.a

extract_sequence: extract_sequence.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o extract_sequence extract_sequence.o librlcsa.a

rlcsa_grep: rlcsa_grep.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o rlcsa_grep rlcsa_grep.o librlcsa.a
	
fmd_grep: fmd_grep.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o fmd_grep fmd_grep.o librlcsa.a

build_plcp: build_plcp.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o build_plcp build_plcp.o librlcsa.a

sample_lcp: sample_lcp.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o sample_lcp sample_lcp.o librlcsa.a

sampler_test: sampler_test.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o sampler_test sampler_test.o librlcsa.a

ss_test: ss_test.o librlcsa.a
	$(CXX) $(CXXFLAGS) -o ss_test ss_test.o librlcsa.a

extract_text: extract_text.o
	$(CXX) $(CXXFLAGS) -o utils/extract_text extract_text.o

convert_patterns: convert_patterns.o utils.o
	$(CXX) $(CXXFLAGS) -o utils/convert_patterns convert_patterns.o misc/utils.o

split_text: split_text.o utils.o
	$(CXX) $(CXXFLAGS) -o utils/split_text split_text.o misc/utils.o

sort_wikipedia: sort_wikipedia.o utils.o
	$(CXX) $(CXXFLAGS) -o utils/sort_wikipedia sort_wikipedia.o misc/utils.o

genpatterns: genpatterns.c
	gcc -O3 -Wall -o utils/genpatterns utils/genpatterns.c

# SWIG C++ file generation and compilation.	
%_wrap.o: %_wrap.cxx
	$(CXX) $(INCLUDES) $(SWIG_FLAGS) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $^
	
%_wrap.cxx: %.i
	mkdir -p java
	swig -c++ -java -outdir java -package $(JAVA_PACKAGE) $(SIZE_FLAGS) $(VECTOR_FLAGS) $^

clean:
	rm -f librlcsa.a
	rm -f rlcsa.so
	rm -f $(PROGRAMS)
	rm -f *.o bits/*.o misc/*.o utils/*.o
	rm -rf java/
	rm -f *_wrap.cxx
	rm -rf jar/
	rm -f rlcsa.jar

package:
	mkdir rlcsa
	mkdir rlcsa/bits rlcsa/misc rlcsa/utils
	cp LICENSE Makefile README dependencies.mk *.cpp *.h *.i *.java rlcsa
	cp bits/*.cpp bits/*.h rlcsa/bits
	cp misc/*.cpp misc/*.h rlcsa/misc
	cp utils/*.cpp utils/*.py rlcsa/utils
	tar cfz rlcsa.tgz rlcsa
	rm -r rlcsa/*
	rmdir rlcsa

deploy:
	find . -type d -exec chmod 770 {} \;
	find . -type f -not \( -perm +111 \) -exec chmod 660 {} \;
	find . -type f -perm +111 -exec chmod 770 {} \;


include dependencies.mk
