#CPPONLY flag removes mex-portions (include and matlab iterface function) to allow C++ independent testing
CXXFLAGS = -j2 -O2 -g -Wall -fmessage-length=0 -D CPPONLY

## NOTE!! make these point to the one in the matlab directories
# linux default preferences
ifeq ($(shell uname),Linux)
	MXX = /usr/local/matlab/bin/mex
	MEXEXT = $(shell /usr/local/matlab/bin/mexext) 
	MXXFLAGS = -I/usr/local/matlab/extern/include
endif
# default darwin preferences
ifeq ($(shell uname),Darwin)
    MXX = /Applications/matlab/bin/mex
    MEXEXT = $(shell /Applications/matlab/bin/mexext)
    MXXFLAGS = -I/Applications/matlab/extern/include
endif
# my cardiff@GrUVi lab machine
ifeq ($(shell uname -n),cardiff)
    MXX = /Applications/MATLAB_R2008a/bin/mex
    MEXEXT = $(shell /Applications/MATLAB_R2008a/bin/mexext)
    MXXFLAGS = -I/Applications/MATLAB_R2008a/extern/include
endif

TARGET = kdtree_build kdtree_nearest_neighbor kdtree_range_query kdtree_ball_query kdtree_k_nearest_neighbors
MEXTARGET = kdtree_build.$(MEXEXT) kdtree_nearest_neighbor.$(MEXEXT) kdtree_range_query.$(MEXEXT) kdtree_ball_query.$(MEXEXT) kdtree_k_nearest_neighbors.$(MEXEXT)
# DEVELOPMENT REDUCED TARGETS
# TARGET = 
# MEXTARGET = kdtree_k_nearest_neighbors.$(MEXEXT)

HDRS = KDTree.h MyHeaps.h

all: $(TARGET) $(MEXTARGET)

# -----------------------------------------------------#
kdtree_build: $(HDRS) kdtree_build.cpp
	$(CXX) $(CXXFLAGS) -o kdtree_build kdtree_build.cpp
kdtree_build.$(MEXEXT): $(HDRS) kdtree_build.cpp
	$(MXX) $(MXXFLAGS) -o kdtree_build.$(MEXEXT) kdtree_build.cpp

# -----------------------------------------------------#
kdtree_nearest_neighbor: $(HDRS) kdtree_nearest_neighbor.cpp
	$(CXX) $(CXXFLAGS) -o kdtree_nearest_neighbor kdtree_nearest_neighbor.cpp
kdtree_nearest_neighbor.$(MEXEXT): $(HDRS) kdtree_nearest_neighbor.cpp
	$(MXX) $(MXXFLAGS) -o kdtree_nearest_neighbor.$(MEXEXT) kdtree_nearest_neighbor.cpp

# -----------------------------------------------------#
kdtree_range_query: $(HDRS) kdtree_range_query.cpp
	$(CXX) $(CXXFLAGS) -o kdtree_range_query kdtree_range_query.cpp
kdtree_range_query.$(MEXEXT): $(HDRS) kdtree_range_query.cpp
	$(MXX) $(MXXFLAGS) -o kdtree_range_query.$(MEXEXT) kdtree_range_query.cpp
	
# -----------------------------------------------------#
kdtree_ball_query: $(HDRS) kdtree_ball_query.cpp
	$(CXX) $(CXXFLAGS) -o kdtree_ball_query kdtree_ball_query.cpp
kdtree_ball_query.$(MEXEXT): $(HDRS) kdtree_ball_query.cpp
	$(MXX) $(MXXFLAGS) -o kdtree_ball_query.$(MEXEXT) kdtree_ball_query.cpp

# -----------------------------------------------------#
kdtree_k_nearest_neighbors: $(HDRS) kdtree_k_nearest_neighbors.cpp
	$(CXX) $(CXXFLAGS) -o kdtree_k_nearest_neighbors kdtree_k_nearest_neighbors.cpp
kdtree_k_nearest_neighbors.$(MEXEXT): $(HDRS) kdtree_k_nearest_neighbors.cpp
	$(MXX) $(MXXFLAGS) -o kdtree_k_nearest_neighbors.$(MEXEXT) kdtree_k_nearest_neighbors.cpp

clean:
	rm -f $(OBJS) $(TARGET) $(MEXTARGET)
