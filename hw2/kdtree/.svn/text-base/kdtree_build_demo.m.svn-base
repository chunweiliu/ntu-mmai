% KDTREE_BUILD_DEMO: illustrate the functionalities of kdtree_build
clc, clear, close all;

mex kdtree_build.cpp
disp('compiled.');

p = rand( 10, 3 );
tree = kdtree_build( p );
disp(sprintf('tree structure pointer: %d',tree ));
