/**********************************************/
/* BEGIN interface for RNAxplorer functions   */
/**********************************************/
%module RNAxplorer

%feature("autodoc", "3");

%{
extern "C" {
#include "../src/barrier_lower_bound.h"
#include "../src/distorted_sampling.h"
#include "../src/distorted_samplingMD.h"
#include "../src/dist_class_sc.h"
#include "../src/meshpoint.h"
#include "../src/paths.h"
#include "../src/RNAwalk.h"
}
%}

%include "typemaps.i"
%include "std_vector.i";
%include "std_string.i";

%template(StringVector) std::vector<std::string>;
%template(DoubleVector) std::vector<double>;

%include barrier_lower_bound.i
%include distorted_sampling.i
%include distorted_samplingMD.i
%include meshpoint.i
%include RNAwalk.i
%include paths.i
