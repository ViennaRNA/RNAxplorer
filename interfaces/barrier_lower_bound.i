/**********************************************/
/* BEGIN interface for barrier_lowerbound functions   */
/**********************************************/

%typemap(out) float barrier_estimate_2D %{
	$result = PyFloat_FromDouble($1);
%}

%{
extern "C" {
#include "../src/barrier_lower_bound.h"
}
%}

%include "../src/barrier_lower_bound.h"
