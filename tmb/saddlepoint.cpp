#include <TMB.hpp>
#include <fenv.h>
#include "utilities.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Survey setup.
  DATA_MATRIX(traps);
  // Activity centres.
  DATA_VECTOR(s);
  // Model parameters.
  DATA_SCALAR(lambda0);
  DATA_SCALAR(sigma);
  // Indicator for nonzero counts.
  DATA_IVECTOR(nonzero);
  // Function to return.
  DATA_STRING(fn);
  // Argument to the KGF.
  PARAMETER_VECTOR(t);
  Type obj = 0;
  if (fn == "kgf_x_given_s"){
    obj += kgf_x_given_s(t, s, traps, lambda0, sigma, nonzero);
  }
  return obj;
}
