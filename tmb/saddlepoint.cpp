#include <TMB.hpp>
#include <fenv.h>
#include "utilities.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Reading in data.
  DATA_IVECTOR(counts);
  DATA_MATRIX(mask);
  DATA_MATRIX(traps);
  Type obj = 0;
  return -obj;
}
