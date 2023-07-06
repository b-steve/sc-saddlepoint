#include <TMB.hpp>
#include <fenv.h>
#include "utilities.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Counts.
  DATA_VECTOR(c);
  // Survey setup.
  DATA_MATRIX(traps);
  // Mask points.
  DATA_MATRIX(m);
  // Model parameters.
  DATA_SCALAR(EN);
  DATA_SCALAR(lambda0);
  DATA_SCALAR(sigma);
  // Indicator for nonzero counts.
  DATA_IVECTOR(nonzero);
  // Indicator for KGF.
  DATA_INTEGER(is_kgf);
  // Argument to the KGF.
  PARAMETER_VECTOR(t);
  // Number of traps.
  int n_traps = traps.rows();
  // Vector of nonzero counts.
  int n_nonzero = sum(nonzero);
  vector<Type> c_nonzero(n_nonzero);
  int k = 0;
  for (int i = 0; i < n_traps; i++){
    if (nonzero(i) == 1){
      c_nonzero(k) = c(i);
      k++;
    }
  }
  Type out;
  out = kgf_c(t, m, traps, lambda0, sigma, EN, nonzero);
  if (is_kgf == 0){
    out -= (t*c_nonzero).sum();
  }
  return out;
}
