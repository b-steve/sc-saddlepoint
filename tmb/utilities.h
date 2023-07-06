#ifndef UTILITIES_H
#define UTILITIES_H

// Stable Poisson PMF.
template<class Type>
Type dpois_stable (const int &x, const Type &lambda, const int &give_log){
  double d_x = x;
  Type out;
  out = pow(lambda, d_x)*exp(-lambda)/exp(lgamma(d_x + 1));
  if (give_log){
    out = log(out + DBL_MIN);
  }
  return out;
}

// KGF for the Poisson distribution.
template<class Type>
Type kgf_pois(const Type &t, const Type &lambda){
  return lambda*(exp(t) - 1);
}

// MGF for the Poisson distribution.
template<class Type>
Type mgf_pois(const Type &t, const Type &lambda){
  return exp(kgf_pois(t, lambda));
}

// KGF for a single individual's counts.
template<class Type>
Type kgf_x_given_s(const vector<Type> &t, const vector<Type> &s, const matrix<Type> &traps, const Type &lambda0, const Type &sigma, const vector<int> nonzero){
  int n_traps = traps.rows();
  Type d;
  Type lambda;
  Type out = 0;
  int k = 0;
  for (int i = 0; i < n_traps; i++){
    d = pow(pow(s(0) - traps(i, 0), 2) +
	    pow(s(1) - traps(i, 1), 2), 0.5);
    lambda = lambda0*exp(-pow(d, 2)/(2*pow(sigma, 2)));
    if (nonzero(i) == 1){
      out += kgf_pois(t(k), lambda);
      k++;
    } else {
      out += dpois_stable(0, lambda, 1);
    }
  }
  return out;
}

// MGF for a single individual's counts.
template<class Type>
Type mgf_x_given_s(const vector<Type> &t, const vector<Type> &s, const matrix<Type> &traps, const Type &lambda0, const Type &sigma, const vector<int> nonzero){
  return exp(mgf_x_given_s(t, s, traps, lambda0, sigma, nonzero));
}

#endif
