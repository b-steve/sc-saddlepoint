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
  return exp(kgf_x_given_s(t, s, traps, lambda0, sigma, nonzero));
}

// MGF for the marginal distribution of a single individual's counts.
template<class Type>
Type mgf_x(const vector<Type> &t, const matrix<Type> &m, const matrix<Type> &traps, const Type &lambda0, const Type &sigma, const vector<int> nonzero){
  int n_mask = m.rows();
  Type out = 0;
  vector<Type> s(2);
  for (int i = 0; i < n_mask; i++){
    s = m.row(i);
    out += mgf_x_given_s(t, s, traps, lambda0, sigma, nonzero);
  }
  return out/n_mask;
}

// KGF for a single individual's counts.
template<class Type>
Type kgf_x(const vector<Type> &t, const matrix<Type> &m, const matrix<Type> &traps, const Type &lambda0, const Type &sigma, const vector<int> nonzero){
  return log(mgf_x(t, m, traps, lambda0, sigma, nonzero));
}

// KGF for the combined counts.
template<class Type>
Type kgf_c(const vector<Type> &t, const matrix<Type> &m, const matrix<Type> &traps, const Type &lambda0, const Type &sigma, const Type &EN, const vector<int> nonzero){
  return kgf_pois(kgf_x(t, m, traps, lambda0, sigma, nonzero), EN);
}

// MGF for the combined counts.
template<class Type>
Type mgf_c(const vector<Type> &t, const matrix<Type> &m, const matrix<Type> &traps, const Type &lambda0, const Type &sigma, const Type &EN, const vector<int> nonzero){
  return exp(kgf_c(t, m, traps, lambda0, sigma, EN, nonzero));
}
  
#endif
