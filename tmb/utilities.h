#ifndef UTILITIES_H
#define UTILITIES_H

// Stable Poisson PMF.
template<class Type>
Type dpois_stable (const Type &x, const Type &lambda, const int &give_log){
  Type out;
  out = pow(lambda, x)*exp(-lambda)/exp(lgamma(x + 1));
  if (give_log){
    out = log(out + DBL_MIN);
  }
  return out;
}
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
Type kgf_pois(const Type &t, const double &lambda){
  lambda*(exp(t) - 1);
}

// MGF for the Poisson distribution.
template<class Type>
Type mgf_pois(const Type &t, const double &lambda){
  exp(kgf_pois(t, lambda));
}

// KGF for a single individual's counts.
template<class Type>
Type kgf_x_given_s(const vector<Type> &t, const matrix<double> &s, const matrix<Type> &traps, const double &lambda0, const double &sigma, const vector<double> &nz){
  int n_s = s.size();
  int n_traps = traps.size();
  double d;
  double lambda;
  Type out = 0;
  for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_traps; j++){
      d = pow(pow(s(i, 1) - traps(j, 1), 2) +
	      pow(s(i, 2) - traps(j, 2), 2), 0.5);
      lambda = lambda0*exp(-pow(d, 2)/(2*pow(sigma, 2)));
      if (nz(i) == 1){
	out += kgf_pois(t, lambda);
      } else {
	out += dpois_stable(0, lambda, 1);
      }
    }
  }
  return out;
}


#endif
