## Create the negative log posterior (U) for use with TMB

U.EB <- "
#include <TMB.hpp>
#include <vector>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR( y_i );  // counts for observation i
  DATA_VECTOR( N_i );  // number of people per county
  DATA_VECTOR( Z_j );  // pollution vector for observation i
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_SPARSE_MATRIX(A);
  
  // Parameters
  PARAMETER( beta_0 );
  PARAMETER( beta_1 );
  PARAMETER_VECTOR( theta );
  
  // Random effects
  PARAMETER_VECTOR( S_j );
  

  // Objective function
  int n_i = y_i.size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero(); 
  
  // Probability of thetas
  // Flat prior on thetas
  
  // Probability of intercept
  // Flat prior on beta_0
  
  // Probability of beta_1
  // Flat prior on beta_1
  
  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(2*theta(0)) * ( exp(4*theta(1))*M0 + Type(2.0)*exp(2*theta(1))*M1 + M2);
  jnll_comp(1) += GMRF(Q) ( S_j );  
  
  // Probability of data conditional on random effects
  vector<Type> lp = S_j + beta_1 * Z_j;
  vector<Type> S_j_star = A * exp(lp);
  for( int i=0; i<n_i; i++){
  jnll_comp(0) -= dpois( y_i(i), N_i(i) * exp(beta_0) * S_j_star(i), true );
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  return jnll;
}
"

setwd("### set working directory")
write(U.EB, file = "TMB/U_EB.cpp")
compile("TMB/U_EB.cpp")
dyn.load(dynlib("TMB/U_EB"))


