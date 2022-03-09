# LOG LIKELIHOOD FUNCTIONS

# negative log maxlike function for the normal distribution
loglike_norm <- function(mu_val, var_val, x_vec){
  x_n = length(x_vec)
  n_ll <- x_n/2*log(2*pi)+x_n/2*log(var_val)+1/(2*var_val)*sum((x_vec-mu_val)^2)
  return(n_ll)
}

# first derivative negative log maxlike function for the normal distribution
loglike_norm_1d <- function(mu_val, var_val, x_vec){
  x_n = length(x_vec)
  n_ll_1d <- -(1/var_val)*(sum(x_vec)-x_n*mu_val)
  return(n_ll_1d)
}

# second derivative negative log maxlike function for the normal distribution
loglike_norm_2d <- function(mu_val, var_val, x_vec){
  x_n = length(x_vec)
  n_ll_2d <- x_n/var_val
  return(n_ll_2d)
}

# negative loglike function for the cauchy distribution
loglike_cauy <- function(loc_val, gam_val, x_vec){
  x_n = length(x_vec)
  n_ll <- -x_n*log(gam_val) +x_n*log(pi) +sum(log(gam_val^2+(x_vec-loc_val)^2))
  ll <- x_n*log(gam_val) -x_n*log(pi) -sum(log(gam_val^2+(x_vec-loc_val)^2))
  return(n_ll)
}

# loglike function for the cauchy distribution
ploglike_cauy <- function(loc_val, gam_val, x_vec){
  x_n = length(x_vec)
  ll <- x_n*log(gam_val) -x_n*log(pi) -sum(log(gam_val^2+(x_vec-loc_val)^2))
  return(ll)
}

# first derivative negative log maxlike function for the cauchy distribution
loglike_cauy_1d <- function(loc_val, gam_val, x_vec){
  x_n = length(x_vec)
  n_ll_1d <- -sum(2*(x_vec-loc_val)/(gam_val^2+(x_vec-loc_val)^2))
  return(n_ll_1d)
}

# second derivative negative log maxlike function for the cauchy distribution
loglike_cauy_2d <- function(loc_val, gam_val, x_vec){
  x_n = length(x_vec)
  n_ll_2d <- -sum(2*(loc_val-gam_val-x_vec)*(loc_val+gam_val-x_vec)
                 /((loc_val-x_vec)^2 + gam_val^2)^2)
  return(n_ll_2d)
}
