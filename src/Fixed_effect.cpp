//#define ARMA_NO_DEBUG
#define STRICT_R_HEADERS // needed on Windows, not on macOS
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <omp.h>
#include <mach/mach.h>



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace arma;

//` @importFrom RcppParallel RcppParallelLibs

arma::vec rep(arma::vec &x, arma::vec &each) {
  arma::vec x_rep(sum(each));
  int ind = 0, m = x.n_elem;
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

struct Info_beta : public Worker {
  const arma::mat input1;
  const arma::vec input2;
  arma::mat output;

  Info_beta(const arma::mat &Z, const arma::vec &pq, arma::mat matoutput) :
    input1(Z), input2(pq), output(matoutput) {}

  void operator()(std::size_t begin, std::size_t end) {
    int p = input1.n_cols;
    for (unsigned int i = begin; i < end; i++) {
      div_t divresult = div(int(i), p);
      output(divresult.quot, divresult.rem) = dot(input1.col(divresult.quot), input1.col(divresult.rem)%input2);
    }
  }
};

arma::mat info_beta_tbb(const arma::mat &Z, const arma::vec &pq) {
  arma::mat output(Z.n_cols, Z.n_cols);
  Info_beta info(Z, pq, output);
  parallelFor(0, Z.n_cols*Z.n_cols, info);
  return info.output;
}

void ind2uppsub(unsigned int index, unsigned int dim, unsigned int &row, unsigned int &col) {
  row = 0, col = dim-1;
  unsigned int n = dim*(dim-1)/2 - (dim-row)*(dim-row-1)/2 + col;
  while (index > n) {
    ++row;
    n = dim*(dim-1)/2 - (dim-row)*(dim-row-1)/2 + col;
  }
  while (index < n) {
    --col;
    --n;
  }
}

arma::mat info_beta_omp(const arma::mat &Z, const vec &pq, const int &threads) {
  omp_set_num_threads(threads);
  unsigned int p = Z.n_cols;
  unsigned int loops = p * (1 + p) / 2;
  arma::mat output(p, p);
  #pragma omp parallel for schedule(static)
  for (unsigned int i = 0; i < loops; i++) {
    unsigned int r, c;
    ind2uppsub(i, p, r, c);
    output(r,c) = dot(Z.col(r), Z.col(c)%pq);
    output(c,r) = output(r,c);
  }
  return(output);
}

double Loglkd(const arma::vec &Y, const arma::vec &Z_beta, const arma::vec &gamma_obs) {
  return sum((gamma_obs+Z_beta)%Y-log(1+exp(gamma_obs+Z_beta)));
}


// [[Rcpp::export]]
void checkMemoryUsage() {
  mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;

  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                (task_info_t)&info, &infoCount) != KERN_SUCCESS) {
    Rcpp::Rcout << "Failed to get memory info." << std::endl;
    return;
  }

  Rcpp::Rcout << "Memory used: " << info.resident_size << " bytes" << std::endl;

  arma::vec testvec = arma::zeros<arma::vec>(100000);
  for (int i=0; i<100000; i++){
    testvec(i) = 1;
  }

  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                (task_info_t)&info, &infoCount) != KERN_SUCCESS) {
    Rcpp::Rcout << "Failed to get memory info." << std::endl;
    return;
  }

  Rcpp::Rcout << "Memory used: " << info.resident_size << " bytes" << std::endl;

}


// [[Rcpp::export]]
List logis_fe_prov(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta, int backtrack=1,
                   int max_iter=10000, double bound=10.0, double tol=1e-5, bool message = true) {

  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  arma::vec gamma_obs(n);
  double crit = 100.0;

  if (message == true) {
    cout << "Implementing BAN algorithm (Rcpp) for fixed provider effects model ..." << endl;
  }

  //double meanratio = 0.0;
  if (backtrack==1) {
    double loglkd, d_loglkd, v, lambda, s = 0.01, t = 0.6;
    arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);
    while (iter < max_iter) {
      if (crit < tol) {
        break;
      }
      iter++;
      // provider effect update
      gamma_obs = rep(gamma, n_prov);
      arma::vec Z_beta = Z * beta;
      arma::vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
      arma::vec Yp = Y - p, pq = p % (1-p);
      arma::vec score_gamma(m), d_gamma(m);
      ind = 0;
      for (int i = 0; i < m; i++) {
        score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1)));
        d_gamma(i) = score_gamma(i) / sum(pq(span(ind,ind+n_prov(i)-1)));
        ind += n_prov(i);
      }
      v = 1.0; // initialize step size
      loglkd = Loglkd(Y, Z_beta, gamma_obs);
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep(gamma_tmp, n_prov);
      d_loglkd = Loglkd(Y, Z_beta, gamma_obs_tmp) - loglkd;
      lambda = dot(score_gamma, d_gamma);
      while (d_loglkd < s*v*lambda) {
        v = t*v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep(gamma_tmp, n_prov);
        d_loglkd = Loglkd(Y, Z_beta, gamma_obs_tmp) - loglkd;
      }
      gamma += v * d_gamma;
      gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
      gamma_obs = rep(gamma, n_prov);

      // regression parameter update
      p = 1/(1+exp(-gamma_obs-Z_beta)); // update p
      pq = p % (1-p);
      arma::vec score_beta = Z.t() * (Y-p);
      // arma::mat info_beta = info_beta_tbb(Z, pq); // tbb
      //t2 = clock();
      arma::mat info_beta = Z.t() * (Z.each_col()%pq); // serial
      //t3 = clock();
      arma::vec d_beta = solve(info_beta, score_beta, solve_opts::fast+solve_opts::likely_sympd);
      v = 1.0; // initialize step size
      loglkd = Loglkd(Y, Z_beta, gamma_obs);
      beta_tmp = beta + v * d_beta;
      d_loglkd = Loglkd(Y, Z*beta_tmp, gamma_obs) - loglkd;
      lambda = dot(score_beta, d_beta);
      while (d_loglkd < s*v*lambda) {
        v = t * v;
        beta_tmp = beta + v * d_beta;
        d_loglkd = Loglkd(Y, Z*beta_tmp, gamma_obs) - loglkd;
      }
      beta += v * d_beta;
      //t4 = clock();
      crit = norm(v*d_beta, "inf");
      if (message == true) {
        cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << setprecision(3) << scientific << crit << ";" << endl;
      }
    }
  } else if (backtrack==0) {
    while (iter < max_iter) {
      if (crit < tol) {
        break;
      }
      iter++;
      // provider effect update
      gamma_obs = rep(gamma, n_prov);
      arma::vec Z_beta = Z * beta;
      arma::vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
      arma::vec Yp = Y - p, pq = p % (1-p);
      ind = 0;
      for (int i = 0; i < m; i++) {
        gamma(i) += sum(Yp(span(ind,ind+n_prov(i)-1))) /
          sum(pq(span(ind,ind+n_prov(i)-1)));
        ind += n_prov(i);
      }
      gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
      gamma_obs = rep(gamma, n_prov);
      // regression parameter update
      p = 1/(1+exp(-gamma_obs-Z_beta)); // update p
      pq = p % (1-p);
      arma::vec score_beta = Z.t() * Yp;
      // info_beta = info_beta_tbb(Z, pq); // tbb
      arma::mat info_beta = Z.t() * (Z.each_col()%pq); // serial
      arma::vec d_beta = solve(info_beta, score_beta, solve_opts::fast+solve_opts::likely_sympd);
      beta += d_beta;
      crit = norm(d_beta, "inf");
      if (message == true) {
        cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << setprecision(3) << scientific << crit << ";" << endl;
      }
    }
  }
  if (message == true) {
    cout << "BAN algorithm (Rcpp) converged after " << iter << " iterations!" << endl;
  }
  //cout << "mean ratio is " << fixed << setprecision(7) << meanratio / iter << endl;
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  return ret;
}

// [[Rcpp::export]]
List logis_BIN_fe_prov(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta,
                       int parallel=1, int threads=1, double tol=1e-8, int max_iter=10000,
                       double bound=10.0, bool message = true, bool backtrack = false, const std::string stop = "beta") {

  Rcout << "start" << endl;

  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  double v;
  arma::vec gamma_obs(n);
  double crit = 100.0;
  if (message == true) {
    Rcout << "Implementing SerBIN algorithm (Rcpp) for fixed provider effects model ..." << endl;
  }

  double s = 0.01, t = 0.6; //only used for "backtrack = true"
  double lambda, d_loglkd, loglkd, loglkd_init;
  arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);

  while (iter <= max_iter) {
    if (crit < tol) {
      break;
    }
    iter++;
    gamma_obs = rep(gamma, n_prov);
    arma::vec Z_beta = Z * beta;
    arma::vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
    arma::vec Yp = Y - p, pq = p % (1-p);
    arma::vec score_gamma(m), info_gamma_inv(m);
    arma::mat info_betagamma(Z.n_cols,m);
    ind = 0;
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1)));
      info_gamma_inv(i) = 1 / sum(pq(span(ind,ind+n_prov(i)-1)));
      info_betagamma.col(i) =
        sum(Z.rows(ind,ind+n_prov(i)-1).each_col()%(p.subvec(ind,ind+n_prov(i)-1)%(1-p.subvec(ind,ind+n_prov(i)-1)))).t();
      ind += n_prov(i);
    }

    arma::vec score_beta = Z.t() * Yp;
    arma::mat info_beta(Z.n_cols, Z.n_cols);
    if (parallel==1) { // parallel
      info_beta = info_beta_omp(Z, pq, threads); // omp
      // info_beta = info_beta_tbb(Z, pq); // tbb
    } else if (parallel==0) { // serial
      info_beta = Z.t() * (Z.each_col()%pq);
    }
    arma::mat mat_tmp1 = trans(info_betagamma.each_row()%info_gamma_inv.t());
    arma::mat schur_inv = inv_sympd(info_beta-mat_tmp1.t()*info_betagamma.t());
    arma::mat mat_tmp2 = mat_tmp1*schur_inv;
    arma::vec d_gamma = info_gamma_inv%score_gamma + mat_tmp2*(mat_tmp1.t()*score_gamma-score_beta);
    arma::vec d_beta = schur_inv*score_beta - mat_tmp2.t()*score_gamma;

    // List ret = List::create(_["score_gamma"]=score_gamma,
    //                         _["info_gamma_inv"]=info_gamma_inv,
    //                         _["d_beta"]=d_beta);
    // return ret;

    v = 1.0; // initialize step size
    if (backtrack == true){
      loglkd = Loglkd(Y, Z * beta, rep(gamma, n_prov));
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep(gamma_tmp, n_prov);
      arma::vec Z_beta_tmp = Z * (beta+v*d_beta);
      lambda = dot(score_gamma, d_gamma) + dot(score_beta, d_beta);
      while (d_loglkd < s*v*lambda) {
        v = t*v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep(gamma_tmp, n_prov);
        Z_beta_tmp = Z * (beta+v*d_beta);
        d_loglkd = Loglkd(Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
      }
    }
    if (iter == 1) loglkd_init = loglkd;
    gamma += v * d_gamma;
    gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
    beta += v * d_beta;
    if (stop == "beta"){
      crit = norm(v*d_beta, "inf");
    }
    else if (stop == "relch") {
      crit = abs(d_loglkd/(d_loglkd+loglkd));
    }
    else if (stop == "ratch") {
      crit = abs(d_loglkd/(d_loglkd+loglkd-loglkd_init));
    }
    else if (stop == "all") {
      arma::vec crits(3);
      crits(0) = norm(v*d_beta, "inf");
      crits(1) = abs(d_loglkd/(d_loglkd+loglkd));
      crits(2) = abs(d_loglkd/(d_loglkd+loglkd-loglkd_init));
      crit = crits.max();
    }

    if (message == true) {
      cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << scientific << setprecision(3) << crit << ";";
    }
  }
  if (message == true) {
    cout << "serBIN (Rcpp) algorithm converged after " << iter << " iterations!" << endl;
  }
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  return ret;
}


// [[Rcpp::export]]
List logis_firth_prov(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta, int n_obs,
                      int m, int parallel=1, int threads=1, double tol=1e-8, int max_iter=10000,
                      double bound=10.0, bool message = true, bool backtrack = false) {

  int iter = 0, n = n_obs, ind;
  double v;
  arma::vec gamma_obs(n);
  double crit = 100.0;
  if (message == true) {
    cout << "Implementing firth-corrected fixed provider effects model (Rcpp) ..." << endl;
  }

  double s = 0.01, t = 0.6; //only used for "backtrack = true"
  double lambda, d_loglkd, loglkd;
  arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);

  while (iter < max_iter) {
    if (crit < tol) {
      break;
    }
    iter++;
    gamma_obs = rep(gamma, n_prov);
    arma::vec Z_beta = Z * beta;
    arma::vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
    arma::vec Yp = Y - p, pq = p % (1-p);
    arma::vec score_gamma(m), info_gamma_inv(m);
    arma::mat info_betagamma(Z.n_cols,m);

    ind = 0;
    for (int i = 0; i < m; i++) {
      info_gamma_inv(i) = 1 / sum(pq(span(ind,ind+n_prov(i)-1)));
      info_betagamma.col(i) =
        sum(Z.rows(ind,ind+n_prov(i)-1).each_col()%(p.subvec(ind,ind+n_prov(i)-1)%(1-p.subvec(ind,ind+n_prov(i)-1)))).t();
      ind += n_prov(i);
    }
    arma::mat info_beta(Z.n_cols, Z.n_cols);
    if (parallel==1) { // parallel
      info_beta = info_beta_omp(Z, pq, threads); // omp
      // info_beta = info_beta_tbb(Z, pq); // tbb
    } else if (parallel==0) { // serial
      info_beta = Z.t() * (Z.each_col()%pq);
    }

    arma::mat mat_tmp1 = trans(info_betagamma.each_row()%info_gamma_inv.t());
    arma::mat schur_inv = inv_sympd(info_beta-mat_tmp1.t()*info_betagamma.t());
    arma::mat mat_tmp2 = mat_tmp1*schur_inv;

    arma::mat prod = mat_tmp1 * schur_inv * mat_tmp1.t();
    arma::vec diag_prod = info_gamma_inv + prod.diag();
    arma::vec c1 = rep(diag_prod, n_prov);

    arma::vec c2(n);
    ind = 0;
    for (int i = 0; i < m; i++) {
      c2.subvec(ind,ind+n_prov(i)-1) = - Z.rows(ind,ind+n_prov(i)-1) * mat_tmp2.t().eval().col(i);
      ind += n_prov(i);
    }

    arma::vec c3(n);
    for (int i = 0; i < n; i++) {
      double a = (Z.row(i) * schur_inv * Z.row(i).t()).eval()(0,0);
      c3(i) = a;
    }

    arma::vec YpA = Yp + pq % (c1 + c2 + c2 + c3) % (0.5 - p);
    arma::vec score_beta = Z.t() * YpA;
    ind = 0;
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(YpA(span(ind,ind+n_prov(i)-1)));
      ind += n_prov(i);
    }

    arma::vec d_gamma = info_gamma_inv%score_gamma + mat_tmp2*(mat_tmp1.t()*score_gamma-score_beta);
    arma::vec d_beta = schur_inv*score_beta - mat_tmp2.t()*score_gamma;
    v = 1.0; // initialize step size
    if (backtrack == true){
      loglkd = Loglkd(Y, Z * beta, rep(gamma, n_prov));
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep(gamma_tmp, n_prov);
      arma::vec Z_beta_tmp = Z * (beta+v*d_beta);
      d_loglkd = Loglkd(Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
      lambda = dot(score_gamma, d_gamma) + dot(score_beta, d_beta);
      while (d_loglkd < s*v*lambda) {
        v = t*v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep(gamma_tmp, n_prov);
        Z_beta_tmp = Z * (beta+v*d_beta);
        d_loglkd = Loglkd(Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
      }
    }
    gamma += v * d_gamma;
    gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
    beta += v * d_beta;
    loglkd += d_loglkd;
    crit = norm(v*d_beta, "inf");

    if (message == true) {
      cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << scientific << setprecision(3) << crit << ";";
    }
  }
  if (message == true) {
    cout << "Algorithm converged after " << iter << " iterations!" << endl;
  }
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  return ret;
}


// [[Rcpp::export]]
List wald_covar(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma, arma::vec &beta, arma::uvec &indices, double null, double alpha) {

  indices -= 1; // switch to C indexing
  arma::vec gamma_obs = rep(gamma, n_prov);
  int m = n_prov.n_elem;
  arma::vec p = 1 / (1 + exp(-gamma_obs-Z*beta));
  p = clamp(p, 1e-10, 1-1e-10);
  int ind = 0;
  arma::vec info_gamma_inv(m);
  arma::mat info_betagamma(Z.n_cols,m);
  for (int i = 0; i < m; i++) {
    info_gamma_inv(i) = 1 / dot(p.subvec(ind,ind+n_prov(i)-1),1-p.subvec(ind,ind+n_prov(i)-1));
    info_betagamma.col(i) =
      sum(Z.rows(ind,ind+n_prov(i)-1).each_col()%(p.subvec(ind,ind+n_prov(i)-1)%(1-p.subvec(ind,ind+n_prov(i)-1)))).t();
    ind += n_prov(i);
  }
  arma::mat info_beta = Z.t()*(Z.each_col()%(p%(1-p)));
  arma::mat info_beta_inv = inv_sympd(info_beta-(info_betagamma.each_row()%info_gamma_inv.t())*info_betagamma.t());
  arma::vec se_beta = sqrt(info_beta_inv.diag());
  arma::vec stat = (beta(indices)-null)/se_beta(indices);
  arma::vec p_lower = normcdf(stat);
  List ret;
  ret["stat"] = stat;
  ret["p"] = 2*min(p_lower,1-p_lower);
  ret["se.beta"] = 1*se_beta(indices);
  ret["beta.lower"] = beta(indices) - R::qnorm5(1-alpha/2,0,1,true,false)*se_beta(indices);
  ret["beta.upper"] = beta(indices) + R::qnorm5(1-alpha/2,0,1,true,false)*se_beta(indices);
  return ret;
}





// double Exp_direct(double gamma, const arma::vec& Z_beta) {
//   arma::vec logist_vals = 1.0 / (1.0 + exp(-(gamma + Z_beta)));
//   return arma::accu(logist_vals);
// }
//
// // [[Rcpp::export]]
// arma::vec computeDirectExp(const arma::vec& gamma_prov, const arma::vec& Z_beta, const int &threads) {
//   omp_set_num_threads(threads);
//   arma::vec results(gamma_prov.n_elem);
//
//   #pragma omp parallel for schedule(static)
//   for (unsigned int i = 0; i < gamma_prov.n_elem; ++i) {
//     results[i] = Exp_direct(gamma_prov[i], Z_beta);
//   }
//
//   return results;
// }
