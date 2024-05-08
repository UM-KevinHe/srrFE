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
