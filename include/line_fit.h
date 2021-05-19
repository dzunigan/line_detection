
#include <cmath>
#include <iostream>
#include <limits>

// eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "polynomial.h"

#include "util/logging.h"
#include "util/random.h"
#include "util/string.h"

#ifndef DLSD_SOLVER_ZERO_TH
  #define DLSD_SOLVER_ZERO_TH 1e-9
#endif

// Solves the right nullspace from QR decomposition,
// returning the size of the kernel
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> solveNullspace(const Eigen::Matrix3d &A) {
    /*
    Eigen::FullPivHouseholderQR<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>> qr(A.transpose());
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ();

    int n = qr.dimensionOfKernel();
    k.resize(Q.rows(), n);

    k = Q.block(0, Q.cols() - n, Q.rows(), n);
    return qr.dimensionOfKernel();
    */
    
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    lu.setThreshold(DLSD_SOLVER_ZERO_TH);
    return lu.kernel();
}

Eigen::VectorXd real_roots(const Eigen::VectorXd &real, const Eigen::VectorXd &imag) {
  CHECK_EQ(real.size(), imag.size());

  Eigen::VectorXd roots(real.size());

	Eigen::VectorXd::Index j = 0;
	for (Eigen::VectorXd::Index i = 0; i < real.size(); ++i) {
	  if (!imag(i)) {
	    roots(j) = real(i);
	    ++j;
	  }
	}

	roots.conservativeResize(j);
	return roots;
}

bool line_fit(const std::vector<Eigen::Vector3d> &observations, Eigen::Vector3d &axis) {  
  //LOG(INFO) << "Number of observations: " << observations.size();
  if (observations.size() < 3) return false;

  Eigen::Matrix3d M;
  M.setZero();

  for (const Eigen::Vector3d &v : observations)
    M += v*v.transpose();
  
  // lambda^3
  //+ (m00 + m11 + m22)*lambda^2
  //+ (- m01^2 - m02^2 - m12^2 + m00*m11 + m00*m22 + m11*m22)*lambda
  //- m22*m01^2 + 2*m01*m02*m12 - m11*m02^2 - m00*m12^2 + m00*m11*m22
  
  Eigen::VectorXd coeffs(4);
  coeffs << 1.,
            M.trace(),
            - std::pow(M(0,1), 2) - std::pow(M(0,2), 2) - std::pow(M(1,2), 2) + M(0,0)*M(1,1) + M(0,0)*M(2,2) + M(1,1)*M(2,2),
            - M(2,2)*std::pow(M(0,1), 2) + 2.*M(0,1)*M(0,2)*M(1,2) - M(1,1)*std::pow(M(0,2), 2) - M(0,0)*std::pow(M(1,2), 2) + M(0,0)*M(1,1)*M(2,2);
  
  Eigen::VectorXd real, imag;
  if (!FindPolynomialRootsCompanionMatrix(coeffs, &real, &imag)) {
    LOG(ERROR) << "Failed to find roots\n"
               << StringPrintf("%.16f %.16f %.16f %.16f",
                               coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
    return 1;
  }

  Eigen::VectorXd lambdas = real_roots(real, imag);
  if (lambdas.size() == 0) {
    LOG(ERROR) << "No real roots found\n"
               << StringPrintf("%.16f %.16f %.16f %.16f",
                               coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
    return 1;
  }
  //LOG(INFO) << "Number of candidate solutions: " << lambdas.size();
  
  bool solved;
  Eigen::Vector3d x; // solution vector
  double min_cost = std::numeric_limits<double>::max();
  for (Eigen::VectorXd::Index i = 0; i < lambdas.size(); ++i) {
    const double lambda = lambdas[i];

    Eigen::MatrixXd kernel = solveNullspace(M + lambda*Eigen::Matrix3d::Identity());
    //LOG(INFO) << "Kernel size: " << kernel.cols();
    
    // The size of the kernel must be 1 by construction
    if (kernel.cols() != 1) {
      LOG(WARNING) << "Kernel must be 1D!";
      //continue;
      solved = false;
      break;
    }
    Eigen::Vector3d candidate_solution = kernel.col(0);
    candidate_solution.normalize();
    
    //LOG(INFO) << "Candidate solution norm: " << candidate_solution.norm();
    if (std::abs(candidate_solution.norm() - 1.) > DLSD_SOLVER_ZERO_TH) {
      //continue;
      LOG(WARNING) << "Candidate solution norm: " << candidate_solution.norm();
      solved = false;
      break;
    }
    
    const double cost = candidate_solution.transpose()*M*candidate_solution;
    //LOG(INFO) << "Candidate solution cost: " << cost;
    if (cost < min_cost) {
      x = candidate_solution;
      min_cost = cost;
      solved = true;
    }
  }
  
  if (!solved) return false;
  
  //LOG(INFO) << "Cost: " << min_cost;
  // TODO Threshold cost?
  
  axis = x;
  return true;
}
