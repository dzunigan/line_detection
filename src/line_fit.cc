
#include <cmath>
#include <iostream>

// eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "polynomial.h"

#include "util/logging.h"
#include "util/random.h"
#include "util/string.h"

// lambda^3
//+ (m00 + m11 + m22)*lambda^2
//+ (- m01^2 - m02^2 - m12^2 + m00*m11 + m00*m22 + m11*m22)*lambda
//- m22*m01^2 + 2*m01*m02*m12 - m11*m02^2 - m00*m12^2 + m00*m11*m22

// Solves the right nullspace from QR decomposition,
// returning the size of the kernel
template<typename Derived, typename Scalar>
int solveNullspace(const Eigen::MatrixBase<Derived> &A, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &k) {
    Eigen::ColPivHouseholderQR<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>> qr(A.transpose());
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ();

    int n = qr.dimensionOfKernel();
    k.resize(Q.rows(), n);

    k = Q.block(0, Q.cols() - n, Q.rows(), n);
    return qr.dimensionOfKernel();
}

Eigen::Vector3d orthogonal(const Eigen::Vector3d &v) {
  int k = 0;
  if (v[1] < v[k]) k = 1;
  if (v[2] < v[k]) k = 2;
  Eigen::Vector3d e = Eigen::Vector3d::Zero();
  e[k] = 1.;
  return v.cross(e).normalized();
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

int main(int argc, char *argv[]) {
  // Initialize Google's logging library.
  InitializeGlog(argv);
  FLAGS_logtostderr = 1;
  
  Eigen::Vector3d axis(RandomReal(-1., 1.), RandomReal(-1., 1.), RandomReal(-1., 1.));
  axis.normalize();
  
  Eigen::Vector3d o = orthogonal(axis);
  
  Eigen::Matrix3d M;
  M.setZero();
  
  int N = 10;
  for (int i = 0; i < N; ++i) {
    const double theta = RandomReal(-EIGEN_PI, EIGEN_PI);
    Eigen::Vector3d v = o*std::cos(theta)
                        + axis.cross(o)*std::sin(theta)
                        + axis*axis.dot(o)*(1.-std::cos(theta));
    M += v*v.transpose();
  }
  
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
  LOG(INFO) << "Number of candidate solutions: " << lambdas.size();
  
  bool solved = false;
  Eigen::Vector3d x; // solution vector
  double min_cost = std::numeric_limits<double>::max();
  for (Eigen::VectorXd::Index i = 0; i < lambdas.size(); ++i) {
    const double lambda = lambdas[i];
    
    Eigen::MatrixXd candidate_solution;
    int kernel_size = solveNullspace(M + lambda*Eigen::Matrix3d::Identity(), candidate_solution);
    LOG(INFO) << "Kernel size: " << kernel_size;
    
    if (kernel_size != 1) continue;
    
    const double cost = candidate_solution.col(0).transpose()*M*candidate_solution.col(0);
    if (cost < min_cost) {
      x = candidate_solution.col(0);
      min_cost = cost;
      solved = true;
    }
  }
  
  if (!solved) {
    LOG(ERROR) << "Solution not found!";
    return 1;
  }
  
  std::cout << "groundtruth axis: " << axis.transpose() << std::endl;
  const double cost = axis.transpose()*M*axis;
  LOG(INFO) << "Cost: " << cost;
  
  std::cout << "estimated axis: " << x.transpose() << std::endl;
  LOG(INFO) << "Cost: " << min_cost;
  
  return 0;
}
