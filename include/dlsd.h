
#ifndef DLSD_H_
#define DLSD_H_

// stl
#include <unordered_set>
#include <vector>

// eigen
#include <Eigen/Core>

// opencv
#include <opencv2/core/mat.hpp>
#include <opencv2/core/types.hpp>

// Parameters
//  Line Segments
#define MIN_LINE_LENGTH 9 // px
#define MAX_LINE_ERROR 1.0 // px
//  Polygonal Chains
#define SMOOTHNESS_TH 0.19634954084936207 // rad (M_PI/16)
#define REPROJECTION_TH 4.0 // px
#define MIN_SINGLE_LINE_LENGTH 25 // px

struct DistortedLineSegment {
  int segmentNo; // Edge segment that this distorted line segment belongs to

  Eigen::Vector3d axis;
  std::vector<cv::Point> segment;
  std::unordered_set<int> invalidPixels;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename CameraModel>
void dlsd(const cv::Mat &image, const double *params, std::vector<DistortedLineSegment> &distortedLineSegments);

#endif  // DLSD_H_
