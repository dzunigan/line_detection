
#define PROGRAM_NAME \
  "evalate_tum"

#define FLAGS_CASES                                                                                \
  FLAG_CASE(double, min_line, 10.0, "Minimum line length [%]")                                     \
  FLAG_CASE(string, output_file, "./results.csv", "Evaluation results csv save file")

#define ARGS_CASES                                                                                 \

#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

// eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>

// nanoflann
#include "nanoflann/nanoflann.hpp"

#include "camera_models.h"
#include "FOVUndistorter.h"

#include "dlsd.h"
#include "EDLines.h"
#include "thirdparty/LSD/include/lsd.h"
#include "thirdparty/AG3line/include/ag3line.h"

#include "sequence.hpp"

#include "util/args.hpp"
#include "util/csv.hpp"
#include "util/logging.h"
#include "util/math.h"
#include "util/random.h"
#include "util/timer.h"

// Parameters
#define SEARCH_RADIUS 4.0 // px
#define MIN_OVERLAP 90.0 // %

// Debug only
#include <opencv2/highgui.hpp>

struct kdtree_data_t
{
  std::vector<cv::Point> points;
  std::vector<int> segments;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return points.size(); }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if (dim == 0) return points[idx].x;
    else return points[idx].y;
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

};

// construct a kd-tree index:
using kdtree_t = nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, kdtree_data_t>,
  kdtree_data_t,
  2>;

void ValidateArgs() { }

void ValidateFlags() {
  CHECK_LE(FLAGS_min_line, 100.0) << "Argument out of percent bounds";
  CHECK_GE(FLAGS_min_line,   0.0) << "Argument out of percent bounds";
}

template <typename CameraModel>
double computeAnglePinhole(const double *params,
                           const double x1, const double y1,
                           const double x2, const double y2) {
  Eigen::Vector3d v1;
  CameraModel::ImageToWorld(params, x1, y1, &v1.x(), &v1.y());
  v1.z() = 1.0;
  v1.normalize();

  Eigen::Vector3d v2;
  CameraModel::ImageToWorld(params, x2, y2, &v2.x(), &v2.y());
  v2.z() = 1.0;
  v2.normalize();

  return std::acos(v1.dot(v2));
}

template double computeAnglePinhole<PinholeCameraModel>(const double *params, const double x1, const double y1, const double x2, const double y2);
template double computeAnglePinhole<FOVCameraModel>(const double *params, const double x1, const double y1, const double x2, const double y2);

double computeMaxViewAngle(const double *params, int width, int height) {
  const double cx = params[2];
  const double cy = params[3];

  double angle, min_angle = computeAnglePinhole<PinholeCameraModel>(params, cx, 0, cx, cy);

  angle = computeAnglePinhole<PinholeCameraModel>(params, width-1, cy, cx, cy);
  if (angle < min_angle) min_angle = angle;

  angle = computeAnglePinhole<PinholeCameraModel>(params, cx, height-1, cx, cy);
  if (angle < min_angle) min_angle = angle;

  angle = computeAnglePinhole<PinholeCameraModel>(params, 0, cy, cx, cy);
  if (angle < min_angle) return angle;

  return 2.*min_angle;
}

int main(int argc, char *argv[]) {
  (void) argc; // unused

  SetPRNGSeed((unsigned) std::time(nullptr));

  // Handle help flag
  if (args::HelpRequired(argc, argv)) {
    args::ShowHelp();
    return 0;
  }

  // Parse input flags
  args::ParseCommandLineNonHelpFlags(&argc, &argv, true);

  // Initialize Google's logging library.
  InitializeGlog(argv);
  FLAGS_logtostderr = 1;

  // Check number of args
  if (argc-1 != args::NumArgs()) {
    args::ShowHelp();
    return -1;
  }

  // Parse input args
  args::ParseCommandLineArgs(argc, argv);

  // Validate input arguments
  ValidateArgs();
  ValidateFlags();

  std::shared_ptr<UndistorterFOV> undistorter = std::make_shared<UndistorterFOV>("./data/camera.txt");

  // https://stackoverflow.com/a/13062069
  std::shared_ptr<double[]> params(new double[5]);
  params[0] = undistorter->getOriginalCalibration()[0];
  params[1] = undistorter->getOriginalCalibration()[1];
  params[2] = undistorter->getOriginalCalibration()[2];
  params[3] = undistorter->getOriginalCalibration()[3];
  params[4] = undistorter->getOriginalCalibration()[4];

  std::shared_ptr<double[]> rectified_params(new double[4]);
  rectified_params[0] = undistorter->getRectifiedCalibration()[0];
  rectified_params[1] = undistorter->getRectifiedCalibration()[1];
  rectified_params[2] = undistorter->getRectifiedCalibration()[2];
  rectified_params[3] = undistorter->getRectifiedCalibration()[3];

  //cv::Point2d center(rectified_params[2], rectified_params[3]);

  // image sizes
  const int widthOrg = undistorter->getInputDims()[0];
  const int heightOrg = undistorter->getInputDims()[1];
  const int width = undistorter->getOutputDims()[0];
  const int height = undistorter->getOutputDims()[1];

  float* internalBuffer = new float[widthOrg*heightOrg];

  // compute evaluation parameters
  const double max_view_angle = computeMaxViewAngle(rectified_params.get(), width, height);
  const double min_line = FLAGS_min_line / 100.0 * max_view_angle;

  LOG(INFO) << StringPrintf("Max view angle: %.2fº", RadToDeg(max_view_angle));
  LOG(INFO) << StringPrintf("Min. line: %.2f%%", FLAGS_min_line);

  // Initialize VideoWriter
  /*
  cv::VideoWriter writer;
  int codec = cv::VideoWriter::fourcc('a', 'v', 'c', '1');
  cv::Size frameSize(width, height);
  writer.open("./whitePaper.mp4", codec, 10.0, frameSize, true);
  */

  Timer timer;
  std::vector<double> rectification_times;
  std::vector<double> lsd_times, ag3line_times, edlines_times;
  std::vector<double> dlsd_times;

  std::vector<std::string> files = getSequence<unsigned>("./datasets/tum-mono/calib_wide_whitePaper/images/");
  Eigen::MatrixXd results(files.size(), 11); // lsd, lsd_min_len, ag3lines, ag3line_min_len, edlines, edlines_min_len, dlsd, dlsd_min_len, lsd_matched, ag3lines_matched, edlines_matched
  results.setZero();

  cv::Ptr<LineSegmentDetector> lsd = createLineSegmentDetector(LSD_REFINE_NONE);

  // kdtree search params
  nanoflann::SearchParams search_params;
  search_params.sorted = true;

  int row = 0;
  for (const std::string &file : files) {
    std::cout << "Processing: " << file << '\r' << std::flush;

    cv::Mat imageRaw;
    imageRaw = cv::imread(file, cv::IMREAD_GRAYSCALE);

    timer.Start();
    undistorter->undistort<unsigned char>(imageRaw.data, internalBuffer, widthOrg*heightOrg, width*height);
    rectification_times.push_back(1e-3*timer.ElapsedMicroSeconds());

    cv::Mat imageRectified(height, width, CV_32F, internalBuffer);

    cv::Mat image;
    imageRectified.convertTo(image, CV_8U);

    // LSD
    timer.Start();
    std::vector<cv::Vec4f> lsdLines;
    lsd->detect(image, lsdLines);
    lsd_times.push_back(1e-3*timer.ElapsedMicroSeconds());
    results(row, 0) = lsdLines.size();

    // AG3line
    timer.Start();
    std::vector<lineag> ag3Lines;
    ag3line(imageRectified, ag3Lines, true);
    ag3line_times.push_back(1e-3*timer.ElapsedMicroSeconds());
    results(row, 2) = ag3Lines.size();

    // EDLines
    timer.Start();
    EDLines edlines = EDLines(image, 1.0, -1, 6.0, 1.3);
    edlines_times.push_back(1e-3*timer.ElapsedMicroSeconds());
    results(row, 4) = edlines.getLines().size();

    // D-LSD
    timer.Start();
    std::vector<DistortedLineSegment> distortedLines;
    dlsd<FOVCameraModel>(imageRaw, params.get(), distortedLines);
    dlsd_times.push_back(1e-3*timer.ElapsedMicroSeconds());
    results(row, 6) = distortedLines.size();

// DEBUG ONLY
/*
    cv::Mat imageDLSD;
    cv::cvtColor(imageRaw, imageDLSD, cv::COLOR_GRAY2RGB);

    cv::Mat imageColor;
    cv::cvtColor(image, imageColor, cv::COLOR_GRAY2RGB);
*/

    // KD Tree
    int dlsd_min_len = 0;
    kdtree_data_t segmentedPixels;
    for (const DistortedLineSegment &distortedLine : distortedLines) {
      const cv::Point &s = distortedLine.segment.front();
      const cv::Point &e = distortedLine.segment.back();
      const double line_angle = computeAnglePinhole<FOVCameraModel>(params.get(),
                                                                    s.x, s.y,
                                                                    e.x, e.y);
      if (line_angle >= min_line) {
        //for (const cv::Point &point : distortedLine.segment)
        //  cv::circle(imageDLSD, point, 1, cv::Scalar(0, 200, 0), -1);

        dlsd_min_len++;
        segmentedPixels.points.insert(segmentedPixels.points.end(), distortedLine.segment.begin(), distortedLine.segment.end());
        segmentedPixels.segments.insert(segmentedPixels.segments.end(), distortedLine.segment.size(), distortedLine.segmentNo);
      } //else {
        //for (const cv::Point &point : distortedLine.segment)
        //  cv::circle(imageDLSD, point, 1, cv::Scalar(200, 0, 0), -1);
      //}
    }
    results(row, 7) = dlsd_min_len;

    kdtree_t index(2, segmentedPixels, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    // LSD
    int lsdlines_min_len = 0, lsdlines_matched = 0;
    for (const cv::Vec4f &line : lsdLines) {
      const cv::Point s(line[0], line[1]);
      const cv::Point e(line[2], line[3]);
      const double line_angle = computeAnglePinhole<PinholeCameraModel>(rectified_params.get(),
                                                                        s.x, s.y,
                                                                        e.x, e.y);
      if (line_angle >= min_line) {
        //cv::LineIterator lit_(image, s, e);
        //for (int i = 0; i < lit_.count; i++, ++lit_)
        //  cv::circle(imageColor, lit_.pos(), 1, cv::Scalar(0, 200, 0), -1);

        lsdlines_min_len++;

        int N = 0;
        cv::LineIterator lit(image, s, e);
        std::unordered_map<int, int> segment_histogram;
        while (N < lit.count) {
          const cv::Point point = lit.pos();
          N++;

          double u, v;
          double point_[2];
          PinholeCameraModel::ImageToWorld<double>(rectified_params.get(), point.x, point.y, &u, &v);
          FOVCameraModel::WorldToImage<double>(params.get(), u, v, &point_[0], &point_[1]);

          std::vector<std::pair<size_t, double>> matches;
          index.radiusSearch(&point_[0], SEARCH_RADIUS*SEARCH_RADIUS, matches, search_params);
          if (!matches.empty())
            segment_histogram[segmentedPixels.segments[matches.front().first]]++;
          ++lit;
        }

        int max_count = 0;
        for (const std::pair<int, int> &entry : segment_histogram)
          if (entry.second > max_count) max_count = entry.second;

        const double ratio = 100.0 * static_cast<double>(max_count) / static_cast<double>(N);
        //LOG(INFO) << StringPrintf("Found %d out of %d (%.2f %%)", max_count, N, ratio);

        if (ratio > MIN_OVERLAP) {
          lsdlines_matched++;
          //cv::line(imageColor, s, e, cv::Scalar(0, 200, 0), 2);
        } //else
          //cv::line(imageColor, s, e, cv::Scalar(0, 0, 200), 2);
      } //else
        //cv::line(imageColor, s, e, cv::Scalar(200, 0, 0), 2);
    }
    results(row, 1) = lsdlines_min_len;
    results(row, 8) = lsdlines_matched;

    //cv::imshow("Raw D-LSD", imageDLSD);
    //cv::imshow("Rectified LSD", imageColor);
    //cv::waitKey(0);

    // AG3line
    int ag3lines_min_len = 0, ag3lines_matched = 0;
    for (const lineag &line : ag3Lines) {
      const cv::Point s(line.x1, line.y1);
      const cv::Point e(line.x2, line.y2);

      if (computeAnglePinhole<PinholeCameraModel>(rectified_params.get(),
                                                  s.x, s.y,
                                                  e.x, e.y) > min_line) {
        ag3lines_min_len++;

        int N = 0;
        cv::LineIterator lit(imageRectified, s, e);
        std::unordered_map<int, int> segment_histogram;
        while (N < lit.count) {
          const cv::Point point = lit.pos();
          N++;

          double u, v;
          double point_[2];
          PinholeCameraModel::ImageToWorld<double>(rectified_params.get(), point.x, point.y, &u, &v);
          FOVCameraModel::WorldToImage<double>(params.get(), u, v, &point_[0], &point_[1]);

          std::vector<std::pair<size_t, double>> matches;
          index.radiusSearch(&point_[0], SEARCH_RADIUS*SEARCH_RADIUS, matches, search_params);
          if (!matches.empty())
            segment_histogram[segmentedPixels.segments[matches.front().first]]++;
          ++lit;
        }

        int max_count = 0;
        for (const std::pair<int, int> &entry : segment_histogram)
          if (entry.second > max_count) max_count = entry.second;

        const double ratio = 100.0 * static_cast<double>(max_count) / static_cast<double>(N);
        //LOG(INFO) << StringPrintf("Found %d out of %d (%.2f %%)", max_count, N, ratio);

        if (ratio > MIN_OVERLAP) {
          ag3lines_matched++;
          //cv::line(imageColor, s, e, cv::Scalar(0, 200, 0), 2);
        } //else
          //cv::line(imageColor, s, e, cv::Scalar(0, 0, 200), 2);
      } //else
        //cv::line(imageColor, s, e, cv::Scalar(200, 0, 0), 2);
    }
    results(row, 3) = ag3lines_min_len;
    results(row, 9) = ag3lines_matched;

    // EDLines
    int edlines_min_len = 0, edlines_matched = 0;
    for (const LineSegment &line : edlines.getLines()) {
      const cv::Point s(line.sx, line.sy);
      const cv::Point e(line.ex, line.ey);

      if (computeAnglePinhole<PinholeCameraModel>(rectified_params.get(),
                                                  s.x, s.y,
                                                  e.x, e.y) > min_line) {
        edlines_min_len++;

        int N = 0;
        cv::LineIterator lit(image, s, e);
        std::unordered_map<int, int> segment_histogram;
        while (N < lit.count) {
          const cv::Point point = lit.pos();
          N++;

          double u, v;
          double point_[2];
          PinholeCameraModel::ImageToWorld<double>(rectified_params.get(), point.x, point.y, &u, &v);
          FOVCameraModel::WorldToImage<double>(params.get(), u, v, &point_[0], &point_[1]);

          std::vector<std::pair<size_t, double>> matches;
          index.radiusSearch(&point_[0], SEARCH_RADIUS*SEARCH_RADIUS, matches, search_params);
          if (!matches.empty())
            segment_histogram[segmentedPixels.segments[matches.front().first]]++;
          ++lit;
        }

        int max_count = 0;
        for (const std::pair<int, int> &entry : segment_histogram)
          if (entry.second > max_count) max_count = entry.second;

        const double ratio = 100.0 * static_cast<double>(max_count) / static_cast<double>(N);
        //LOG(INFO) << StringPrintf("Found %d out of %d (%.2f %%)", max_count, N, ratio);

        if (ratio > MIN_OVERLAP) {
          edlines_matched++;
          //cv::line(imageColor, s, e, cv::Scalar(0, 200, 0), 2);
        } //else
          //cv::line(imageColor, s, e, cv::Scalar(0, 0, 200), 2);
      } //else
        //cv::line(imageColor, s, e, cv::Scalar(200, 0, 0), 2);
    }
    results(row, 5) = edlines_min_len;
    results(row, 10) = edlines_matched;

    //writer.write(imageLines);

    row++;
  }
  std::cout << std::endl;

  std::cout << StringPrintf("Rectification: %.3fms (%.3fms)", Mean(rectification_times), StdDev(rectification_times)) << std::endl;
  std::cout << StringPrintf("LSD: %.3fms (%.3fms)", Mean(lsd_times), StdDev(lsd_times)) << std::endl;
  std::cout << StringPrintf("AG3Line: %.3fms (%.3fms)", Mean(ag3line_times), StdDev(ag3line_times)) << std::endl;
  std::cout << StringPrintf("EDLines: %.3fms (%.3fms)", Mean(edlines_times), StdDev(edlines_times)) << std::endl;
  std::cout << StringPrintf("D-LSD: %.3fms (%.3fms)", Mean(dlsd_times), StdDev(dlsd_times)) << std::endl;

  if (!FLAGS_output_file.empty()) {
    csv::write<>(results, FLAGS_output_file, "# lsd, lsd_min_len, ag3lines, ag3line_min_len, edlines, edlines_min_len, dlsd, dlsd_min_len, lsd_matched, ag3lines_matched, edlines_matched");
    std::cout << "Saving file " << FLAGS_output_file << std::endl;
  }

  delete[] internalBuffer;

  return 0;
}
