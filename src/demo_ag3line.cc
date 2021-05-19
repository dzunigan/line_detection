
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// eigen
#include <Eigen/Core>

// gflags
#include <gflags/gflags.h>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "FOVUndistorter.h"

#include "thirdparty/AG3line/include/ag3line.h"

#include "sequence.hpp"

#include "util/alignment.h"
#include "util/logging.h"
#include "util/math.h"
#include "util/string.h"
#include "util/timer.h"

#define MIN_SINGLE_LINE_LENGTH 100 // px

int main(int argc, char *argv[]) {
  (void) argc; // unused

  std::vector<double> rectification;
  std::vector<double> line_extraction;

  // Initialize Google's logging library.
  InitializeGlog(argv);
  FLAGS_logtostderr = 1;

  std::shared_ptr<UndistorterFOV> undistorter = std::make_shared<UndistorterFOV>("./data/camera.txt");

  // get image widths.
  int widthOrg = undistorter->getInputDims()[0];
  int heightOrg = undistorter->getInputDims()[1];
  int width = undistorter->getOutputDims()[0];
  int height = undistorter->getOutputDims()[1];

  float* internalBuffer = new float[widthOrg*heightOrg];

  // Initialize VideoWriter
  /*
  cv::VideoWriter writer;
  int codec = cv::VideoWriter::fourcc('a', 'v', 'c', '1');
  cv::Size frameSize(width, height);
  writer.open("./whitePaper_rectified.mp4", codec, 10.0, frameSize, true);
  */

  Timer timer;
  for (const std::string &file : getSequence<unsigned>("./datasets/tum/calib_wide_whitePaper/images/")) {
    cv::Mat imageRaw;
    imageRaw = cv::imread(file, cv::IMREAD_GRAYSCALE);

    timer.Start();
    undistorter->undistort<unsigned char>(imageRaw.data, internalBuffer, widthOrg*heightOrg, width*height);
    rectification.push_back(timer.ElapsedMicroSeconds());

    cv::Mat imageRectified(height, width, CV_32F, internalBuffer);

    std::vector<lineag> lines;

    timer.Start();
    ag3line(imageRectified, lines, true);
    line_extraction.push_back(timer.ElapsedMicroSeconds());
    
    /*
    cv::Mat image;
    imageRectified.convertTo(image, CV_8U);

    cv::Mat imageLines;
    cv::cvtColor(image, imageLines, cv::COLOR_GRAY2RGB);
    for (const lineag &line : lines)
      cv::line(imageLines, cv::Point(line.x1, line.y1), cv::Point(line.x2, line.y2), cv::Scalar(0, 200, 0), 1);
    writer.write(imageLines);
    */
  }

  LOG(INFO) << StringPrintf("Rectification: %.1f ms", Mean(rectification)*1e-3);
  LOG(INFO) << StringPrintf("AG3Line: %.1f ms", Mean(line_extraction)*1e-3);

  delete[] internalBuffer;

  return 0;
}
