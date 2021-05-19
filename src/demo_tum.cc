
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// gflags
#include <gflags/gflags.h>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "camera_models.h"
#include "FOVUndistorter.h"

#include "EDLines.h"

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

    cv::Mat image;
    imageRectified.convertTo(image, CV_8U);

    timer.Start();
    EDLines edlines = EDLines(image, 1.0, -1, 6.0, 1.3);
    line_extraction.push_back(timer.ElapsedMicroSeconds());

    std::vector<std::vector<cv::Point>> segments = edlines.getSegments();
    std::vector<LineSegment> lines = edlines.getLines();

    /*
    cv::Mat imageLines;
    cv::cvtColor(image, imageLines, cv::COLOR_GRAY2RGB);
    for (const LineSegment &line : lines) {
      if (line.len < MIN_SINGLE_LINE_LENGTH) continue;
      const std::vector<cv::Point> &segment = segments[line.segmentNo];
      for (int k = 0; k < line.len; k++)
        cv::circle(imageLines, segment[line.firstPixelIndex + k], 1, cv::Scalar(0, 200, 0), -1);
    }
    writer.write(imageLines);
    */
  }

  LOG(INFO) << StringPrintf("Rectification: %.1f ms", Mean(rectification)*1e-3);
  LOG(INFO) << StringPrintf("EDLines: %.1f ms", Mean(line_extraction)*1e-3);

  delete[] internalBuffer;

  return 0;
}
