
#include <memory>
#include <vector>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>

#include "FOVUndistorter.h"

#include "thirdparty/LSD/include/lsd.h"

#include "sequence.hpp"

#include "util/logging.h"
#include "util/math.h"
#include "util/random.h"
#include "util/string.h"
#include "util/timer.h"

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
  cv::VideoWriter writer;
  int codec = cv::VideoWriter::fourcc('a', 'v', 'c', '1');
  cv::Size frameSize(width, height);
  writer.open("./whitePaper_rectified.mp4", codec, 10.0, frameSize, true);

  cv::Ptr<LineSegmentDetector> lsd = createLineSegmentDetector(LSD_REFINE_NONE);

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
    std::vector<cv::Vec4f> lines;
    lsd->detect(image, lines);
    line_extraction.push_back(timer.ElapsedMicroSeconds());

    cv::Mat imageLines;
    cv::cvtColor(image, imageLines, cv::COLOR_GRAY2RGB);
    for (const cv::Vec4f &line : lines)
      cv::line(imageLines, cv::Point(line[0], line[1]), cv::Point(line[2], line[3]), cv::Scalar(0, 200, 0), 2);
    writer.write(imageLines);
  }

  LOG(INFO) << StringPrintf("Rectification: %.1f ms", Mean(rectification)*1e-3);
  LOG(INFO) << StringPrintf("LSD: %.1f ms", Mean(line_extraction)*1e-3);

  delete[] internalBuffer;

  return 0;
}
