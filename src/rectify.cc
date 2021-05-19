
#define PROGRAM_NAME                                                                               \
  "rectify"

#define FLAGS_CASES                                                                                \
  FLAG_CASE(string, output_file, "./rectified.png", "Rectified output image [PNG]")

#define ARGS_CASES                                                                                 \
  ARG_CASE(image)

#include <memory>
#include <iostream>
#include <vector>

// boost
#include <boost/filesystem.hpp>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>

#include "FOVUndistorter.h"

#include "util/args.hpp"
#include "util/logging.h"

namespace fs = boost::filesystem;

void ValidateArgs() { }

void ValidateFlags() {
  CHECK_EQ(fs::path(FLAGS_output_file).extension().string(), ".png");
}

int main(int argc, char *argv[]) {
  
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
  
  // image sizes
  const int widthOrg = undistorter->getInputDims()[0];
  const int heightOrg = undistorter->getInputDims()[1];
  const int width = undistorter->getOutputDims()[0];
  const int height = undistorter->getOutputDims()[1];
  
  float* internalBuffer = new float[widthOrg*heightOrg];
  
  cv::Mat imageRaw;
  imageRaw = cv::imread(ARGS_image, cv::IMREAD_GRAYSCALE);
  
  undistorter->undistort<unsigned char>(imageRaw.data, internalBuffer, widthOrg*heightOrg, width*height);
  
  cv::Mat imageRectified(height, width, CV_32F, internalBuffer);

  cv::Mat image;
  imageRectified.convertTo(image, CV_8U);
  
  // cv::IMWRITE_PNG_STRATEGY, [cv::IMWRITE_PNG_STRATEGY_DEFAULT, cv::IMWRITE_PNG_STRATEGY_FILTERED, cv::cv::IMWRITE_PNG_STRATEGY_RLE]
  std::vector<int> compression_params = {cv::IMWRITE_PNG_COMPRESSION, 9};
  
  cv::imwrite(FLAGS_output_file, image, compression_params);
  
  delete[] internalBuffer;
  
  return 0;
}
