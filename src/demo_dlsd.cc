
#include "util/alignment.h"
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(Eigen::Vector3d)

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <utility>

// eigen
#include <Eigen/Core>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>

#include "camera_models.h"

#include "ED.h"
#include "EDLines.h"

#include "line_fit.h"

#include "sequence.hpp"

#include "util/logging.h"
#include "util/math.h"
#include "util/random.h"
#include "util/string.h"
#include "util/timer.h"

// Parameters
//  Line Segments
#define MIN_LINE_LENGTH 9 // px
#define MAX_LINE_ERROR 1.0 // px
//  Polygonal Chains
#define SMOOTHNESS_TH 0.19634954084936207 // rad (M_PI/16)
#define REPROJECTION_TH 4.0 // px
#define MIN_SINGLE_LINE_LENGTH 100 // px

struct GapSegment {
  // TODO segmentNo can be ommited
  int segmentNo;       // Edge segment that this gap belongs to
  int firstPixelIndex; // Index of the first pixel of this gap
  int len;             // No of pixels making up the gap segment

  GapSegment()
    : segmentNo(0), firstPixelIndex(0), len(0) { }

  GapSegment(int _segmentNo, int _firstPixelIndex, int _len)
    : segmentNo(_segmentNo), firstPixelIndex(_firstPixelIndex), len(_len){ }
};

struct PolygonalChain {
  int turn;           // Turn direction: 1 or -1
  int segmentNo;      // Edge segment that this polygonal path belongs to

  std::vector<int> lines;
  std::vector<GapSegment> gaps;
};

struct DistortedLineSegment {
  int segmentNo;       // Edge segment that this polygonal path belongs to

  Eigen::Vector3d axis;
  std::vector<cv::Point> segment;
  std::unordered_set<int> invalidPixels;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

//-----------------------------------------------------------------------------------------
// Computes the minimum line length using the NFA formula given width & height values
int ComputeMinLineLength(int width, int height) {
	// The reason we are dividing the theoretical minimum line length by 2 is because
	// we now test short line segments by a line support region rectangle having width=2.
	// This means that within a line support region rectangle for a line segment of length "l"
	// there are "2*l" many pixels. Thus, a line segment of length "l" has a chance of getting
	// validated by NFA.

	const double logNT = 2.0*(std::log10((double)width) + std::log10((double)height));
	return (int) std::round((-logNT / std::log10(0.125))*0.5);
} //end-ComputeMinLineLength

//-----------------------------------------------------------------------------------
// Fits a line of the form y=a+bx (invert == 0) OR x=a+by (invert == 1)
//
void LineFit(double * x, double * y, int count, double &a, double &b, double &e, int &invert) {
	if (count<2) return;

	double S = count, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
	for (int i = 0; i<count; i++) {
		Sx += x[i];
		Sy += y[i];
	} //end-for

	double mx = Sx / count;
	double my = Sy / count;

	double dx = 0.0;
	double dy = 0.0;
	for (int i = 0; i < count; i++) {
		dx += (x[i] - mx)*(x[i] - mx);
		dy += (y[i] - my)*(y[i] - my);
	} //end-for

	if (dx < dy) {
		// Vertical line. Swap x & y, Sx & Sy
		invert = 1;
		double *t = x;
		x = y;
		y = t;

		double d = Sx;
		Sx = Sy;
		Sy = d;

	}
	else {
		invert = 0;
	} //end-else

	  // Now compute Sxx & Sxy
	for (int i = 0; i<count; i++) {
		Sxx += x[i] * x[i];
		Sxy += x[i] * y[i];
	} //end-for

	double D = S*Sxx - Sx*Sx;
	a = (Sxx*Sy - Sx*Sxy) / D;
	b = (S  *Sxy - Sx* Sy) / D;

	if (b == 0.0) {
		// Vertical or horizontal line
		double error = 0.0;
		for (int i = 0; i<count; i++) {
			error += fabs((a) - y[i]);
		} //end-for
		e = error / count;

	}
	else {
		double error = 0.0;
		for (int i = 0; i<count; i++) {
			// Let the line passing through (x[i], y[i]) that is perpendicular to a+bx be c+dx
			double d = -1.0 / (b);
			double c = y[i] - d*x[i];
			double x2 = ((a) - c) / (d - (b));
			double y2 = (a) + (b)*x2;

			double dist = (x[i] - x2)*(x[i] - x2) + (y[i] - y2)*(y[i] - y2);
			error += dist;
		} //end-for

		e = std::sqrt(error / count);
	} //end-else
}

//-----------------------------------------------------------------------------------
// Fits a line of the form y=a+bx (invert == 0) OR x=a+by (invert == 1)
// Assumes that the direction of the line is known by a previous computation
//
void LineFit(double * x, double * y, int count, double &a, double &b, int invert) {
	if (count<2) return;

	double S = count, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
	for (int i = 0; i<count; i++) {
		Sx += x[i];
		Sy += y[i];
	} //end-for

	if (invert) {
		// Vertical line. Swap x & y, Sx & Sy
		double *t = x;
		x = y;
		y = t;

		double d = Sx;
		Sx = Sy;
		Sy = d;
	} //end-if

	  // Now compute Sxx & Sxy
	for (int i = 0; i<count; i++) {
		Sxx += x[i] * x[i];
		Sxy += x[i] * y[i];
	} //end-for

	double D = S*Sxx - Sx*Sx;
	a = (Sxx*Sy - Sx*Sxy) / D;
	b = (S  *Sxy - Sx* Sy) / D;
}

double ComputeMinDistance(double x1, double y1, double a, double b, int invert) {
	double x2, y2;

	if (invert == 0) {
		if (b == 0) {
			x2 = x1;
			y2 = a;

		}
		else {
			// Let the line passing through (x1, y1) that is perpendicular to a+bx be c+dx
			double d = -1.0 / (b);
			double c = y1 - d*x1;

			x2 = (a - c) / (d - b);
			y2 = a + b*x2;
		} //end-else

	}
	else {
		/// invert = 1
		if (b == 0) {
			x2 = a;
			y2 = y1;

		}
		else {
			// Let the line passing through (x1, y1) that is perpendicular to a+by be c+dy
			double d = -1.0 / (b);
			double c = x1 - d*y1;

			y2 = (a - c) / (d - b);
			x2 = a + b*y2;
		} //end-else
	} //end-else

	return std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

//---------------------------------------------------------------------------------
// Given a point (x1, y1) and a line equation y=a+bx (invert=0) OR x=a+by (invert=1)
// Computes the (x2, y2) on the line that is closest to (x1, y1)
//
void ComputeClosestPoint(double x1, double y1, double a, double b, int invert, double &xOut, double &yOut) {
	double x2, y2;

	if (invert == 0) {
		if (b == 0) {
			x2 = x1;
			y2 = a;

		}
		else {
			// Let the line passing through (x1, y1) that is perpendicular to a+bx be c+dx
			double d = -1.0 / (b);
			double c = y1 - d*x1;

			x2 = (a - c) / (d - b);
			y2 = a + b*x2;
		} //end-else

	}
	else {
		/// invert = 1
		if (b == 0) {
			x2 = a;
			y2 = y1;

		}
		else {
			// Let the line passing through (x1, y1) that is perpendicular to a+by be c+dy
			double d = -1.0 / (b);
			double c = x1 - d*y1;

			y2 = (a - c) / (d - b);
			x2 = a + b*y2;
		} //end-else
	} //end-else

	xOut = x2;
	yOut = y2;
}

//-----------------------------------------------------------------
// Given a full segment of pixels, splits the chain to lines
// This code is used when we use the whole segment of pixels
//
void SplitSegment2Lines(double * x, double * y, int noPixels, int segmentNo, std::vector<LineSegment> &lines, int min_line_len = MIN_LINE_LENGTH, double line_error = MAX_LINE_ERROR) {

  lines.clear();

	// First pixel of the line segment within the segment of points
	int firstPixelIndex = 0;

	while (noPixels >= min_line_len) {
		// Start by fitting a line to MIN_LINE_LEN pixels
		bool valid = false;
		double lastA, lastB, error;
		int lastInvert;

		while (noPixels >= min_line_len) {
			LineFit(x, y, min_line_len, lastA, lastB, error, lastInvert);
			if (error <= 0.5) { valid = true; break; }

#if 1
			noPixels -= 1;   // Go slowly
			x += 1; y += 1;
			firstPixelIndex += 1;
#else
			noPixels -= 2;   // Go faster (for speed)
			x += 2; y += 2;
			firstPixelIndex += 2;
#endif
		} //end-while

		if (valid == false) return;

		// Now try to extend this line
		int index = min_line_len;
		int len = min_line_len;

		while (index < noPixels) {
			int startIndex = index;
			int lastGoodIndex = index - 1;
			int goodPixelCount = 0;
			int badPixelCount = 0;
			while (index < noPixels) {
				double d = ComputeMinDistance(x[index], y[index], lastA, lastB, lastInvert);

				if (d <= line_error) {
					lastGoodIndex = index;
					goodPixelCount++;
					badPixelCount = 0;

				}
				else {
					badPixelCount++;
					if (badPixelCount >= 5) break;
				} //end-if

				index++;
			} //end-while

			if (goodPixelCount >= 2) {
				len += lastGoodIndex - startIndex + 1;
				LineFit(x, y, len, lastA, lastB, lastInvert);  // faster LineFit
				index = lastGoodIndex + 1;
			} // end-if

			if (goodPixelCount < 2 || index >= noPixels) {
				// End of a line segment. Compute the end points
				double sx, sy, ex, ey;

				int index = 0;
				while (ComputeMinDistance(x[index], y[index], lastA, lastB, lastInvert) > line_error) index++;
				ComputeClosestPoint(x[index], y[index], lastA, lastB, lastInvert, sx, sy);
				int noSkippedPixels = index;

				index = lastGoodIndex;
				while (ComputeMinDistance(x[index], y[index], lastA, lastB, lastInvert) > line_error) index--;
				ComputeClosestPoint(x[index], y[index], lastA, lastB, lastInvert, ex, ey);

				// Add the line segment to lines
				lines.push_back(LineSegment(lastA, lastB, lastInvert, sx, sy, ex, ey, segmentNo, firstPixelIndex + noSkippedPixels, index - noSkippedPixels + 1));
				//linesNo++;
				len = index + 1;
				break;
			} //end-else
		} //end-while

		noPixels -= len;
		x += len;
		y += len;
		firstPixelIndex += len;
	} //end-while
}

// TODO description
void GroupLines2PolygonalChains(const double *x, const double *y, int noPixels, int segmentNo, const std::vector<LineSegment> &lines, std::vector<PolygonalChain> &chains) {

  chains.clear();

  int noLines = lines.size();
  if (noLines == 0) return; // no lines, no polygons...

  int firstPixelIndex = 0;

  PolygonalChain chain; // current chain
  chain.turn = 0; // 0: not set; 1: counter clockwise; -1: clockwise;
  chain.segmentNo = segmentNo;
  chain.lines.push_back(0);
  chain.gaps.emplace_back(segmentNo, firstPixelIndex, lines[0].firstPixelIndex - firstPixelIndex); // pregap

  if (noLines == 1) {
    firstPixelIndex = lines[0].firstPixelIndex + lines[0].len;
    chain.gaps.emplace_back(segmentNo, firstPixelIndex, noPixels - firstPixelIndex); // postgap
    return;
  }

  // start line grouping
  GapSegment lastGap;
  for (int index = 1; index < noLines; index++) {
    const LineSegment *l1 = &lines[index-1];
    const LineSegment *l2 = &lines[index];

    // compute the angle between the lines & their turn direction
		double v1x = l1->ex - l1->sx;
		double v1y = l1->ey - l1->sy;
		double v1Len = std::sqrt(v1x*v1x + v1y*v1y);

		double v2x = l2->ex - l2->sx;
		double v2y = l2->ey - l2->sy;
		double v2Len = std::sqrt(v2x*v2x + v2y*v2y);

		// fill gap info between l1 and l2
		firstPixelIndex = l1->firstPixelIndex + l1->len;
		lastGap.segmentNo = segmentNo;
		lastGap.firstPixelIndex = firstPixelIndex;
		lastGap.len = l2->firstPixelIndex - firstPixelIndex;
		chain.gaps.push_back(lastGap); // copy

		double dotProduct = (v1x*v2x + v1y*v2y) / (v1Len*v2Len);
		if (dotProduct > 1.0) dotProduct = 1.0;
		else if (dotProduct < -1.0) dotProduct = -1.0;

		double angle = std::acos(dotProduct);
		int sign = (v1x*v2y - v2x*v1y) >= 0 ? 1 : -1;  // cross product

		//LOG(INFO) << "Angle: " << RadToDeg(angle);
		//LOG(INFO) << "Sign: " << sign;

		if (angle <= SMOOTHNESS_TH && chain.turn == sign) {
      //LOG(INFO) << "append line to current polygonal chain";

      // continue polyginal chain
      chain.lines.push_back(index);
		} else if (angle <= SMOOTHNESS_TH && chain.turn == 0) {
      //LOG(INFO) << "start new polygonal chain";

      // new polygonal chain
      chain.turn = sign;
      chain.lines.push_back(index);
		} else { // angle > SMOOTHNESS_TH || (chain.turn != sign && chain.turn != 0)
		  //LOG(INFO) << "end current polygonal chain";

		  // save current polygonal chain
		  chains.push_back(std::move(chain));

		  // re-initialize current chain
		  chain.turn = 0;
		  chain.segmentNo = segmentNo;
		  chain.lines.clear();
		  chain.lines.push_back(index);
		  chain.gaps.clear();
		  chain.gaps.push_back(std::move(lastGap));
		}
  }

  // postgap
  firstPixelIndex = lines[noLines-1].firstPixelIndex + lines[noLines-1].len;
	lastGap.segmentNo = segmentNo;
	lastGap.firstPixelIndex = firstPixelIndex;
	lastGap.len = noPixels - firstPixelIndex;
	chain.gaps.push_back(std::move(lastGap));

  // check closed contour
  double dx = x[0] - x[noPixels-1];
	double dy = y[0] - y[noPixels-1];
	double d = std::sqrt(dx*dx + dy*dy);
	//LOG(INFO) << "Distance for closed contour: " << d;

	if (d < 1.5) { // 1.5 ~ sqrt(2)
	  // closed
	  const LineSegment *l1 = &lines[noLines-1];
	  const LineSegment *l2 = &lines[0];

    // compute the angle between the lines & their turn direction
	  double v1x = l1->ex - l1->sx;
	  double v1y = l1->ey - l1->sy;
	  double v1Len = std::sqrt(v1x*v1x + v1y*v1y);

	  double v2x = l2->ex - l2->sx;
	  double v2y = l2->ey - l2->sy;
	  double v2Len = std::sqrt(v2x*v2x + v2y*v2y);

	  double dotProduct = (v1x*v2x + v1y*v2y) / (v1Len*v2Len);
	  if (dotProduct > 1.0) dotProduct = 1.0;
	  else if (dotProduct < -1.0) dotProduct = -1.0;

	  double angle = std::acos(dotProduct);
	  int sign = (v1x*v2y - v2x*v1y) >= 0 ? 1 : -1;  // cross product

	  //LOG(INFO) << "Angle: " << RadToDeg(angle);
	  //LOG(INFO) << "Sign: " << sign;

	  if (angle <= SMOOTHNESS_TH && (chain.turn == sign || chain.turn == 0)) {
	    //LOG(INFO) << "append line to first polygonal chain";

      // extend first polyginal chain
      chain.turn = sign;
      chains[0].turn = sign;

      chain.lines.insert(chain.lines.end(), chains[0].lines.begin(), chains[0].lines.end());
      chain.gaps.insert(chain.gaps.end(), chains[0].gaps.begin(), chains[0].gaps.end());
      chains[0].lines = std::move(chain.lines);
      chains[0].gaps = std::move(chain.gaps);

      chain.lines.clear();
	  }
	}

	if (!chain.lines.empty()) // unsaved chain
    chains.push_back(std::move(chain)); // save current polygonal chain
}

void ComputeClosestPoint(double x, double y, const double params[5], const Eigen::Vector3d &axis, double &xOut, double &yOut) {
  Eigen::Vector3d p;
  FOVCameraModel::ImageToWorld<double>(params, x, y, &p.x(), &p.y());
  p.z() = 1.;

  Eigen::Vector3d v = p - p.dot(axis)*axis; // direction vector projected onto fitted plane
  v /= v.z(); // normalized image plane projection

  FOVCameraModel::WorldToImage<double>(params, v.x(), v.y(), &xOut, &yOut);
}

// TODO Description
void ValidatePolygonalChain(const double *x, const double *y, int segmentNo, const double params[5], const std::vector<PolygonalChain> &chains, const std::vector<LineSegment> &lines, std::vector<DistortedLineSegment> &distortedLineSegments, int min_single_line_len = MIN_SINGLE_LINE_LENGTH) {
  // Validation
  // TODO Single lines (i.e. not assigned to a cluster) can be accepted without validation
  // Line fitting (3D) with line segment end points:
  //  1. Undistort end points
  //  2. Perform line fitting (3D)
  //  3. Project 3D end points to the fitted plane
  //  4. Project "model" end points to image and accept if all end points are within a reprojection error

  for (const PolygonalChain &chain : chains) {
    int noLines = chain.lines.size();

    // filter out small single line polygonal chains
    if (noLines == 1 && lines[chain.lines[0]].len < min_single_line_len) continue;

    std::vector<Eigen::Vector3d> endpoints(2*noLines);
    std::vector<Eigen::Vector3d> allLinePixels;
    for (int index = 0; index < noLines; index++) {
      const LineSegment &line = lines[chain.lines[index]];

      Eigen::Vector3d s; // start point
      FOVCameraModel::ImageToWorld<double>(params, line.sx, line.sy, &s.x(), &s.y());
      s.z() = 1.;
      endpoints[2*index] = s.normalized();

      Eigen::Vector3d e; // end point
      FOVCameraModel::ImageToWorld<double>(params, line.ex, line.ey, &e.x(), &e.y());
      e.z() = 1.;
      endpoints[2*index+1] = e.normalized();

      if (noLines > 1) {
        for (int k = 0; k < line.len; k++) {
          Eigen::Vector3d v;
          FOVCameraModel::ImageToWorld<double>(params, x[line.firstPixelIndex + k], y[line.firstPixelIndex + k], &v.x(), &v.y());
          v.z() = 1.;
          allLinePixels.push_back(v.normalized());
        }
      }
    }

    Eigen::Vector3d axis;
    if (noLines > 1) {
      bool success = line_fit(allLinePixels, axis);
      if (!success) {
        LOG(INFO) << "Solver failed!";
        continue;
      }
    } else // Single line polygonal chains do not need validation
      axis = endpoints[0].cross(endpoints[1]).normalized();
    //LOG(INFO) << "Axis: " << axis.transpose();

    bool valid = true;
    if (noLines > 1) {
      for (int index = 0; index < noLines; index++) {
        double x, y;
        double dx, dy, d;
        const LineSegment &line = lines[chain.lines[index]];

        //cv::circle(imageValidation, cv::Point2d(line.sx, line.sy), 2, cv::Scalar(0, 200, 0), -1);

        Eigen::Vector3d s = endpoints[2*index] - endpoints[2*index].dot(axis)*axis; // start point projected onto fitted plane
        s /= s.z(); // normalized image plane projection
        FOVCameraModel::WorldToImage<double>(params, s.x(), s.y(), &x, &y);
        dx = line.sx - x;
        dy = line.sy - y;
        d = std::sqrt(dx*dx + dy*dy);

        if (d > REPROJECTION_TH) {
          //LOG(WARNING) << "Reprojection error: " << d;
          valid = false;
          break;
        }

        Eigen::Vector3d e = endpoints[2*index+1] - endpoints[2*index+1].dot(axis)*axis; // end point projected onto fitted plane
        e /= e.z(); // normalized image plane projection
        FOVCameraModel::WorldToImage<double>(params, e.x(), e.y(), &x, &y);
        dx = line.ex - x;
        dy = line.ey - y;
        d = std::sqrt(dx*dx + dy*dy);

        if (d > REPROJECTION_TH) {
          //LOG(WARNING) << "Reprojection error: " << d;
          valid = false;
          break;
        }
      }
    }

    if (!valid) continue;
    //LOG(INFO) << "Polygonal chain validated!";

    DistortedLineSegment distortedLineSegment;
    distortedLineSegment.segmentNo = segmentNo;
    distortedLineSegment.axis = axis;

    // segment
    for (int line_index = 0; line_index <= noLines; line_index++) {
      // validate pregaps
      const GapSegment &gap = chain.gaps[line_index];
      bool validPixelInGap = false;
      for (int k = 0; k < gap.len; k++) {
        const int index = gap.firstPixelIndex + k;

        double _x, _y;
        ComputeClosestPoint(x[index], y[index], params, axis, _x, _y);

        double dx = x[index] - _x;
        double dy = y[index] - _y;
        double d = std::sqrt(dx*dx + dy*dy);

        if (line_index == 0) { // pregap
          if (d > REPROJECTION_TH && validPixelInGap)
            distortedLineSegment.invalidPixels.insert(index);
          else if (d <= REPROJECTION_TH) {
            validPixelInGap = true;
            distortedLineSegment.segment.emplace_back(cv::Point(x[index], y[index]));
          }
        } else if (line_index == noLines) { // postgap
          if (d > REPROJECTION_TH) break;
          else distortedLineSegment.segment.emplace_back(cv::Point(x[index], y[index]));
        } else { // intra chain gaps
          if (d > REPROJECTION_TH) distortedLineSegment.invalidPixels.insert(index);
          else distortedLineSegment.segment.emplace_back(cv::Point(x[index], y[index]));
        }
      }

      if (line_index < noLines) {
        // lines validated by default, so just add them to the distorted line segment
        const LineSegment &line = lines[chain.lines[line_index]];
        for (int k = 0; k < line.len; k++)
          distortedLineSegment.segment.emplace_back(cv::Point(x[line.firstPixelIndex + k], y[line.firstPixelIndex + k]));
      }
    }

    distortedLineSegments.push_back(std::move(distortedLineSegment));
  }
}

int main(int argc, char *argv[]) {
  (void) argc; // unused

  std::vector<double> edge_drawing;
  std::vector<double> line_extraction;
  std::vector<double> line_grouping;
  std::vector<double> validation;

  // Initialize Google's logging library.
  InitializeGlog(argv);
  FLAGS_logtostderr = 1;

  // https://stackoverflow.com/a/13062069
  std::shared_ptr<double[]> params(new double[5]);
  params[0] = 446.915863;
  params[1] = 447.071228;
  params[2] = 631.219238;
  params[3] = 510.997498;
  params[4] = 0.933271;

  int width = 1280;
  int height = 1024;

  // Temporary buffers used during line fitting
  double *x = new double[(width + height) * 8];
  double *y = new double[(width + height) * 8];

  // Initialize VideoWriter
  /*
  cv::VideoWriter writer;
  int codec = cv::VideoWriter::fourcc('a', 'v', 'c', '1');
  cv::Size frameSize(width, height);
  writer.open("./whitePaper.mp4", codec, 10.0, frameSize, true);
  */

  Timer timer;
  for (const std::string &file : getSequence<unsigned>("./datasets/tum/calib_wide_whitePaper/images/")) {
    cv::Mat imageRaw;
    imageRaw = cv::imread(file, cv::IMREAD_GRAYSCALE);

    timer.Start();
    //Call ED constructor
    ED ed = ED(imageRaw, SOBEL_OPERATOR, 36, 8, 1, 10, 1.0, true); // apply ED algorithm
    edge_drawing.push_back(timer.ElapsedMicroSeconds());

    //Output number of segments
    //int noSegments = ed.getSegmentNo();
    //LOG(INFO) << "Number of edge segments: " << noSegments;

    double total_line_extraction = 0.0;
    double total_line_grouping = 0.0;
    double total_validation = 0.0;

    int segmentNo = 0;
    std::vector<DistortedLineSegment> distortedLineSegments;
    for (const std::vector<cv::Point> &segment : ed.getSortedSegments()) {
      int noPixels = segment.size();

      // Lines from segment
      //double line_error = 1.0;
      //int min_line_len = std::max(ComputeMinLineLength(width, height), 10);
      //LOG(INFO) << "Min line length: " << min_line_len;

	    CHECK_LT(segment.size(), (width + height) * 8);
	    //double *x = new double[noPixels];
	    //double *y = new double[noPixels];
	    for (int k = 0; k < noPixels; k++) {
		    x[k] = segment[k].x;
		    y[k] = segment[k].y;
	    }

      timer.Start();
	    // SplitSegment2Lines
	    std::vector<LineSegment> lines;
	    SplitSegment2Lines(x, y, noPixels, segmentNo, lines, MIN_LINE_LENGTH, MAX_LINE_ERROR);
      total_line_extraction += timer.ElapsedMicroSeconds();

	    //int noLines = lines.size();
	    //LOG(INFO) << "Number of lines: " << noLines;

	    //cv::Mat imageColor;
      //cv::cvtColor(imageRaw, imageColor, cv::COLOR_GRAY2RGB);
	    //for (const LineSegment &line : lines) {
	    //  cv::Scalar color(RandomInteger(0, 255), RandomInteger(0, 255), RandomInteger(0, 255));
	    //  for (int k = 0; k < line.len; k++)
	    //    cv::circle(imageColor, cv::Point2d(x[line.firstPixelIndex + k], y[line.firstPixelIndex + k]), 2, color, -1);
	    //}
	    //cv::imshow("Lines", imageColor);
      //cv::waitKey();

      timer.Start();
      std::vector<PolygonalChain> chains;
      GroupLines2PolygonalChains(x, y, noPixels, segmentNo, lines, chains);
      total_line_grouping += timer.ElapsedMicroSeconds();

      //int noChains = chains.size();
      //LOG(INFO) << "Number of polygonal chains: " << noChains;

      //cv::Mat imageChains;
      //cv::cvtColor(imageRaw, imageChains, cv::COLOR_GRAY2RGB);
	    //for (const PolygonalChain &chain : chains) {
	    //  cv::Scalar color(RandomInteger(0, 255), RandomInteger(0, 255), RandomInteger(0, 255));
	    //  for (int index : chain.lines) {
	    //    const LineSegment &line = lines[index];
	    //    for (int k = 0; k < line.len; k++)
	    //      cv::circle(imageChains, cv::Point2d(x[line.firstPixelIndex + k], y[line.firstPixelIndex + k]), 2, color, -1);
	    //    for (const GapSegment &gap : chain.gaps) {
	    //      for (int k = 0; k < gap.len; k++)
	    //        cv::circle(imageChains, cv::Point2d(x[gap.firstPixelIndex + k], y[gap.firstPixelIndex + k]), 2, cv::Scalar(0, 200, 0), -1);
	    //    }
	    //  }
	    //}
	    //cv::imshow("Polygonal chains", imageChains);
      //cv::waitKey();

      timer.Start();
      ValidatePolygonalChain(x, y, segmentNo, params.get(), chains, lines, distortedLineSegments);
	    total_validation += timer.ElapsedMicroSeconds();

	    segmentNo++;
	  }

	  line_extraction.push_back(total_line_extraction);
	  line_grouping.push_back(total_line_grouping);
	  validation.push_back(total_validation);

    /*
    cv::Mat imageDistortedLines;
    cv::cvtColor(imageRaw, imageDistortedLines, cv::COLOR_GRAY2RGB);
    for (const DistortedLineSegment &distortedLine : distortedLineSegments) {
      for (const cv::Point point : distortedLine.segment)
        cv::circle(imageDistortedLines, point, 1, cv::Scalar(0, 200, 0), -1);
    }
    writer.write(imageDistortedLines);
    */
  }

  LOG(INFO) << StringPrintf("Edge Drawing: %.1f ms", Mean(edge_drawing)*1e-3);
  LOG(INFO) << StringPrintf("Line Extraction: %.1f ms", Mean(line_extraction)*1e-3);
  LOG(INFO) << StringPrintf("Line Grouping: %.1f ms", Mean(line_grouping)*1e-3);
  LOG(INFO) << StringPrintf("Line Validation: %.1f ms", Mean(validation)*1e-3);

  //Show resulting edge image
  //cv::Mat edgeImg = testED.getEdgeImage();
  //cv::imshow("Edges", edgeImg);
  //cv::waitKey();

/*
  EDLines edlines = EDLines(ed);

  //Acquiring line information, i.e. start & end points
  std::vector<LS> lines = edlines.getLines();
  int noLines = edlines.getLinesNo();
  std::cout << "Number of line segments: " << noLines << std::endl;

  cv::Mat lineImg = edlines.getLineImage(); //draws on an empty image
  cv::imshow("Line Image", lineImg);
  cv::waitKey();
*/

  delete[] x;
  delete[] y;

  return 0;
}

  // TODO Line Segments discontinuity:
  //        1. Automatically validate line segments discontinuity less that or equal to 3 px (assume line connecting end points)
  //        2. Assume line (connect discontinuity end points) and validate it (the discontinuity if usually caused by the min_lin_len and 1 px fit threshold)
  //        3. Validate all segmented pixels are within 1 px (thus ignoring min_lin_len)
  //      [SplitSegment2Lines should produce line segments discontinuity paths for the above]
/*
  std::vector<std::vector<int>> lineClusters;
  std::vector<int> clusterSign;
  std::vector<int> currentCluster;
  std::vector<int> info(noLines, -1);
  int currentDirection = 0;
  for (int k = 1; k < noLines; k++) {
    LineSegment *l1 = &lines[k-1];
    LineSegment *l2 = &lines[k];

    // Compute the angle between the lines & their turn direction
		double v1x = l1->ex - l1->sx;
		double v1y = l1->ey - l1->sy;
		double v1Len = std::sqrt(v1x*v1x + v1y*v1y);

		double v2x = l2->ex - l2->sx;
		double v2y = l2->ey - l2->sy;
		double v2Len = std::sqrt(v2x*v2x + v2y*v2y);

		double dotProduct = (v1x*v2x + v1y*v2y) / (v1Len*v2Len);
		if (dotProduct > 1.0) dotProduct = 1.0;
		else if (dotProduct < -1.0) dotProduct = -1.0;

		double angle = std::acos(dotProduct);
		int sign = (v1x*v2y - v2x*v1y) >= 0 ? 1 : -1;  // compute cross product
		//info[j].taken = false;

		LOG(INFO) << "Angle: " << RadToDeg(angle);
		LOG(INFO) << "Sign: " << sign;

		if (angle < SMOOTHNESS_TH && sign == currentDirection) {
		  LOG(INFO) << "append line to current cluster";
		  currentCluster.push_back(k); // append to current cluster
		  info[k] = lineClusters.size(); // assign cluster index
		} else if (currentDirection == 0 && angle < SMOOTHNESS_TH) {
	    // New cluster
	    LOG(INFO) << "start new cluster";
	    currentDirection = sign;
	    currentCluster.push_back(k-1);
	    currentCluster.push_back(k);
	    // assign cluster index
	    info[k-1] = lineClusters.size();
	    info[k] = lineClusters.size();
		} else if (!currentCluster.empty()) {
		  LOG(INFO) << "end current cluster";
		  // Save current cluster
		  lineClusters.push_back(currentCluster);
		  clusterSign.push_back(currentDirection);
		  currentCluster.clear();
		  currentDirection = 0;
		}
  }
  // Check closed loop
  double dx = lines[0].sx - lines[noLines - 1].ex;
	double dy = lines[0].sy - lines[noLines - 1].ey;
	double d = std::sqrt(dx*dx + dy*dy);
	LOG(INFO) << "Distance: " << d;
	if (d <= 10) { // min_lin_len + 1
	  // Close enough
	  LineSegment *l1 = &lines[noLines-1];
	  LineSegment *l2 = &lines[0];

    // Compute the angle between the lines & their turn direction
	  double v1x = l1->ex - l1->sx;
	  double v1y = l1->ey - l1->sy;
	  double v1Len = std::sqrt(v1x*v1x + v1y*v1y);

	  double v2x = l2->ex - l2->sx;
	  double v2y = l2->ey - l2->sy;
	  double v2Len = std::sqrt(v2x*v2x + v2y*v2y);

	  double dotProduct = (v1x*v2x + v1y*v2y) / (v1Len*v2Len);
	  if (dotProduct > 1.0) dotProduct = 1.0;
	  else if (dotProduct < -1.0) dotProduct = -1.0;

	  double angle = std::acos(dotProduct);
	  int sign = (v1x*v2y - v2x*v1y) >= 0 ? 1 : -1;  // compute cross product
	  //info[j].taken = false;

	  LOG(INFO) << "Angle: " << RadToDeg(angle);
	  LOG(INFO) << "Sign: " << sign;

	  if (angle < SMOOTHNESS_TH && sign == currentDirection) {
	    LOG(INFO) << "append line to current cluster";
	    if (info[0] == -1) {
	      currentCluster.push_back(0); // append to current cluster
	      info[0] = lineClusters.size(); // assign cluster index
	    } else {
	      // Join current cluster with cluster l2's cluster
	      currentCluster.insert(currentCluster.end(), lineClusters[info[0]].begin(), lineClusters[info[0]].end());
	      lineClusters[info[0]] = currentCluster;
	      currentCluster.clear();
	    }
	  } else if (currentDirection == 0 && angle < SMOOTHNESS_TH) {
      // New cluster
      LOG(INFO) << "start new cluster";
      currentDirection = sign;
      if (info[0] == -1) {
        // New cluster
        currentCluster.push_back(noLines-1);
        currentCluster.push_back(0);
        // assign cluster index
        info[0] = lineClusters.size();
        info[noLines-1] = lineClusters.size();
      } else if (currentDirection == clusterSign[info[0]]) {
        // Join current cluster with cluster l2's cluster
        currentCluster.push_back(noLines-1);
        currentCluster.insert(currentCluster.end(), lineClusters[info[0]].begin(), lineClusters[info[0]].end());
        lineClusters[info[0]] = currentCluster;
        currentCluster.clear();
      }
	  }

	  if (!currentCluster.empty()) {
	    LOG(INFO) << "end current cluster";
	    // Save current cluster
	    lineClusters.push_back(currentCluster);
	    clusterSign.push_back(currentDirection);
	    currentCluster.clear();
	    currentDirection = 0;
	  }
	}
*/