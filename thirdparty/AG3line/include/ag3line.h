
#ifndef AG3LINE_H_
#define AG3LINE_H_

#include <vector>

// opencv
#include <opencv2/core/mat.hpp>

/*----------------------------------------------------------------------------*/
/** Rectangle structure: line segment with width.
*/
struct lineag {
	float x1, y1, x2, y2; /* first and second Point3i of the line segment */
};

void ag3line(const cv::Mat &image, std::vector<lineag> &lines, bool control);

#endif  // AG3LINE_H_
