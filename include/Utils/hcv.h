#pragma once

#include "Defines/h_math.h"
#include "Defines/h_defines.h"
#include "Utils/hvector.h"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "hmath.h"
#include <Defines/h_print.h>

namespace hcv
{
	static void VERSION()
	{
		std::cout << "Using OpenCV version " << CV_VERSION << "\n";
	}

	///------------------------------------------------------------------------------------
	///------   Image functions
	///------------------------------------------------------------------------------------
	static bool isImgGood(cv::Mat& image)
	{
		return !image.empty();
	}

	////-------------------------------------------
	static void checkImg(cv::Mat& image)
	{
		if (isImgGood(image))
		{
			std::cout << "Image is good. \n";
		}
		else
		{
			std::cout << " File not loaded \n";
		}
	}

	////-------------------------------------------
	static void printImgSize(cv::Mat& inImg, const std::string s = "")
	{
		if (s != "") std::cout << s << "   ";
		std::cout
			<< "w["
			<< inImg.size().width
			<< "], h[" << inImg.size().height
			<< "], channels[" << inImg.channels()
			<< "]\n";
	}

	////-------------------------------------------
	// keep = 0 dont keep default
	//		> 1  keep until key pressed
	//		< 0  keep for abs(keep) sseconds
	static void showImg(std::string name, cv::Mat image, int keep = 0)
	{
		imshow(name.c_str(), image);
		if (keep > 0)
		{
			while (true)
			{
				int c = cv::waitKey();
				//cout << c << endl;
				if (c == -1 || c == 27 || c == 32 || c == 13) break;
			}
		}
		else
		{
			cv::waitKey(abs(keep));
		}
	}

	////-------------------------------------------
	static void cropImg(cv::Mat& inImg, cv::Mat& outImg, int w1, int w2, int h1, int h2)
	{
		int ww = (w1 < w2) ? w1 : w2;
		int hh = (h1 < h2) ? h1 : h2;
		int w_ = abs(w1 - w2);
		int h_ = abs(h1 - h2);
		if (ww < 0)
		{
			w_ = (w1 > w2) ? w1 : w2;
			ww = 0;
		}
		if (hh < 0)
		{
			h_ = (h1 > h2) ? h1 : h2;
			hh = 0;
		}
		if (ww + w_ > inImg.size().width) { w_ = inImg.size().width - ww; }
		if (hh + h_ > inImg.size().height) { h_ = inImg.size().height - hh; }
		outImg = inImg(cv::Rect(ww, hh, w_, h_)).clone();
	}

	////-------------------------------------------
	static void scaleImg(cv::Mat& inImg, cv::Mat& outImg, double scale)
	{
		resize(inImg, outImg, cv::Size(), scale, scale, cv::INTER_LANCZOS4);
	}

	////-------------------------------------------
	static void scaleImg(cv::Mat& inImg, cv::Mat& outImg, double scaleW, double scaleH)
	{
		resize(inImg, outImg, cv::Size(), scaleW, scaleH, cv::INTER_LANCZOS4);
	}

	////-------------------------------------------
	static void scaleImgDim(cv::Mat& inImg, cv::Mat& outImg, int i_width, int i_height)
	{
		double scale;
		if (i_height == 0)
		{
			scale = double(i_width) / double(inImg.size().width);
			scaleImg(inImg, outImg, scale);
			return;
		}
		if (i_width == 0)
		{
			scale = double(i_height) / double(inImg.size().height);
			scaleImg(inImg, outImg, scale);
		}
		else
		{
			double scaleW = double(i_width) / double(inImg.size().width);
			double scaleH = double(i_height) / double(inImg.size().height);
			scaleImg(inImg, outImg, scaleW, scaleH);
		}
	}

	////-------------------------------------------
	static void scaleImgFit(cv::Mat& inImg, int i_width, int i_height, cv::Mat* outImg = nullptr)
	{
		if (outImg == nullptr) outImg = &inImg;

		const double scale1 = double(i_width) / double(inImg.size().width);
		const double scale2 = double(i_height) / double(inImg.size().height);
		if (scale1 > scale2)
		{
			i_width = int(round(inImg.size().width * scale2));
		}
		else
		{
			i_height = int(round(inImg.size().height * scale1));
		}
		resize(inImg, *outImg, cv::Size(i_width, i_height), 0, 0, cv::INTER_LANCZOS4);
	}

	////-------------------------------------------
	static void blurImg(cv::Mat& inImg, cv::Mat& outImg, int isz = 5)
	{
		// GaussianBlur (   InputArray src, OutputArray dst, Size ksize,
		//                  double sigmaX, double sigmaY=0,
		//                  int borderType=BORDER_DEFAULT )
		GaussianBlur(inImg, outImg, cv::Size(isz, isz), 0.0, 0.0);
	}

	////-------------------------------------------
	static void CannyB(
		cv::Mat& image,
		cv::Mat& edges,
		int gaussian_blur_window_size,
		double threshold1,
		double threshold2,
		int apertureSize,
		bool L2gradient)
	{
		cv::Mat tmp = image.clone();
		blurImg(image, tmp, gaussian_blur_window_size);
		Canny(tmp, edges, threshold1, threshold2, apertureSize, L2gradient);
	}

	////-------------------------------------------
	/// MASKS
	static cv::Mat createMask(cv::Mat& inImg, int border = 0)
	{
		int w = inImg.size().width + 2 * border;
		int h = inImg.size().height + 2 * border;
		cv::Mat mask(cv::Size(w, h), CV_8UC1);
		return mask;
	}

	////-------------------------------------------
	static void invMask(cv::Mat& inMask, cv::Mat* outMask = nullptr)
	{
		if (outMask == nullptr)
		{
			inMask = 255 - inMask;
			//bitwise_not(inMask, inMask);
		}
		else
		{
			*outMask = inMask.clone();
			*outMask = 255 - *outMask;
			//bitwise_not(inMask, *outMask);
		}
	}

	static void addMask(cv::Mat& inImg, cv::Mat& mask, int R = 255, int G = 255, int B = 0)
	{
		std::vector<cv::Mat> rgbChannels(3);
		split(inImg, rgbChannels);
		rgbChannels[0] += mask * B / 255;
		rgbChannels[1] += mask * G / 255;
		rgbChannels[2] += mask * R / 255;
		merge(rgbChannels, inImg);
	}

	static int threshold(cv::Mat& inImg, cv::Mat& outImg, int lower_than)
	{
		int cnt = 0;
		int w = inImg.size().width;
		int h = inImg.size().height;

		outImg = inImg.clone() * 0;

		cv::Point p;
		int v;

		for (auto i = 0; i < w; i++)
		{
			p.x = i;
			for (auto j = 0; j < h; j++)
			{
				p.y = j;
				v = inImg.at<uchar>(p);

				if (v < lower_than)
				{
					outImg.at<uchar>(p) = uchar(0);
					cnt++;
				}
				else outImg.at<uchar>(p) = uchar(v);
			}
		}
		return cnt;
	}

	///------------------------------------------------------------------------------------
	///------   Point functions
	///------------------------------------------------------------------------------------
	static double dist(cv::Point& a, cv::Point& b)
	{
		return pow(pow(b.x - a.x, 2) + pow(b.y - a.y, 2), 0.5);
	}

	////-------------------------------------------
	static double dist(double& x1, double& y1, double& x2, double& y2)
	{
		return pow(pow(x1 - x2, 2) + pow(y1 - y2, 2), 0.5);
	}

	////-------------------------------------------
	static double angle(cv::Point& a, cv::Point& b)
	{
		double x1 = a.x;
		double x2 = b.x;
		double y1 = a.y;
		double y2 = b.y;
		double ang = atan2(y2 - y1, x2 - x1) * 180.0 / CV_PI;
		if (ang < 0.0) ang = 360.0 + ang;
		return ang;
	}

	///------------------------------------------------------------------------------------
	///------   Polygon functions
	///------------------------------------------------------------------------------------

	////-------------------------------------------
	static cv::Point getCenter(std::vector<cv::Point>& pts)
	{
		const size_t c = pts.size();
		double x = 0;
		double y = 0;
		for (size_t i = 0; i < c; i++)
		{
			x += pts[i].x;
			y += pts[i].y;
		}
		return { static_cast<int>(round(x / c)), static_cast<int>(round(y / c)) };
	}

	////-------------------------------------------
	static double getAverageRadius(std::vector<cv::Point>& pts, cv::Point* center = nullptr)
	{
		const size_t c = pts.size();
		cv::Point ctr;
		if (center == nullptr)
			ctr = getCenter(pts);
		else
			ctr = *center;
		double radius = 0.0;
		for (size_t i = 0; i < c; i++)
		{
			radius += dist(pts[i], ctr);
		}
		return radius / c;
	}

	////-------------------------------------------
	static double getMaxRadius(std::vector<cv::Point>& pts, cv::Point* center = nullptr)
	{
		const size_t c = pts.size();
		cv::Point ctr;
		if (center == nullptr)
			ctr = getCenter(pts);
		else
			ctr = *center;

		double radius = 0.0;
		double tr;
		for (size_t i = 0; i < c; i++)
		{
			tr = dist(pts[i], ctr);
			if (tr > radius) radius = tr;
		}
		return radius;
	}

	////-------------------------------------------
	static double getMinRadius(std::vector<cv::Point>& pts, cv::Point* center = nullptr)
	{
		const size_t c = pts.size();
		cv::Point ctr;
		if (center == nullptr)
			ctr = getCenter(pts);
		else
			ctr = *center;

		double radius = 9999999999999.0;
		double tr;
		for (size_t i = 0; i < c; i++)
		{
			tr = dist(pts[i], ctr);
			if (tr < radius) radius = tr;
		}
		return radius;
	}

	////-------------------------------------------
	static void getRadii(std::vector<cv::Point>& pts, std::vector<double>& radii, cv::Point* center = nullptr)
	{
		radii.clear();
		const size_t c = pts.size();
		cv::Point ctr;
		if (center == nullptr)
			ctr = getCenter(pts);
		else
			ctr = *center;

		double radius = 0.0;
		for (size_t i = 0; i < c; i++)
		{
			radius = dist(pts[i], ctr);
			radii.push_back(radius);
		}
	}
	////-------------------------------------------

	static void diffRadii(std::vector<double>& radii, std::vector<double>& diff)
	{
		diff.clear();
		const size_t c = radii.size();

		for (size_t i = 0; i < c - 1; i++)
		{
			diff.push_back(radii[i + 1] - radii[i]);
		}
		diff.push_back(radii[0] - radii[c - 1]);
	}
	////-------------------------------------------

	static void diffPolygon(std::vector<cv::Point>& pts, std::vector<double>& diffX, std::vector<double>& diffY)
	{
		diffX.clear();
		diffY.clear();
		const size_t c = pts.size();

		for (size_t i = 0; i < c - 1; i++)
		{
			diffX.push_back(pts[i + 1].x - pts[i].x);
			diffY.push_back(pts[i + 1].y - pts[i].y);
		}
		diffX.push_back(pts[0].x - pts[c - 1].x);
		diffY.push_back(pts[0].y - pts[c - 1].y);
	}

	////-------------------------------------------
	static double getArea(std::vector<cv::Point>& pts)
	{
		const int n = int(pts.size());
		double area = 0.0;
		int j = n - 1;
		for (int i = 0; i < n; i++)
		{
			area += (pts[j].x + pts[i].x) * (pts[j].y - pts[i].y);
			j = i; // j is previous vertex to i
		}
		// Return absolute value
		return abs(area / 2.0);
	}
	////-------------------------------------------
		static double getPerimeter(std::vector<cv::Point>& pts)
	{
		const int n = int(pts.size());
		double per = 0.0;
		for (int i = 0; i < n-1; i++)
		{
			per += dist(pts[i],pts[i+1]);
		}
		per += dist(pts[0],pts[n-1]);
		return per;
	}

	////-------------------------------------------
	static void featherPoints(std::vector<cv::Point>& pts, cv::Point& center, int border)
	{
		size_t c = pts.size();
		double ang;
		double radius;
		int r;
		cv::Point p;

		for (size_t i = 0; i < c; i++)
		{
			ang = round(angle(center, pts[i]));
			radius = dist(center, pts[i]);
			r = int(round(radius)) + border;

			p.x = int(center.x + r * dcos(ang));
			p.y = int(center.y + r * dsin(ang));
			pts[i] = p;
		}
	}

	////-------------------------------------------
	// LinearCurveLength calculates the total length of the linear
	// interpolation through a vector of Points.  It is the sum of
	// the Euclidean distances between all consecutive points in
	// the vector.
	static double getLength(std::vector<cv::Point>& points, bool closed = false)
	{
		size_t c = points.size() - 1;
		double sum = 0.0;
		for (size_t i = 0; i < c; ++i)
		{
			sum += dist(points[i], points[i + 1]);
		}
		if (closed) sum += dist(points[0], points[c]);
		return sum;
	}

	////-------------------------------------------
	static void interpolatePolygon(
		std::vector<cv::Point>& source,
		std::size_t result_count,
		bool close = true,
		std::vector<cv::Point>* result = nullptr)
	{
		bool flag = false;
		std::vector<cv::Point> tmp;
		cv::Point p;
		if (result == nullptr)
		{
			result = &tmp;
			flag = true;
		}
		if (close) source.push_back(source[0]);
		result->clear();

		if (source.size() < 2 || result_count < 2)
		{
			// degenerate source vector or result_count value
			// for simplicity, this returns an empty result
			// but special cases may be handled when appropriate for the application
			return;
		}

		const double total_length = getLength(source, close);
		const double segment_length = total_length / (result_count - 1);

		// start and finish are the current source segment's endpoints
		auto start = source.begin();
		auto finish = start + 1;

		double src_segment_offset = 0;
		double src_segment_length = dist(*start, *finish);

		// The first point in the result is the same as the first point
		// in the source.
		result->push_back(*start);
		double dd;

		for (std::size_t i = 1; i < result_count - 1; ++i)
		{
			const double next_offset = segment_length * i;

			while (src_segment_offset + src_segment_length < next_offset)
			{
				src_segment_offset += src_segment_length;
				start = finish++;
				src_segment_length = dist(*start, *finish);
			}

			const double part_offset = next_offset - src_segment_offset;
			const double part_ratio = part_offset / src_segment_length;

			dd = finish->x - start->x;
			p.x = start->x + int(round(part_ratio * dd));
			dd = finish->y - start->y;
			p.y = start->y + int(round(part_ratio * dd));
			result->push_back(p);
		}
		result->push_back(source.back());

		if (close)
		{
			source.pop_back();
			result->pop_back();
		}

		if (flag)
		{
			std::vector<cv::Point> out;
			source = tmp;
		}
	}

	////-------------------------------------------
	static bool pointInPolygon(cv::Point& point, std::vector<cv::Point>& polygon)
	{
		std::vector<cv::Point> points = polygon;
		size_t i, j;
		size_t nvert = points.size();
		bool c = false;
		for (i = 0, j = nvert - 1; i < nvert; j = i++)
		{
			if (((points[i].y >= point.y) != (points[j].y >= point.y)) &&
				(point.x <= (points[j].x - points[i].x) *
				(point.y - points[i].y) / (points[j].y - points[i].y)
					+ points[i].x))
				c = !c;
		}
		return c;
	}

	//////-------------------------------------------
	//static void fillPolygon(cv::Mat& img, std::vector<cv::Point>& polygon)
	//{
	//	int xmin = 9999999;
	//	int ymin = 9999999;
	//	int xmax = -9999999;
	//	int ymax = -9999999;
	//	size_t n = polygon.size();
	//	for (size_t i = 0; i < n; i++)
	//	{
	//		if (polygon[i].x < xmin) xmin = polygon[i].x;
	//		if (polygon[i].y < ymin) ymin = polygon[i].y;
	//		if (polygon[i].x > xmax) xmax = polygon[i].x;
	//		if (polygon[i].y > ymax) ymax = polygon[i].y;
	//	}
	//	cv::Point p;
	//	for (int x = xmin; x <= xmax; x++)
	//		for (int y = ymin; y <= ymax; y++)
	//		{
	//			p.x = x;
	//			p.y = y;
	//			if (pointInPolygon(p, polygon)) img.at<uchar>(p) = uchar(255);
	//		}
	//}

	////-------------------------------------------
	static void shrinkLine(cv::Point& start, cv::Point& end, double ratio)
	{
		end.x = int(start.x + (end.x - start.x) * ratio);
		end.y = int(start.y + (end.y - start.y) * ratio);
	}

	////-------------------------------------------
	static double roundPolygon(
		std::vector<cv::Point>& pts, cv::Point& center,
		std::vector<double>& radii_output, double median = 0.9)
	{
		getRadii(pts, radii_output, &center);
		double m = vMedian(radii_output, median);

		size_t c = pts.size();
		int cc = 0;
		for (size_t i = 0; i < c; ++i)
		{
			if (radii_output[i] > m)
			{
				double ratio = m / radii_output[i];
				shrinkLine(center, pts[i], ratio);
				cc++;
			}
		}
		return vStdevRatio(radii_output);
	}

	///------------------------------------------------------------------------------------
	///------   Draw functions
	///------------------------------------------------------------------------------------

	static void drawLine(cv::Mat& inImg, int x1, int y1, int x2, int y2,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		line(inImg, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
	}

	////-------------------------------------------
	static void drawLine(cv::Mat& inImg, cv::Point& p1, cv::Point& p2,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		line(inImg, p1, p2, cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
	}

	////-------------------------------------------
	static void drawLineNAA(cv::Mat& inImg, cv::Point& p1, cv::Point& p2,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		line(inImg, p1, p2, cv::Scalar(B, G, R), i_line_width);
	}

	////-------------------------------------------
	static void drawCircle(cv::Mat& inImg, int x, int y, int radius,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		circle(inImg, cv::Point(x, y), radius, cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
	}

	////-------------------------------------------
	static void drawCircle(cv::Mat& inImg, cv::Point& center, int radius,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		circle(inImg, center, radius, cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
	}

	////-------------------------------------------
	static void drawPoint(cv::Mat& inImg, cv::Point& point,
		int R = 255, int G = 255, int B = 0)
	{
		inImg.at<cv::Vec3b>(point) = cv::Vec3b(B, G, R);
		//inImg.at<uchar>(r.y, r.x, 0) = B;
		//inImg.at<uchar>(r.y, r.x, 1) = B;
		//inImg.at<uchar>(r.x, r.y, 2) = R;
	}

	////-------------------------------------------
	static void drawLines(cv::Mat& inImg, cv::Point& center, std::vector<cv::Point>& points,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		size_t c = points.size();
		for (size_t i = 0; i < c; i++)
		{
			drawLine(inImg, center, points[i], i_line_width, R, G, B);
		}
	}

	////-------------------------------------------
	static void drawPointsC(cv::Mat& inImg, std::vector<cv::Point>& points,
		int radius = 1,
		int i_line_width = 1, int R = 0, int G = 0, int B = 255)
	{
		size_t c = points.size();
		for (size_t i = 0; i < c; i++)
			drawCircle(inImg, points[i], radius, i_line_width, R, G, B);
	}

	////-------------------------------------------
	static void drawPoints(cv::Mat& inImg, std::vector<cv::Point>& points,
		int R = 255, int G = 255, int B = 0)
	{
		size_t c = points.size();
		for (size_t i = 0; i < c; i++)
			inImg.at<cv::Vec3b>(points[i]) = cv::Vec3b(B, G, R);
	}

	////-------------------------------------------
	static void drawPolygon(cv::Mat& inImg, std::vector<cv::Point>& points,
		bool closed = true,
		int i_line_width = 1, int R = 0, int G = 255, int B = 0)
	{
		size_t c = points.size();
		if (c > 0)
		{
			c = c - 1;
			for (size_t i = 0; i < c; i++)
				drawLine(inImg, points[i], points[i + 1], i_line_width, R, G, B);
			if (closed) drawLine(inImg, points[c], points[0], i_line_width, R, G, B);
		}
	}

	////-------------------------------------------
	static void drawPolygonNAA(cv::Mat& inImg, std::vector<cv::Point>& points,
		bool closed = true,
		int i_line_width = 1, int R = 0, int G = 255, int B = 0)
	{
		size_t c = points.size();
		if (c > 0)
		{
			c = c - 1;
			for (size_t i = 0; i < c; i++)
				drawLineNAA(inImg, points[i], points[i + 1], i_line_width, R, G, B);
			if (closed) drawLine(inImg, points[c], points[0], i_line_width, R, G, B);
		}
	}

	////-------------------------------------------

	static void drawText(cv::Mat& inImg, std::string txt,
		int x, int y, double font_size,
		int i_line_width = 1, int R = 0, int G = 0, int B = 0)
	{
		putText(
			inImg, //target image
			txt.c_str(), //text
			cv::Point(x, y), //top-left position
			cv::FONT_HERSHEY_DUPLEX,
			font_size,
			CV_RGB(R, G, B), //font color
			i_line_width,
			cv::LINE_AA);
	}

	static void fillPolygon(cv::Mat& inImg, std::vector<cv::Point>& pts,
		int R = 255, int G = 255, int B = 255)
	{
		cv::fillConvexPoly(inImg, pts.data(), int(pts.size()), cv::Scalar(B, G, R));
	}
	///------------------------------------------------------------------------------------
	///------   Utilities
	///------------------------------------------------------------------------------------

	static void printPointVector(std::vector<cv::Point>& pts)
	{
		for (cv::Point p : pts)
		{
			pr2c(p.x, p.y);
		}
	}

	/// average (2*dist+1) points
	static void movingAveragePointVector(std::vector<cv::Point>& pts, int dist, std::vector<cv::Point>& out)
	{
		int n = int(pts.size());
		out.clear();

		for (int i = 0; i < n; i++)
		{
			cv::Point   pcum(0, 0);

			for (int j = -dist; j <= dist; j++)
			{
				auto k = i + j;
				citr(k, n);
				pcum += pts[k];
			}
			pcum = pcum / (2 * dist + 1);
			out.push_back(pcum);
		}
	}

	/// average (2*dist+1) points
	static std::vector<cv::Point> movingAveragePointVector(std::vector<cv::Point>& pts, int dist)
	{
		int n = int(pts.size());
		std::vector<cv::Point> out;

		for (int i = 0; i < n; i++)
		{
			cv::Point   pcum(0, 0);

			for (int j = -dist; j <= dist; j++)
			{
				auto k = i + j;
				citr(k, n);
				pcum += pts[k];
			}
			pcum = pcum / (2 * dist + 1);
			out.push_back(pcum);
		}
		return out;
	}
	///------------------------------------------------------------------------------------
	///------   ray
	///------------------------------------------------------------------------------------

	class hcvRay
	{
	public:
		std::vector<int> ray;
		int im{}, jm{};
		int pixel_count{};
		cv::Point ray_end_point;
		int ray_inc, angle_inc;
		int tolerance;

		std::vector<cv::Point> grain_pts;
		std::vector<double> grain_radii;
		cv::Point grain_center;
		stats grain_stats;
		std::vector<stats> grains;

		// no. of points to look back at
		// back_span = 5;
		int back_span;

		cv::Mat* img{};
		cv::Mat img_tmp;
		cv::Mat mask;
		cv::Mat grain_mask;

		std::vector<cv::Point> pointVector;

		hcvRay()
		{
			ray.reserve(500);
			tolerance = 25;
			back_span = 0;
			pointVector.reserve(1000000);
			ray_inc = 1;
			angle_inc = 5;
		}

		~hcvRay() = default;

		void addGrain()
		{
			grain_stats.area = getArea(grain_pts);
			grain_stats.perimeter = getPerimeter(grain_pts);
			grains.push_back(grain_stats);
		}

		void printGrains()
		{
			for (stats s : grains)
			{
				pr6c(s.area,s.perimeter,s.mean,s.min,s.max,s.stdev_ratio);
			}
		}




		void setImage(cv::Mat* inImg)
		{
			img = inImg;
			im = inImg->size().width - 1;
			jm = inImg->size().height - 1;
			//hcv::getMask(*img, mask, 0);
			mask = createMask(*img, 0);

			//hcv::getMask(*img, grain_mask, 0);
			grain_mask = createMask(*img, 0);

			grain_mask = grain_mask * 0;
			mask = mask * 0;
		}

		bool checkPointInMask(cv::Point& p)
		{
			return mask.at<uchar>(p) != uchar(255);
		}

		bool pointInsideImage(cv::Point& p) const
		{
			if (p.x <= 0) return false;
			if (p.x > im) return false;
			if (p.y <= 0) return false;
			if (p.y > jm) return false;
			return true;
		}

		void featherGrainMask(std::vector<cv::Point>& pts, cv::Point& center, int border = 3)
		{
			featherPoints(pts, center, border);
			grain_mask = grain_mask * 0;
			//drawPolygonNAA(grain_mask, pts, true, 1, 255, 255, 255);
			//floodFill(center, grain_mask, grain_mask);
			//cv::fillConvexPoly(grain_mask, pts.data(), pts.size(), cv::Scalar(255,255,255));
			fillPolygon(grain_mask, pts);
		}

		void mergeGrainMask()
		{
			mask = mask + grain_mask;
		}

		///------------------------------------------------------------------------------
#define TOLORENCE_START_RANGE 10
#define TOLORENCE_STEP 5

		bool nfill2(cv::Point p, double threshhold, int Min_GRAIN_SIZE = 4, int Max_GRAIN_SIZE = 400)
		{
			int opc = 0; // old pixel count
			int cnt = 0;

			for (int i = TOLORENCE_START_RANGE; i < 255; i = i + TOLORENCE_STEP)
			{
				cv::Mat mask_tmp = mask.clone();
				img_tmp = img->clone();
				img_tmp += mask;
				tolerance = i;
				floodFill(p, img_tmp, mask_tmp);
				cnt++;
				if (opc != 0)
				{
					if (pixel_count <= 4)
					{
						//pr5("fill cycles = ", cnt, "       Pixel count = ", pixel_count, "  (EMPTY)");
						return false;
					}

					const double ratio = double(pixel_count - opc) / double(opc);
					if (ratio > threshhold)
					{
						tolerance = tolerance - TOLORENCE_STEP;
						img_tmp = img->clone();
						img_tmp += mask;
						grain_mask = grain_mask * 0;
						floodFill(p, img_tmp, grain_mask);
						cnt++;

						if (pixel_count < Min_GRAIN_SIZE || pixel_count > Max_GRAIN_SIZE)
						{
							grain_mask = grain_mask * 0;
#ifdef HDEBUG
							if (pixel_count < Min_GRAIN_SIZE)
							{
								pr5("fill cycles = ", cnt, "       Pixel count = ", pixel_count,
									"  (GRAIN is too small)");
							}
							else
							{
								pr5("fill cycles = ", cnt, "       Pixel count = ", pixel_count,
									"  (GRAIN is too large)");
							}
#endif
							return false;
						}
#ifdef HDEBUG
						pr4("fill cycles = ", cnt, "       Pixel count = ", pixel_count);
#endif
						return true;
					}
				}
				opc = pixel_count;
			}
			return false;
		}

		///------------------------------------------------------------------------------
		bool nfill(cv::Point p, int i_tolerance, int Min_GRAIN_SIZE = 4, int Max_GRAIN_SIZE = 400)
		{
			cv::Mat mask_tmp = mask.clone();
			img_tmp = img->clone();
			img_tmp += mask;

			if (i_tolerance >= 0) tolerance = i_tolerance;

			tolerance = 2;

			grain_mask = grain_mask * 0;
			floodFill(p, img_tmp, grain_mask);

			if (pixel_count < Min_GRAIN_SIZE || pixel_count > Max_GRAIN_SIZE)
			{
				grain_mask = grain_mask * 0;
#ifdef HDEBUG
				if (pixel_count < Min_GRAIN_SIZE)
				{
					pr5("fill cycles = ", cnt, "       Pixel count = ", pixel_count,
						"  (GRAIN is too small)");
				}
				else
				{
					pr5("fill cycles = ", cnt, "       Pixel count = ", pixel_count,
						"  (GRAIN is too large)");
				}
#endif
				return false;
			}
#ifdef HDEBUG
			pr4("fill cycles = ", cnt, "       Pixel count = ", pixel_count);
#endif
			return true;
		}
		///------------------------------------------------------------------------------
		void floodFill(cv::Point p, cv::Mat& img_source, cv::Mat& mask)
		{
			int pxl;
			cv::Point pr;
			pointVector.clear();
			pointVector.push_back(p);
			size_t flag = pointVector.size();
			pixel_count = 0;
			while (flag)
			{
				p = pointVector[flag - 1];
				pointVector.pop_back();
				if (pointInsideImage(p))
				{
					img_source.at<uchar>(p) = uchar(255);
					if (mask.at<uchar>(p) < 10)
					{
						pixel_count++;
						mask.at<uchar>(p) = uchar(255);
					}
					////----------------------------------------------------------------------------------------------
					pr.x = p.x - 1;					pr.y = p.y;
					if (pointInsideImage(pr))
					{
						pxl = img_source.at<uchar>(pr); if (pxl < tolerance) pointVector.push_back(pr);
					}
					////----------------------------------------------------------------------------------------------
					pr.x = p.x + 1;					pr.y = p.y;
					if (pointInsideImage(pr))
					{
						pxl = img_source.at<uchar>(pr); if (pxl < tolerance) pointVector.push_back(pr);
					}
					//------------------------------------------------------------------------------------------------
					pr.x = p.x;					pr.y = p.y - 1;
					if (pointInsideImage(pr))
					{
						pxl = img_source.at<uchar>(pr); if (pxl < tolerance) pointVector.push_back(pr);
					}
					////----------------------------------------------------------------------------------------------
					pr.x = p.x;					pr.y = p.y + 1;
					if (pointInsideImage(pr))
					{
						pxl = img_source.at<uchar>(pr); if (pxl < tolerance) pointVector.push_back(pr);
					}
					//------------------------------------------------------------------------------------------------
				}
				flag = pointVector.size();
			}
		}
		///------------------------------------------------------------------------------

			//void floodFill(cv::Point p, cv::Mat& img_source);
			//void floodFill(cv::Point p);
			//void rfill(cv::Point p);
			//void hfill(cv::Point& p, cv::Point& pl, cv::Point& pr, int& l, cv::Mat* mask = nullptr);
			//void fill(cv::Mat& mask, cv::Point p);
			//bool nfill3(cv::Point& p, double threshhold, std::vector<cv::Point>& pts, int Min_GRAIN_SIZE = 4);

		bool getRay(cv::Point center, double angle, cv::Mat* inImg = nullptr)
		{
			double c = 0;
			int cc = -1;
			if (inImg == nullptr) inImg = img;

			double di = dcos(angle);
			double dj = dsin(angle);
			double ic = double(center.x);
			double jc = double(center.y);
			// clear back propagation vector
			if (back_span != 0) ray.clear();

			for (;;)
			{
				cc++;
				c += ray_inc;
				int i = int(round(ic + c * di));
				int j = int(round(jc + c * dj));

				if (i < 0 || i >= im || j < 0 || j >= jm)
				{
					if (i < 0) i = 0;
					if (i >= im) i = im - 1;
					if (j < 0)j = 0;
					if (j >= jm) j = jm - 1;
					ray_end_point = cv::Point(i, j);
					return true;
				}
				int pxl = inImg->at<uchar>(j, i);

				//// check back propagation vector
				//if (back_span != 0) {
				//	ray.push_back(pxl);
				//	int cpxl = getBackAverage(cc);
				//	if (pxl - cpxl > tolerance) {
				//		ray_end_point = cv::Point(i, j);
				//		return true;
				//	}
				//}
				// no back propagation
				//else {
				if (pxl > tolerance)
				{
					ray_end_point = cv::Point(i, j);
					return true;
				}
				//}
			}
		}

		bool getRayInv(cv::Point center, double angle, cv::Mat* inImg = nullptr)
		{
			double c = 0;
			int cc = -1;
			if (inImg == nullptr) inImg = &grain_mask;

			double di = dcos(angle);
			double dj = dsin(angle);
			double ic = double(center.x);
			double jc = double(center.y);
			// clear back propagation vector
			if (back_span != 0) ray.clear();

			for (;;)
			{
				cc++;
				c += ray_inc;
				int i = int(round(ic + c * di));
				int j = int(round(jc + c * dj));

				if (i < 0 || i >= im || j < 0 || j >= jm)
				{
					if (i < 0) i = 0;
					if (i >= im) i = im - 1;
					if (j < 0)j = 0;
					if (j >= jm) j = jm - 1;
					ray_end_point = cv::Point(i, j);
					return true;
				}
				const int pxl = inImg->at<uchar>(j, i);
				if (pxl <= 0)
				{
					ray_end_point = cv::Point(i, j);
					return true;
				}
			}
		}

		void getRaysInv(const cv::Point& center, std::vector<cv::Point>& points,
			cv::Mat* inImg = nullptr)
		{
			points.clear();
			for (int angle = 0; angle < 360; angle += angle_inc)
			{
				if (getRayInv(center, angle, inImg))
				{
					points.push_back(ray_end_point);
				}
			}
		}

		void getRays(cv::Point center, std::vector<cv::Point>& points,
			cv::Mat* inImg = nullptr)
		{
			points.clear();
			for (auto angle = 0; angle < 360; angle += angle_inc)
			{
				if (getRay(center, angle, inImg))
				{
					points.push_back(ray_end_point);
				}
			}
		}
		///------------------------------------------------------------------------------
		bool getGrain(
			cv::Point& p,
			double fill_tol = 0.3,
			int min_grain_size = 200,
			int max_grain_size = 400,
			int feather_grain_pixels = 0,
			int moving_average_distance = 3,
			int interpolate_polygon_count = 12)
		{
			if (checkPointInMask(p)) {
				const bool fill_result = nfill2(p, fill_tol, min_grain_size, max_grain_size);
				if (fill_result)
				{
					getRaysInv(p, grain_pts);
					hcv::interpolatePolygon(grain_pts, interpolate_polygon_count, true);
					p = hcv::getCenter(grain_pts);

					///---------- Moving Average
					if (moving_average_distance != 0)
					{
						grain_pts = movingAveragePointVector(grain_pts, moving_average_distance);
						p = hcv::getCenter(grain_pts);
					}

					///---------- Feather Grain
					if (feather_grain_pixels != 0)
					{
						featherPoints(grain_pts, p, feather_grain_pixels);
					}

					grain_mask = grain_mask * 0;
					fillPolygon(grain_mask, grain_pts);
					grain_center = p;
					getRadii(grain_pts,grain_radii,&p);
					grain_stats = vStats(grain_radii);
				}
				return fill_result;
			}
			return false;
		}

		//int getBackAverage(int index);
		//int getBackMax(int index);
		//int getBackMin(int index);
	};
}
