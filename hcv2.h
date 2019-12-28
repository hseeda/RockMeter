#include <iostream>
#include <initializer_list>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include "Defines/h_print.h"
#include "Defines/h_defines.h"
#include "Defines/h_math.h"
#include "Utils/hutil.h"


namespace py = pybind11;

#define pyImg  py::array_t<uint8_t>
#define npImg  np<uint8_t>
#define Img  np<uint8_t>

#define vPoint  std::vector<Point>
#define vDouble  std::vector<double>
#define vInt  std::vector<int>
#define PI 3.141592653589793238

///=========================================================================================
///-----------------------------------------------------------------------------------------
struct stats {double mean, min, max, stdev, median, stdev_ratio, area, perimeter; };

std::string vAsString(vDouble &v, int digits = 16, char separator = 0)
{
	std::stringstream s;
	const auto c = v.size() - 1;
	if (separator == 0)
	{
		for (size_t i = 0; i < c; ++i)
			s << std::setw(digits) << v[i];
		s << std::setw(digits) << v[c];
	}
	else
	{
		for (size_t i = 0; i < c; ++i)
			s /*<< std::showpos*/ << std::setw(digits) << v[i] << separator;
		s /*<< std::showpos*/ << std::setw(digits) << v[c];
	}
	return s.str();
}
///---------------------------------------------------------------------------
void vSet(vDouble &v, std::initializer_list<double> values) {
	// for (auto it = values.begin(); it != values.end(); ++it)
	for (double value : values) v.push_back(value);
}
///---------------------------------------------------------------------------
double vMin(vDouble &v) { return  *std::min_element(v.begin(), v.end()); }
///---------------------------------------------------------------------------
static double vMax(vDouble v) { return  *std::max_element(v.begin(), v.end()); }
///---------------------------------------------------------------------------
static double vAverage(vDouble &v)
{
	double c = 0.0;
	for (double value : v) c += value;
	return  c / v.size();
}
///---------------------------------------------------------------------------
static double vNorm(vDouble &v)
{
	double c = 0.0;
	for (double value : v) c += value * value;
	return  c / v.size();
}

// result = a X + b Y
static void vAdd(
	vDouble &result,
	vDouble &X,
	vDouble &Y,
	double a = 1.0, double b = 1.0)
{
	size_t c = X.size();
	result.resize(c);
	for (size_t i = 0; i < c; ++i) result[i] = a * X[i] + b * Y[i];
}

// return = a X * Y
static double vMult(
	vDouble &X, vDouble &Y,
	double a = 1.0)
{
	size_t c = X.size();
	double result = 0.0;
	for (size_t i = 0; i < c; ++i) result += a * X[i] * Y[i];
	return result;
}
// return = yData(x)
static double vInterp(
	vDouble &xData,
	vDouble &yData,
	double x, bool extrapolate = false)
{
	size_t size = xData.size();
	size_t i = 0;                                      // find left end of interval for interpolation
	double xL, yL, xR, yR, dydx;
	if (x <= xData[0])
	{
		if (extrapolate)
		{
			xL = xData[0];	yL = yData[0];
			xR = xData[1];	yR = yData[1];
			dydx = (yR - yL) / (xR - xL);
			return yL + dydx * (x - xL);
		}
		else
		{
			return yData[0];
		}
	}
	else if (x >= xData[size - 1])
	{
		if (extrapolate)
		{
			xL = xData[size - 2];	yL = yData[size - 2];
			xR = xData[size - 1];	yR = yData[size - 1];
			dydx = (yR - yL) / (xR - xL);
			return yL + dydx * (x - xL);
		}
		else
		{
			return yData[size - 1];
		}
	}
	else
	{
		i = 1;
		while (x > xData[i]) i++;
		xL = xData[i - 1];		yL = yData[i - 1];
		xR = xData[i];	yR = yData[i];
		dydx = (yR - yL) / (xR - xL);
		return yL + dydx * (x - xL);
	}
}

static double vStdev(vDouble const &func)
{
	double mean = std::accumulate(func.begin(), func.end(), 0.0) / func.size();
	double sq_sum = std::inner_product(func.begin(), func.end(), func.begin(), 0.0,
		[](double const & x, double const & y) { return x + y; },
		[mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
	return std::sqrt(sq_sum / ( func.size() - 1 ));
}

static double vStdevRatio(vDouble const &func)
{
	double mean = std::accumulate(func.begin(), func.end(), 0.0) / func.size();
	double sq_sum = std::inner_product(func.begin(), func.end(), func.begin(), 0.0,
		[](double const & x, double const & y) { return x + y; },
		[mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
	return std::sqrt(sq_sum / ( func.size() - 1 ))/mean;
}

static double vMean(vDouble const &func)
{
	return std::accumulate(func.begin(), func.end(), 0.0) / func.size();
}


static double vMedian(vDouble &v, double position)
{
	size_t n = int(double(v.size()) * position);
	if (n >= v.size()) n = v.size() - 1;
	nth_element(v.begin(), v.begin() + n, v.end());
	return v[n];
}

// y[] = yData( x[] )
static void vInterp(
	vDouble &xData,
	vDouble &yData,
	vDouble &x,
	vDouble &y,
	bool extrapolate = false)
{
	for (double vv : x)
	{
		double yy = vInterp(xData, yData, vv, extrapolate);
		y.push_back(yy);
	}
}


static void vThreshold(vDouble &in, vDouble &out, double i_threshold)
{
	out.clear();
	int c = int(in.size());
	for (int i = 0; i < c; i++) {
		out.push_back((abs(in[i]) > i_threshold)?1:0);
	}
}

static stats vStats(vDouble &v)
{
	stats s;
	s.mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
	double mean = s.mean;
	//pr(mean);
	//prvs(v,v.size(),"\n");

	double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0,
		[](double const & x, double const & y) { return x + y; },
		[mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
	s.stdev = std::sqrt(sq_sum / ( v.size() - 1 ));

	s.min = *std::min_element(v.begin(), v.end());
	s.max = *std::max_element(v.begin(), v.end());
	s.median = vMedian(v,0.5);
	s.stdev_ratio = s.stdev / s.mean;
	return s;
}

///=========================================================================================
///-----------------------------------------------------------------------------------------
struct Point{
    Point(){}
    Point(int ix, int iy):x(ix), y(iy) {}
    int x, y;
    Point& operator = (const Point& pt){x = pt.x; y = pt.y; return *this; }
    Point& operator += (const Point& pt){x += pt.x; y += pt.y; return *this; }
    Point& operator / (const double n){ 
         x = (int)( (double) x / n );
         y = (int)( (double) y / n );
         return *this;
         }
    };
///=========================================================================================
///-----------------------------------------------------------------------------------------
template<typename T>
class np
{
    public:
        T *ptr;
        struct pybind11::buffer_info buffer;
        int sz;
        np(){}

        np(py::array_t<T> &arr){
            init(arr);
        }

        void init(py::array_t<T> &arr){
            buffer = arr.request();
            ptr = (T*) buffer.ptr;
            switch (buffer.ndim)
            { case 1: sz = int (buffer.shape[0]); break;
              case 2: sz = int (buffer.shape[0] * buffer.shape[1]); break;
              case 3: sz = int (buffer.shape[0] * buffer.shape[1] * buffer.shape[2]) ; break;
            }
        }

        ~np(){}
        void info(){
                std::cout << "NDim (" << buffer.ndim << ")";
                switch(buffer.ndim)
                {
                    case 1: std::cout << " Size(" << buffer.shape[0] << ")"; break;
                    case 2: std::cout << " Size(" << buffer.shape[0] << "," << buffer.shape[1] << ")"; break;
                    case 3: std::cout << " Size(" << buffer.shape[0] << "," << buffer.shape[1] << "," << buffer.shape[2] << ")"; break;
                }
                switch(buffer.ndim)
                {
                    case 1: std::cout << " Strides(" << buffer.strides[0] << ")\n"; break;
                    case 2: std::cout << " Strides(" << buffer.strides[0] << "," << buffer.strides[1] << ")\n"; break;
                    case 3: std::cout << " Strides(" << buffer.strides[0] << "," << buffer.strides[1] << "," << buffer.strides[2] << ")\n"; break;
                }

        }
        T get (int i){ return *(ptr + i * buffer.strides[0]);}
        T get (int i, int j){ return *(ptr + i * buffer.strides[0] + j * buffer.strides[1]);}
        T get (int i, int j, int k){ return *(ptr + i * buffer.strides[0] + j * buffer.strides[1] + k * buffer.strides[2]);}
        T get (Point p){ return *(ptr + p.y * buffer.strides[0] + p.x * buffer.strides[1]); }
        T get (Point p, int z){ return *(ptr + p.y * buffer.strides[0] + p.x * buffer.strides[1] + z * buffer.strides[2]); }

        void set (T value, int i){ *(ptr + i * buffer.strides[0]) = value; }
        void set (T value, int i, int j){ *(ptr + i * buffer.strides[0] + j * buffer.strides[1]) = value; }
        void set (T value, int i, int j, int k){ *(ptr + i * buffer.strides[0] + j * buffer.strides[1] + k * buffer.strides[2]) = value; }
        void set (T value, Point p)        { *(ptr + p.y * buffer.strides[0] + p.x * buffer.strides[1]) = value; }
        void set (T value, Point p, int z) { *(ptr + p.y * buffer.strides[0] + p.x * buffer.strides[1] + z * buffer.strides[2]) = value; }


        bool in(Point p){
            if( p.y < 0 || p.y >= buffer.shape[0]) return false;
            if( p.x < 0 || p.x >= buffer.shape[1]) return false;
            return true;
        }
//        int size()
//        {
//            switch (buffer.ndim)
//            { case 1: return d0; case 2: return d0 * d1; case 3: return d0 * d1 * d2 ; }
//            return 0;
//        }
        void copy(np<T> &source)  { for (int i=0; i< sz; ++i) *(ptr+i) = *(source.ptr+i); }
        void or(np<T> &source)   { for (int i=0; i< sz; ++i)  *(ptr+i) |= *(source.ptr+i);}
        void add(np<T> &source)   { for (int i=0; i< sz; ++i)  *(ptr+i) += *(source.ptr+i);}
                // pr2c( int (*(ptr+i)) , int(*(source.ptr+i)) )


        void mult (np<T> &source) { for (int i=0; i< sz; ++i) *(ptr+i) = *(ptr+i) * *(source.ptr+i); }
        void mult (int n)         { for (int i=0; i< sz; ++i) *(ptr+i) = *(ptr+i) * n; }
        void zero()               { for (int i=0; i< sz; ++i) *(ptr+i) = 0; }
};
///=========================================================================================
///-----------------------------------------------------------------------------------------
#define INF 10000
bool onSegment(Point &p, Point &q, Point &r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
            q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;
    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point &p, Point &q, Point &r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // co-linear
    return (val > 0)? 1: 2; // clock or counter-clock wise
}

// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(Point &p1, Point &q1, Point &p2, Point &q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

// Returns true if the point p lies inside the polygon[] with n vertices
bool isInside(Point &p , std::vector<Point> &polygon)
{
    size_t n = polygon.size();
    // There must be at least 3 vertices in polygon[]
    if (n < 3)  return false;

    // Create a point for line segment from p to infinite
    Point extreme = {INF, p.y};

    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;
    do
    {
        int next = (i+1)%n;

        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(polygon[i], polygon[next], p, extreme))
        {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (orientation(polygon[i], p, polygon[next]) == 0)
               return onSegment(polygon[i], p, polygon[next]);

            count++;
        }
        i = next;
    } while (i != 0);

    // Return true if count is odd, false otherwise
    return count&1;  // Same as (count%2 == 1)
}
//=========================================================================
int mask_floodFill(
    Point p,
    npImg &mask,
    int tolerance = 64)
{
	int pxl;
	Point pr;
    std::vector<Point> pointVector;
	pointVector.reserve(10000);
	pointVector.clear();
	pointVector.push_back(p);
	size_t flag = pointVector.size();
	int pixel_count = 0;
	while (flag)
	{
		p = pointVector[flag - 1];
		pointVector.pop_back();
		if (mask.in(p))
		{
			mask.set(255,p);

			////----------------------------------------------------------------------------------------------
			pr.x = p.x - 1;					pr.y = p.y;
			if (mask.in(pr))
			{
				pxl = mask.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			////----------------------------------------------------------------------------------------------
			pr.x = p.x + 1;					pr.y = p.y;
			if (mask.in(pr))
			{
				pxl = mask.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			//------------------------------------------------------------------------------------------------
			pr.x = p.x;					pr.y = p.y - 1;
			if (mask.in(pr))
			{
				pxl = mask.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			////----------------------------------------------------------------------------------------------
			pr.x = p.x;					pr.y = p.y + 1;
			if (mask.in(pr))
			{
				pxl = mask.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			//------------------------------------------------------------------------------------------------
		}
		flag = pointVector.size();
	}
	return pixel_count;
}

int floodFill(
    Point p,
    Img &img_source, Img &mask,
    int tolerance)
{
	int pxl;
	Point pr;
    vPoint pointVector;
	pointVector.reserve(1000000);
	pointVector.clear();
	pointVector.push_back(p);
	size_t flag = pointVector.size();
	int pixel_count = 0;
	while (flag)
	{
		p = pointVector[flag - 1];
		pointVector.pop_back();
		if (img_source.in(p))
		{
			img_source.set(255,p);
			if (mask.get(p) < 64)
			{
				pixel_count++;
				mask.set(255, p);
			}
			////----------------------------------------------------------------------------------------------
			pr.x = p.x - 1;					pr.y = p.y;
			if (img_source.in(pr))
			{
				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			////----------------------------------------------------------------------------------------------
			pr.x = p.x + 1;					pr.y = p.y;
			if (img_source.in(pr))
			{
				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			//------------------------------------------------------------------------------------------------
			pr.x = p.x;					pr.y = p.y - 1;
			if (img_source.in(pr))
			{
				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			////----------------------------------------------------------------------------------------------
			pr.x = p.x;					pr.y = p.y + 1;
			if (img_source.in(pr))
			{
				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
			}
			//------------------------------------------------------------------------------------------------
		}
		flag = pointVector.size();
	}
	return pixel_count;
}

///------------------------------------------------------------------------------------
///------   Point functions
///------------------------------------------------------------------------------------
static double dist(Point& a, Point& b)
{    return pow(pow(b.x - a.x, 2) + pow(b.y - a.y, 2), 0.5);}
////-------------------------------------------
double dist(double& x1, double& y1, double& x2, double& y2)
{    return pow(pow(x1 - x2, 2) + pow(y1 - y2, 2), 0.5);}
////-------------------------------------------
double angle(Point& a, Point& b)
{
    double ang = atan2(b.y - a.y, b.x - a.x) * 180.0 / PI;
    if (ang < 0.0) ang = 360.0 + ang;
    return ang;
}

///------------------------------------------------------------------------------------
///------   Polygon functions
///------------------------------------------------------------------------------------
Point getCenter(vPoint& pts)
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
double getAverageRadius(vPoint& pts, Point& ctr)
{
    const size_t c = pts.size();
    double radius = 0.0;
    for (size_t i = 0; i < c; i++)
    {
        radius += dist(pts[i], ctr);
    }
    return radius / c;
}
////-------------------------------------------
double getMaxRadius(vPoint& pts, Point& ctr)
{
    const size_t c = pts.size();
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
double getMinRadius(vPoint& pts, Point& ctr)
{
    const size_t c = pts.size();
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
void getRadii(vPoint& pts, vDouble &radii, Point& ctr)
{
    radii.clear();
    const size_t c = pts.size();
    double radius = 0.0;
    for (size_t i = 0; i < c; i++)
    {
        radius = dist(pts[i], ctr);
        radii.push_back(radius);
    }
}
////-------------------------------------------
void diffRadii(vDouble& radii, vDouble& diff)
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
void diffPolygon(vPoint& pts, vDouble& diffX, vDouble& diffY)
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
double getArea(vPoint& pts)
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
double getPerimeter(vPoint& pts)
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
void featherPoints(vPoint& pts, Point& center, int border)
{
    size_t c = pts.size();
    double ang, radius;
    int r;
    Point p;
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
///-------------------------------------------
// LinearCurveLength calculates the total length of the linear
// interpolation through a vector of Points.  It is the sum of
// the Euclidean distances between all consecutive points in
// the vector.
double getLength(vPoint& points, bool close = false)
{
    size_t c = points.size() - 1;
    double sum = 0.0;
    for (size_t i = 0; i < c; ++i)
    {
        sum += dist(points[i], points[i + 1]);
    }
    if (close) sum += dist(points[0], points[c]);
    return sum;
}
////-------------------------------------------
void interpolatePolygon(
    vPoint& source,
    int result_count,
    vPoint &result,
    bool close = true)
{
    Point p;

    if (close) source.push_back(source[0]);

    result.clear();

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
    result.push_back(*start);
    double dd;

    for (int i = 1; i < result_count - 1; ++i)
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
        result.push_back(p);
    }
    result.push_back(source.back());

    if (close)
    {
        source.pop_back();
        result.pop_back();
    }
}
////-------------------------------------------
void interpolatePolygon(vPoint& source, int result_count, bool close = true)
{
    std::vector<Point> result;
    interpolatePolygon(source, result_count, result, close);
    source = result;
}
////-------------------------------------------
bool pointInPolygon(Point& point, vPoint& polygon)
{
    std::vector<Point> points = polygon;
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
////-------------------------------------------
static void shrinkLine(Point& start, Point& end, double ratio)
{
    end.x = int(start.x + (end.x - start.x) * ratio);
    end.y = int(start.y + (end.y - start.y) * ratio);
}

////-------------------------------------------
double roundPolygon(
    vPoint& pts, Point& center,
    vDouble& radii_output, double median = 0.9)
{
    getRadii(pts, radii_output, center);
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
////-------------------------------------------
void fillPolygon(Img &img, vPoint& pts)
    {
        Point p;
        p.x= (pts[0].x + pts[1].x) / 2 + 2;
        p.y= (pts[0].y + pts[1].y) / 2 ;
        if (!isInside (p, pts)) p.x -= 4;
        if (isInside (p, pts))
        {
            mask_floodFill(p,img);
        }
    }

///------------------------------------------------------------------------------------
///------   Utilities
///------------------------------------------------------------------------------------
void printPointVector(vPoint& pts)
{
    for (Point p : pts)
    {
        std::cout << p.x << " , " << p.y;
    }
}
////-------------------------------------------
/// average (2*dist+1) points
void movingAveragePointVector(vPoint& pts, int dist, vPoint& out)
{
    int n = int(pts.size());
    out.clear();
    for (int i = 0; i < n; i++)
    {
        Point pcum(0, 0);

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
////-------------------------------------------
/// average (2*dist+1) points
static vPoint movingAveragePointVector(std::vector<Point>& pts, int dist)
{
    int n = int(pts.size());
    std::vector<Point> out;

    for (int i = 0; i < n; i++)
    {
        Point   pcum(0, 0);

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
