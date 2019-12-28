#pragma once
#include "Defines/h_print.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include  <numeric>
#include <vector>

static std::string vAsString(std::vector<double>& v, int digits = 16, char separator = 0)
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
static void vSet(std::vector<double>& v,std::initializer_list<double> values) {
	// for (auto it = values.begin(); it != values.end(); ++it)
	for (double value : values)
	{
		// v.push_back(*it);
		v.push_back(value);
	}
}
///---------------------------------------------------------------------------
static double vMin(std::vector<double>& v) { return  *std::min_element(v.begin(), v.end()); }
///---------------------------------------------------------------------------
static double vMax(std::vector<double>& v) { return  *std::max_element(v.begin(), v.end()); }
///---------------------------------------------------------------------------
static double vAverage(std::vector<double>& v)
{
	double c = 0.0;
	for (double value : v)
	{
		c += value;
	}
	return  c / v.size();
}
///---------------------------------------------------------------------------
static double vNorm(std::vector<double>& v)
{
	double c = 0.0;
	for (double value : v)
	{
		c += value * value;
	}
	return  c / v.size();
}


// result = a X + b Y
static void vAdd(
	std::vector<double>& result,
	std::vector<double>& X, std::vector<double>& Y,
	double a = 1.0, double b = 1.0)
{
	size_t c = X.size();
	result.resize(c);
	for (size_t i = 0; i < c; ++i)
		result[i] = a * X[i] + b * Y[i];
}

// return = a X * Y
static double vMult(
	std::vector<double>& X, std::vector<double>& Y,
	double a = 1.0)
{
	size_t c = X.size();
	double result = 0.0;
	for (size_t i = 0; i < c; ++i)
		result += a * X[i] * Y[i];
	return result;
}
// return = yData(x)
static double vInterp(
	std::vector<double>& xData, std::vector<double>& yData,
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


static double vStdev(std::vector<double> const & func)
{
	double mean = std::accumulate(func.begin(), func.end(), 0.0) / func.size();
	double sq_sum = std::inner_product(func.begin(), func.end(), func.begin(), 0.0,
		[](double const & x, double const & y) { return x + y; },
		[mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
	return std::sqrt(sq_sum / ( func.size() - 1 ));
}

static double vStdevRatio(std::vector<double> const & func)
{
	double mean = std::accumulate(func.begin(), func.end(), 0.0) / func.size();
	double sq_sum = std::inner_product(func.begin(), func.end(), func.begin(), 0.0,
		[](double const & x, double const & y) { return x + y; },
		[mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
	return std::sqrt(sq_sum / ( func.size() - 1 ))/mean;
}

static double vMean(std::vector<double> const & func)
{
	return std::accumulate(func.begin(), func.end(), 0.0) / func.size();
}


static double vMedian(std::vector<double> v, double position)
{
	size_t n = int(double(v.size()) * position);
	if (n >= v.size()) n = v.size() - 1;
	nth_element(v.begin(), v.begin() + n, v.end());
	return v[n];
}

// y[] = yData( x[] )
static void vInterp(
	std::vector<double>& xData,
	std::vector<double>& yData,
	std::vector<double>& x,
	std::vector<double>& y,
	bool extrapolate = false)
{
	for (double vv : x)
	{
		double yy = vInterp(xData, yData, vv, extrapolate);
		y.push_back(yy);
	}
}


static void vThreshhold(std::vector<double>& in, std::vector<double>& out, double i_threshold)
{
	out.clear();
	int c = int(in.size());
	for (int i = 0; i < c; i++) {
		out.push_back((abs(in[i]) > i_threshold)?1:0);
	}
}

static struct stats{
	double mean;
	double min;
	double max;
	double stdev;
	double median;
	double stdev_ratio;
	double area;
	double perimeter;
};


static stats vStats(std::vector<double> &v)
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


///===================================================================================
///===================================================================================
///===================  Spline class
///===================================================================================
///===================================================================================
///===================================================================================
template <typename T>
class hspline
{
public:
	hspline() {}

	/** A spline with x and y values */
	hspline(const std::vector<T>& x, const std::vector<T>& y) {
		if (x.size() != y.size()) {
			std::cerr << "X and Y must be the same size " << std::endl;
			return;
		}

		if (x.size() < 3) {
			std::cerr << "Must have at least three points for interpolation" << std::endl;
			return;
		}

		//typedef typename std::vector<T>::difference_type size_type;

		size_t n = y.size() - 1;

		std::vector<T> b(n), d(n), a(n), c(n + 1), l(n + 1), u(n + 1), z(n + 1);
		std::vector<T> h(n + 1);

		l[0] = T(1);
		u[0] = T(0);
		z[0] = T(0);
		h[0] = x[1] - x[0];

		for (size_t i = 1; i < n; i++) {
			h[i] = x[i + 1] - x[i];
			l[i] = T(2 * (x[i + 1] - x[i - 1])) - T(h[i - 1]) * u[i - 1];
			u[i] = T(h[i]) / l[i];
			a[i] = (T(3) / T(h[i])) * (y[i + 1] - y[i]) - (T(3) / T(h[i - 1])) * (y[i] - y[i - 1]);
			z[i] = (a[i] - T(h[i - 1]) * z[i - 1]) / l[i];
		}

		l[n] = T(1);
		z[n] = c[n] = T(0);

		for (size_t j = n - 1; j >= 0; j--) {
			c[j] = z[j] - u[j] * c[j + 1];
			b[j] = (y[j + 1] - y[j]) / T(h[j]) - (T(h[j]) * (c[j + 1] + T(2) * c[j])) / T(3);
			d[j] = (c[j + 1] - c[j]) / T(3 * h[j]);
		}

		for (size_t i = 0; i < n; i++) {
			mElements.push_back(Element(x[i], y[i], b[i], c[i], d[i]));
		}
	}
	virtual ~hspline() {}

	T operator[](const T& x) const {
		return interpolate(x);
	}

	T interpolate(const T& x) const {
		if (mElements.size() == 0) return T();

		//typename std::vector<element_type>::const_iterator it;
		auto it = std::lower_bound(mElements.begin(), mElements.end(), element_type(x));
		if (it != mElements.begin()) {
			--it;
		}
		return it->eval(x);
	}

	std::vector<T> operator[](const std::vector<T>& xx) const {
		return interpolate(xx);
	}

	/* Evaluate at multiple locations, assuming xx is sorted ascending */
	std::vector<T> interpolate(const std::vector<T>& xx) const {
		if (mElements.size() == 0) return std::vector<T>(xx.size());

		typename std::vector<element_type>::const_iterator it2 = mElements.begin();
		std::vector<T> ys;
		for (typename std::vector<T>::const_iterator it = xx.begin(); it != xx.end(); ++it) {
			it2 = std::lower_bound(it2, mElements.end(), element_type(*it));
			if (it2 != mElements.begin()) {
				--it2;
			}
			ys.push_back(it2->eval(*it));
		}

		return ys;
	}

	void interpolate(const std::vector<T>& xx, std::vector<T>& ys) const {
		typename std::vector<element_type>::const_iterator it2 = mElements.begin();
		for (typename std::vector<T>::const_iterator it = xx.begin(); it != xx.end(); ++it) {
			it2 = std::lower_bound(it2, mElements.end(), element_type(*it));
			if (it2 != mElements.begin()) {
				--it2;
			}
			ys.push_back(it2->eval(*it));
		}
		return;
	}

protected:

	class Element {
	public:
		Element(T _x) : x(_x) {}
		Element(T _x, T _a, T _b, T _c, T _d)
			: x(_x), a(_a), b(_b), c(_c), d(_d) {}

		T eval(const T& xx) const {
			T xix(xx - x);
			return a + b * xix + c * (xix * xix) + d * (xix * xix * xix);
		}

		bool operator<(const Element& e) const {
			return x < e.x;
		}
		bool operator<(const T& xx) const {
			return x < xx;
		}

		T x;
		T a, b, c, d;
	};

	typedef Element element_type;
	std::vector<element_type> mElements;
};

///---------------------------------------------------------------------------
///---------------------------------------------------------------------------
///---------------------------------------------------------------------------
// y[] = yData( x[] )
static void vInterpS(
	std::vector<double>& xData,
	std::vector<double>& yData,
	std::vector<double>& x,
	std::vector<double>& y,
	bool extrapolate = false)
{
	hspline<double> s(xData, yData);
	if (extrapolate) {
		s.interpolate(x, y);
	}
	else
	{
		const size_t n = xData.size() - 1;
		const double xmin = xData[0];
		const double xmax = xData[n];

		for (auto v : x)
		{
			if (v <= xmin)	y.push_back(yData[0]);
			else if (v >= xmax)	y.push_back(yData[n]);
			else y.push_back(s.interpolate(v));
		}
	}
}
///---------------------------------------------------------------------------



