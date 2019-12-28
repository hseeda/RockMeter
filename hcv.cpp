/*cppimport
<%
# cls; p -c "import cppimport; cppimport.imp('hcv')"
import pybind11
import os
cvl = r"C:\Program Files\OpenCV\x64\vc16\lib"
qtl=r"C:\Dev_Programs\Qt5.12.3\5.12.3\msvc2017_64\lib"

libs = []
#for file in os.listdir(cvl):
#    if file.endswith(".lib"):
#        libs.append(os.path.splitext(file)[0])
#for file in os.listdir(qtl):
#    if file.endswith(".lib"):
#        libs.append(os.path.splitext(file)[0])
#libs.append('python37')
#l=[	"kernel32", "user32" ,"gdi32", "winspool", "comdlg32",
#	"advapi32", "shell32", "ole32", "oleaut32",
#	"uuid", "odbc32", "odbccp32" ]
#libs += l
cfg['compiler_args'] = ['/std:c++latest']
cfg['include_dirs'] = [
	pybind11.get_include(),
	pybind11.get_include(True),
	r"C:\_Hassan\_Dev\_Lib\pybind11\include",
	r"C:\Users\hseed\Anaconda3\pkgs\python-3.7.3-h8c8aaf0_1\include",
	r"C:\_Hassan\_Dev\_HLibrary",
	r"C:\Program Files\OpenCV\include",
	r"C:\_Hassan\_Dev\_Lib\vcpkg\installed\x64-windows\include"
 ]

cfg['linker_args'] = [
#	'/LIBPATH:'+cvl,
#	'/LIBPATH:'+qtl,
#	r'/LIBPATH:C:\Users\hseed\Anaconda3\libs'
 ]

#cfg['libraries'] = libs

%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include <iostream>
#include <string>
#include <windows.h>
#include "hcv2.h"
#include "hcv.h"

#define LIB_NAME hcv

#define PASTE(a, b) a##b
/// pyptr(d,double)  --> auto d_ = d.request(); double* d__ = (double*) d_.ptr;;
#define pyptr(X,T) auto PASTE(X,_) = X.request(); PASTE(T, *) PASTE(X,__) = (PASTE(T, *)) PASTE(X,_).ptr
/// pyarr(pp,int,9) -->  py::array_t<int> pp = py::array_t<int>(9);
#define pyarr(X,T,sz)	py::array_t<T> X = py::array_t<T>(sz);

//#define pymat(in, out) \
//	pyptr(in,uint8_t);cv::Mat out; \
//	switch (PASTE(in,_).shape[2]) \
//	{ case 1: out = cv::Mat(PASTE(in,_).shape[0], PASTE(in,_).shape[1], CV_8UC1, PASTE(in,__));break; \
//	case 3: 	out = cv::Mat(PASTE(in,_).shape[0], PASTE(in,_).shape[1], CV_8UC3, PASTE(in,__));break;}

namespace py = pybind11;
using namespace pybind11::literals;
///-----------------------------------------------------------------------------------------
///-----------------------------------------------------------------------------------------
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
	switch (fdwReason)
	{case DLL_PROCESS_ATTACH:
	    {std::cout << "hello from hcv.\n";break;}
	case DLL_PROCESS_DETACH:
	    {std::cout << "bye from hcv.\n";break;}
	}
	return true;
}

///-----------------------------------------------------------------------------------------
void print(const std::string& name) {
	pr2("name = ",name);
}
///-----------------------------------------------------------------------------------------
void dtest(py::dict d)
{
	d["999"] = 000;
}
///-----------------------------------------------------------------------------------------
//void test(py::list l) {
//	//l.attr("pop")();
//	std::cout << "List has length " << l.size() << std::endl;
//
//	for (py::handle obj : l) {  // iterators!
//		std::cout << "  - " << obj.attr("__str__")().cast<std::string>() << std::endl;
//	}
//	l.append("ppp");  // automatic casting (through templating)!
//	l.append(l[3] + l[2]);  // automatic casting (through templating)!
//	l[0] = "999";
//}
///-----------------------------------------------------------------------------------------
void float_cast(py::list l) {
	l[0] = 0.0;
}
///-----------------------------------------------------------------------------------------
int add(py::int_& i, py::int_& j) {
	i = 99;
	return (int)i + (int)j;
}
///-----------------------------------------------------------------------------------------
py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2)
{
	auto buf1 = input1.request(), buf2 = input2.request();

	if (buf1.size != buf2.size)
		throw std::runtime_error("Input shapes must match");

	/*  allocate the buffer */
	//py::array_t<double> result = py::array_t<double>(buf1.size);
	pyarr(result,double,buf1.size);

	auto buf3 = result.request();

	double* ptr1 = (double*)buf1.ptr,
		* ptr2 = (double*)buf2.ptr,
		* ptr3 = (double*)buf3.ptr;

	int X = buf1.shape[0];
	int Y = buf1.shape[1];

	for (size_t idx = 0; idx < X; idx++)
		for (size_t idy = 0; idy < Y; idy++) {
			ptr3[idx * Y + idy] = ptr1[idx * Y + idy] + ptr2[idx * Y + idy];

			ptr1[idx * Y + idy] = 1234;
		}
	// reshape array to match input shape
	result.resize({ X,Y });
	return result;
}
///-----------------------------------------------------------------------------------------
//void cvshow(py::array_t<uint8_t> img)
//{
//	pymat(img,im);
//	QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
//	hcv::showImg("im2", im, 0);
//}

///==========================================================================================
///==========================================================================================
///==========================================================================================
//int floodFill(
//    int x, int y,
//    pyImg i_img_source,
//    pyImg i_mask,
//    int tolerance)
//{
//    np<uint8_t> img_source = np<uint8_t>(i_img_source);
//	np<uint8_t> mask = np<uint8_t>(i_mask);
//	int pxl;
//	Point p;
//	p.x = x; p.y = y;
//	Point pr;
//    std::vector<Point> pointVector;
//	pointVector.reserve(1000000);
//	pointVector.clear();
//	pointVector.push_back(p);
//	size_t flag = pointVector.size();
//	int pixel_count = 0;
//	while (flag)
//	{
//		p = pointVector[flag - 1];
//		pointVector.pop_back();
//		if (img_source.in(p))
//		{
//			img_source.set(255,p);
//			if (mask.get(p) < 10)
//			{
//				pixel_count++;
//				mask.set(255, p);
//			}
//			////----------------------------------------------------------------------------------------------
//			pr.x = p.x - 1;					pr.y = p.y;
//			if (img_source.in(pr))
//			{
//				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
//			}
//			////----------------------------------------------------------------------------------------------
//			pr.x = p.x + 1;					pr.y = p.y;
//			if (img_source.in(pr))
//			{
//				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
//			}
//			//------------------------------------------------------------------------------------------------
//			pr.x = p.x;					pr.y = p.y - 1;
//			if (img_source.in(pr))
//			{
//				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
//			}
//			////----------------------------------------------------------------------------------------------
//			pr.x = p.x;					pr.y = p.y + 1;
//			if (img_source.in(pr))
//			{
//				pxl = img_source.get(pr); if (pxl < tolerance) pointVector.push_back(pr);
//			}
//			//------------------------------------------------------------------------------------------------
//		}
//		flag = pointVector.size();
//	}
//	return pixel_count;
//}

///==========================================================================================
///==========================================================================================
///==========================================================================================
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<Point>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MODULE(LIB_NAME, m) {
	m.doc() = "hcv plugin"; // optional module docstring
	m.def("add", &add, "A function which adds two numbers", py::arg("i"), py::arg("j"));
	m.def("float_cast", &float_cast);

	m.def("print", &print);
	m.def("floodFill", &floodFill, "Flood Fill");



	m.def("dist", py::overload_cast<Point&, Point&>(&dist));
	m.def("dist", py::overload_cast<double&, double&, double&, double&>(&dist));
	m.def("angle", py::overload_cast<Point&, Point&>(&angle));
	m.def("getCenter", &getCenter, "Point getCenter(vPoint& pts)");
	m.def("getAverageRadius", &getAverageRadius, "double getAverageRadius(vPoint& pts, Point ctr)");
	m.def("getMaxRadius", &getMaxRadius, "double getMaxRadius(vPoint& pts, Point& ctr)");
	m.def("getMinRadius", &getMinRadius, "double getMinRadius(vPoint& pts, Point& ctr)");
	m.def("getRadii", &getRadii, "void getRadii(vPoint& pts, vDouble &radii, Point& ctr)");
	m.def("diffRadii", &diffRadii, "void diffRadii(vDouble& radii, vDouble& diff)");
	m.def("diffPolygon", &diffPolygon, "void diffPolygon(vPoint& pts, vDouble& diffX, vDouble& diffY)");
	m.def("getArea", &getArea, "double getArea(vPoint& pts)");
	m.def("getPerimeter", &getPerimeter, "double getPerimeter(vPoint& pts)");
	m.def("featherPoints", &featherPoints, "void featherPoints(vPoint& pts, Point& center, int border)");

	m.def("getLength",
	py::overload_cast<vPoint&, bool> (&getLength),
	"points"_a , "close"_a = false,
	"double getLength(vPoint& points, bool close = false)");

	m.def("interpolatePolygon",
	py::overload_cast<vPoint&, int, vPoint&, bool> (&interpolatePolygon),
	"source"_a , "result_count"_a, "result"_a, "close"_a = true,
    "void interpolatePolygon(vPoint& source, int result_count, vPoint &result, bool close = true)");

	m.def("interpolatePolygon",
	py::overload_cast<vPoint&, int, bool> (&interpolatePolygon),
	"source"_a , "result_count"_a, "close"_a = true,
    "void interpolatePolygon(vPoint& source, int result_count, bool close = true)");

    m.def("pointInPolygon", &pointInPolygon, "bool pointInPolygon(Point& point, vPoint& polygon)");

	m.def("roundPolygon",
	py::overload_cast<vPoint&, Point&, vDouble&, double> (&roundPolygon),
	"pts"_a , "center"_a, "radii_output"_a, "median"_a = 0.9,
    "double roundPolygon(vPoint& pts, Point& center, vDouble& radii_output, double median = 0.9)");

	m.def("fillPolygon", &fillPolygon,"void fillPolygon(Img &img, vPoint& pts)");
	m.def("printPointVector", &printPointVector,"void printPointVector(vPoint& pts)");

	m.def("add_arrays", &add_arrays, "Add two NumPy arrays");







	py::class_<Point>(m, "Point")
		.def(py::init<const int, const int>())
		.def_readwrite("x", &Point::x)
		.def_readwrite("y", &Point::y);

	py::class_<Ray>(m, "Ray")
		.def(py::init<>())

        .def_readwrite("img", &Ray::img)
        .def_readwrite("img_tmp", &Ray::img_tmp)
        .def_readwrite("mask", &Ray::mask)
        .def_readwrite("mask_tmp", &Ray::mask_tmp)
        .def_readwrite("grain_mask", &Ray::grain_mask)
        .def_readwrite("grain_pts", &Ray::grain_pts)

        .def_readwrite("im", &Ray::im)
        .def_readwrite("jm", &Ray::jm)

        .def_readwrite("ray_inc", &Ray::ray_inc)
        .def_readwrite("angle_inc", &Ray::angle_inc)
        .def_readwrite("tolerance", &Ray::tolerance)
        .def_readwrite("fill_tol", &Ray::fill_tol)
        .def_readwrite("min_grain_size", &Ray::min_grain_size)
        .def_readwrite("max_grain_size", &Ray::max_grain_size)
        .def_readwrite("feather_grain_pixels", &Ray::feather_grain_pixels)
        .def_readwrite("moving_average_distance", &Ray::moving_average_distance)
        .def_readwrite("interpolate_polygon_count", &Ray::interpolate_polygon_count)
        .def_readwrite("back_span", &Ray::back_span)
        .def_readwrite("ray_end_point", &Ray::ray_end_point)
        .def_readwrite("grain_center", &Ray::grain_center)
        .def_readwrite("grain_pts", &Ray::grain_pts)
        .def_readwrite("grain_radii", &Ray::grain_radii)
        .def_readwrite("TOL_START_RANGE", &Ray::TOL_START_RANGE)
        .def_readwrite("TOL_STEP", &Ray::TOL_STEP)
		.def("info", &Ray::info)
		.def("init", &Ray::init)
		.def("getGrain", &Ray::getGrain,"getGrain(Point p)")
		.def("test", &Ray::test);

	py::class_<stats>(m, "stats")
        .def_readwrite("mean", &stats::mean)
        .def_readwrite("min", &stats::min)
        .def_readwrite("max", &stats::max)
        .def_readwrite("stdev", &stats::stdev)
        .def_readwrite("median", &stats::median)
        .def_readwrite("stdev_ratio", &stats::stdev_ratio)
        .def_readwrite("area", &stats::area)
        .def_readwrite("perimeter", &stats::perimeter);

	py::class_<Img>(m, "Img")
		.def(py::init<pyImg&>())
		.def("init",  &Img::init)
		.def("add",   &Img::add)
		.def("copy",  py::overload_cast<Img&>(&Img::copy))
		.def("mult", py::overload_cast<int>(&Img::mult))
		.def("mult", py::overload_cast<Img&>(&Img::mult))
		.def("zero",  &Img::zero)
        .def("get", py::overload_cast<int>(&Img::get))
        .def("get", py::overload_cast<int,int>(&Img::get))
        .def("get", py::overload_cast<int,int,int>(&Img::get))
        .def("get", py::overload_cast<Point>(&Img::get))
        .def("get", py::overload_cast<Point,int>(&Img::get))
        .def("set", py::overload_cast<uint8_t,int>(&Img::set))
        .def("set", py::overload_cast<uint8_t,int,int>(&Img::set))
        .def("set", py::overload_cast<uint8_t,int,int,int>(&Img::set))
        .def("set", py::overload_cast<uint8_t, Point>(&Img::set))
        .def("set", py::overload_cast<uint8_t, Point,int>(&Img::set))
		.def("info",  &Img::info);


//    .def("set", py::overload_cast<int>(&Pet::set), "Set the pet's age")
//    .def("set", py::overload_cast<const std::string &>(&Pet::set), "Set the pet's name");
//clsExampleOne.def(py::init<>());
//clsExampleOne.def(py::init<std::string const&, ExampleOne::State>(), "fileName"_a, "state"_a=ExampleOne::State::RED);
//clsExampleOne.def(py::init<ExampleOne const&, bool>(), "other"_a, "deep"_a=true); // Copy constructor


    m.def("isInside", &isInside, "( Point p, vPoint vp) return if point inside polygon");
    m.def("angle", &angle, "( Point p, vPoint vp)");

///---------------------------------------------------------------------------
    m.def("vAsString", &vAsString, "std::string vAsString(vDouble &v, int digits = 16, char separator = 0)");
///---------------------------------------------------------------------------
    m.def("vSet", &vSet, "void vSet(vDouble &v, std::initializer_list<double> values)");
    m.def("vMin", &vMin, "double vMin(vDouble &v) ");
    m.def("vMax", &vMax, "double vMax(vDouble v)");
    m.def("vAverage", &vAverage, "double vAverage(vDouble &v)");
    m.def("vNorm", &vNorm, "double vNorm(vDouble &v)");
    m.def("vAdd", &vAdd, "result = a X + b Y void vAdd(vDouble &result, vDouble &X, vDouble &Y, double a = 1.0, double b = 1.0)");
    m.def("vMult", &vMult, "return = a X * Y  double vMult(vDouble &X, vDouble &Y,	double a = 1.0)");
    m.def("vThreshold", &vThreshold, "void vThreshold(vDouble &in, vDouble &out, double i_threshold)");
    m.def("vStats", &vStats, "stats vStats(vDouble &v)");

    m.def("vInterp", py::overload_cast<vDouble &,vDouble &,double,bool>(&vInterp),
     "return = yData(x) double vInterp(vDouble &xData, vDouble &yData, double x, bool extrapolate = false)");

    m.def("vInterp", py::overload_cast<vDouble &,vDouble &,vDouble &,vDouble &,bool>(&vInterp),
     "y[] = yData( x[] ) void vInterp(	vDouble &xData,	vDouble &yData,	vDouble &x,	vDouble &y,	bool extrapolate = false)");


	py::bind_vector<std::vector<int>>    (m, "vInt");
	py::bind_vector<std::vector<double>> (m, "vDouble");
	py::bind_vector<std::vector<Point>>  (m, "vPoint");
}