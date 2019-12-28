#include <iostream>


//#include <QtWidgets/QApplication>
namespace py = pybind11;


///------------------------------------------------------------------------------------
///------   Image functions
///------------------------------------------------------------------------------------

//static void cropImg(cv::Mat& inImg, cv::Mat& outImg, int w1, int w2, int h1, int h2)
//{
//    int ww = (w1 < w2) ? w1 : w2;
//    int hh = (h1 < h2) ? h1 : h2;
//    int w_ = abs(w1 - w2);
//    int h_ = abs(h1 - h2);
//    if (ww < 0)
//    {
//        w_ = (w1 > w2) ? w1 : w2;
//        ww = 0;
//    }
//    if (hh < 0)
//    {
//        h_ = (h1 > h2) ? h1 : h2;
//        hh = 0;
//    }
//    if (ww + w_ > inImg.size().width) { w_ = inImg.size().width - ww; }
//    if (hh + h_ > inImg.size().height) { h_ = inImg.size().height - hh; }
//    outImg = inImg(cv::Rect(ww, hh, w_, h_)).clone();
//}
//
//////-------------------------------------------
//static void scaleImg(cv::Mat& inImg, cv::Mat& outImg, double scale)
//{
//    resize(inImg, outImg, cv::Size(), scale, scale, cv::INTER_LANCZOS4);
//}
//
//////-------------------------------------------
//static void scaleImg(cv::Mat& inImg, cv::Mat& outImg, double scaleW, double scaleH)
//{
//    resize(inImg, outImg, cv::Size(), scaleW, scaleH, cv::INTER_LANCZOS4);
//}
//
//////-------------------------------------------
//static void scaleImgDim(cv::Mat& inImg, cv::Mat& outImg, int i_width, int i_height)
//{
//    double scale;
//    if (i_height == 0)
//    {
//        scale = double(i_width) / double(inImg.size().width);
//        scaleImg(inImg, outImg, scale);
//        return;
//    }
//    if (i_width == 0)
//    {
//        scale = double(i_height) / double(inImg.size().height);
//        scaleImg(inImg, outImg, scale);
//    }
//    else
//    {
//        double scaleW = double(i_width) / double(inImg.size().width);
//        double scaleH = double(i_height) / double(inImg.size().height);
//        scaleImg(inImg, outImg, scaleW, scaleH);
//    }
//}
//
//////-------------------------------------------
//static void scaleImgFit(cv::Mat& inImg, int i_width, int i_height, cv::Mat* outImg = nullptr)
//{
//    if (outImg == nullptr) outImg = &inImg;
//
//    const double scale1 = double(i_width) / double(inImg.size().width);
//    const double scale2 = double(i_height) / double(inImg.size().height);
//    if (scale1 > scale2)
//    {
//        i_width = int(round(inImg.size().width * scale2));
//    }
//    else
//    {
//        i_height = int(round(inImg.size().height * scale1));
//    }
//    resize(inImg, *outImg, cv::Size(i_width, i_height), 0, 0, cv::INTER_LANCZOS4);
//}

//////-------------------------------------------
//static void blurImg(cv::Mat& inImg, cv::Mat& outImg, int isz = 5)
//{
//    // GaussianBlur (   InputArray src, OutputArray dst, Size ksize,
//    //                  double sigmaX, double sigmaY=0,
//    //                  int borderType=BORDER_DEFAULT )
//    GaussianBlur(inImg, outImg, cv::Size(isz, isz), 0.0, 0.0);
//}
//
//////-------------------------------------------
//static void CannyB(
//    cv::Mat& image,
//    cv::Mat& edges,
//    int gaussian_blur_window_size,
//    double threshold1,
//    double threshold2,
//    int apertureSize,
//    bool L2gradient)
//{
//    cv::Mat tmp = image.clone();
//    blurImg(image, tmp, gaussian_blur_window_size);
//    Canny(tmp, edges, threshold1, threshold2, apertureSize, L2gradient);
//}


////-------------------------------------------
/// MASKS
//static cv::Mat createMask(cv::Mat& inImg, int border = 0)
//{
//    int w = inImg.size().width + 2 * border;
//    int h = inImg.size().height + 2 * border;
//    cv::Mat mask(cv::Size(w, h), CV_8UC1);
//    return mask;
//}
//
//
//////-------------------------------------------
//static void invMask(cv::Mat& inMask, cv::Mat* outMask = nullptr)
//{
//    if (outMask == nullptr)
//    {
//        inMask = 255 - inMask;
//        //bitwise_not(inMask, inMask);
//    }
//    else
//    {
//        *outMask = inMask.clone();
//        *outMask = 255 - *outMask;
//        //bitwise_not(inMask, *outMask);
//    }
//}
//
//static void addMask(cv::Mat& inImg, cv::Mat& mask, int R = 255, int G = 255, int B = 0)
//{
//    std::vector<cv::Mat> rgbChannels(3);
//    split(inImg, rgbChannels);
//    rgbChannels[0] += mask * B / 255;
//    rgbChannels[1] += mask * G / 255;
//    rgbChannels[2] += mask * R / 255;
//    merge(rgbChannels, inImg);
//}
//
//static int threshold(cv::Mat& inImg, cv::Mat& outImg, int lower_than)
//{
//    int cnt = 0;
//    int w = inImg.size().width;
//    int h = inImg.size().height;
//
//    outImg = inImg.clone() * 0;
//
//    Point p;
//    int v;
//
//    for (auto i = 0; i < w; i++)
//    {
//        p.x = i;
//        for (auto j = 0; j < h; j++)
//        {
//            p.y = j;
//            v = inImg.at<uchar>(p);
//
//            if (v < lower_than)
//            {
//                outImg.at<uchar>(p) = uchar(0);
//                cnt++;
//            }
//            else outImg.at<uchar>(p) = uchar(v);
//        }
//    }
//    return cnt;
//}

///------------------------------------------------------------------------------------
///------   Draw functions
///------------------------------------------------------------------------------------

//static void drawLine(cv::Mat& inImg, int x1, int y1, int x2, int y2,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    line(inImg, Point(x1, y1), Point(x2, y2), cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
//}

////-------------------------------------------
//static void drawLine(cv::Mat& inImg, Point& p1, Point& p2,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    line(inImg, p1, p2, cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
//}

////-------------------------------------------
//static void drawLineNAA(cv::Mat& inImg, Point& p1, Point& p2,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    line(inImg, p1, p2, cv::Scalar(B, G, R), i_line_width);
//}

////-------------------------------------------
//static void drawCircle(cv::Mat& inImg, int x, int y, int radius,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    circle(inImg, Point(x, y), radius, cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
//}

////-------------------------------------------
//static void drawCircle(cv::Mat& inImg, Point& center, int radius,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    circle(inImg, center, radius, cv::Scalar(B, G, R), i_line_width, cv::LINE_AA);
//}

////-------------------------------------------
//static void drawPoint(cv::Mat& inImg, Point& point,
//    int R = 255, int G = 255, int B = 0)
//{
//    inImg.at<cv::Vec3b>(point) = cv::Vec3b(B, G, R);
//    //inImg.at<uchar>(r.y, r.x, 0) = B;
//    //inImg.at<uchar>(r.y, r.x, 1) = B;
//    //inImg.at<uchar>(r.x, r.y, 2) = R;
//}

////-------------------------------------------
//static void drawLines(cv::Mat& inImg, Point& center, std::vector<Point>& points,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    size_t c = points.size();
//    for (size_t i = 0; i < c; i++)
//    {
//        drawLine(inImg, center, points[i], i_line_width, R, G, B);
//    }
//}

////-------------------------------------------
//static void drawPointsC(cv::Mat& inImg, std::vector<Point>& points,
//    int radius = 1,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 255)
//{
//    size_t c = points.size();
//    for (size_t i = 0; i < c; i++)
//        drawCircle(inImg, points[i], radius, i_line_width, R, G, B);
//}

////-------------------------------------------
//static void drawPoints(cv::Mat& inImg, std::vector<Point>& points,
//    int R = 255, int G = 255, int B = 0)
//{
//    size_t c = points.size();
//    for (size_t i = 0; i < c; i++)
//        inImg.at<cv::Vec3b>(points[i]) = cv::Vec3b(B, G, R);
//}

////-------------------------------------------
//static void drawPolygon(cv::Mat& inImg, std::vector<Point>& points,
//    bool closed = true,
//    int i_line_width = 1, int R = 0, int G = 255, int B = 0)
//{
//    size_t c = points.size();
//    if (c > 0)
//    {
//        c = c - 1;
//        for (size_t i = 0; i < c; i++)
//            drawLine(inImg, points[i], points[i + 1], i_line_width, R, G, B);
//        if (closed) drawLine(inImg, points[c], points[0], i_line_width, R, G, B);
//    }
//}

////-------------------------------------------
//static void drawPolygonNAA(cv::Mat& inImg, std::vector<Point>& points,
//    bool closed = true,
//    int i_line_width = 1, int R = 0, int G = 255, int B = 0)
//{
//    size_t c = points.size();
//    if (c > 0)
//    {
//        c = c - 1;
//        for (size_t i = 0; i < c; i++)
//            drawLineNAA(inImg, points[i], points[i + 1], i_line_width, R, G, B);
//        if (closed) drawLine(inImg, points[c], points[0], i_line_width, R, G, B);
//    }
//}

////-------------------------------------------
//static void drawText(cv::Mat& inImg, std::string txt,
//    int x, int y, double font_size,
//    int i_line_width = 1, int R = 0, int G = 0, int B = 0)
//{
//    putText(
//        inImg, //target image
//        txt.c_str(), //text
//        Point(x, y), //top-left position
//        cv::FONT_HERSHEY_DUPLEX,
//        font_size,
//        CV_RGB(R, G, B), //font color
//        i_line_width,
//        cv::LINE_AA);
//}

//static void fillPolygon(cv::Mat& inImg, std::vector<Point>& pts,
//    int R = 255, int G = 255, int B = 255)
//{
//    cv::fillConvexPoly(inImg, pts.data(), int(pts.size()), cv::Scalar(B, G, R));
//}

///====================================================================================
///------------------------------------------------------------------------------------
///------   Ray class   ---------------------------------------------------------------
///------------------------------------------------------------------------------------
///====================================================================================
class Ray
{
public:
    vInt ray;
    int im, jm;
    int pixel_count;
    int ray_inc;
    int angle_inc;
    int tolerance;
    int min_grain_size;
    int max_grain_size;
    int feather_grain_pixels;
    int moving_average_distance;
    int interpolate_polygon_count;

    int TOL_START_RANGE;
    int TOL_STEP;

    double fill_tol;

    // no. of points to look back at
    int back_span;

    Point ray_end_point;
    Point grain_center;
    stats grain_stats;
    std::vector<stats> grains;

    vPoint grain_pts;
    vDouble grain_radii;

    np<uint8_t> np_img;
    np<uint8_t> np_img_tmp;
    np<uint8_t> np_mask;
    np<uint8_t> np_mask_tmp;
    np<uint8_t> np_grain_mask;

    pyImg img, img_tmp, mask, mask_tmp, grain_mask;

    Ray()
    {
    }

    void init(pyImg &i_img, pyImg &i_img_tmp, pyImg &i_mask, pyImg &i_mask_tmp, pyImg &i_grain_mask)
    {
        img = i_img;
        img_tmp = i_img_tmp;
        mask = i_mask;
        mask_tmp = i_mask_tmp;
        grain_mask = i_grain_mask;

        fill_tol = 0.3;
        min_grain_size = 200;
        max_grain_size = 9999999;
        feather_grain_pixels = 3;
        moving_average_distance = 4;
        interpolate_polygon_count = 64;

        ray.reserve(500);
        tolerance = 25;
        back_span = 0;
        ray_inc = 1;
        angle_inc = 1;

        np_img.init(img);
        np_img_tmp.init(img_tmp);
        np_mask.init(mask);
        np_mask_tmp.init(mask_tmp);
        np_grain_mask.init(grain_mask);

        im = (int) np_img. buffer.shape[0];
        jm = (int) np_img. buffer.shape[1];

        TOL_START_RANGE =  10;
        TOL_STEP  = 5;

    }
///------------------------------------------------------------------------------
//    Ray(pyImg &img, pyImg &img_tmp, pyImg &mask, pyImg &mask_tmp, pyImg &grain_mask)
//    {
//        init(img, img_tmp, mask, mask_tmp, grain_mask);
//    }
///------------------------------------------------------------------------------
    void info(){
         np_img.info(); np_img_tmp.info();
         np_mask.info(); np_grain_mask.info();
    }
///------------------------------------------------------------------------------
    ~Ray() = default;
///------------------------------------------------------------------------------
    void addGrain()
    {
        grain_stats.area = getArea(grain_pts);
        grain_stats.perimeter = getPerimeter(grain_pts);
        grains.push_back(grain_stats);
    }
///------------------------------------------------------------------------------
    void printGrains()
    {
        for (stats s : grains)
        {
            pr6c(s.area,s.perimeter,s.mean,s.min,s.max,s.stdev_ratio);
        }
    }
///------------------------------------------------------------------------------
    bool checkPointInMask(Point& p)
    {
        return np_mask.get(p) != 255;
    }
///------------------------------------------------------------------------------
    bool pointInsideImage(Point& p) const
    {
        if (p.x <= 0) return false;
        if (p.x > im) return false;
        if (p.y <= 0) return false;
        if (p.y > jm) return false;
        return true;
    }
///------------------------------------------------------------------------------
    void featherGrainMask(vPoint &pts, Point& center, int border = 3)
    {
        featherPoints(pts, center, border);
        np_grain_mask.zero();
        fillPolygon(np_grain_mask, pts);
    }
///------------------------------------------------------------------------------
    bool nfill3(Point p, double fill_threshold, int Min_GRAIN_SIZE = 4, int Max_GRAIN_SIZE = 400)
    {
        int opc = 0;
        int cnt = 0;
        np_mask_tmp.copy(np_mask); np_img_tmp.copy(np_img); np_img_tmp.add(np_mask);
        pixel_count = floodFill(p, np_img_tmp, np_mask_tmp, TOL_START_RANGE);
        if (pixel_count <= 4) return false;
        for (int i = TOL_START_RANGE + TOL_STEP; i < 255; i = i + TOL_STEP)
        {
            opc = pixel_count;
            tolerance = i;
            np_mask_tmp.copy(np_mask); np_img_tmp.copy(np_img); np_img_tmp.add(np_mask);
            pixel_count = floodFill(p, np_img_tmp, np_mask_tmp, tolerance);
            double ratio = double(pixel_count - opc) / double(opc);

            if (ratio > fill_threshold || cnt > 20 /*no of fills with relatively no change*/)
            {
                if (pixel_count < Min_GRAIN_SIZE || pixel_count > Max_GRAIN_SIZE)
                {np_grain_mask.zero(); return false;}
                tolerance = tolerance - TOL_STEP;
                if(tolerance < TOL_START_RANGE) tolerance = TOL_START_RANGE;
                np_img_tmp.copy(np_img);
                np_img_tmp.add(np_mask);
                np_grain_mask.zero();
                floodFill(p, np_img_tmp, np_grain_mask, tolerance);
                return true;
            }
            else
                cnt++;
        }
        return false;
    }
///------------------------------------------------------------------------------
    bool getRay(Point &center, double angle, Img &inImg)
    {
        double c = 0;
        int cc = -1;
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
                ray_end_point = Point(i, j);
                return true;
            }
            int pxl = inImg.get(j, i);
            if (pxl > tolerance)
            {
                ray_end_point = Point(i, j);
                return true;
            }
            //}
        }
    }
///------------------------------------------------------------------------------
    bool getRayInv(Point &center, double angle, Img &inImg)
    {
        double c = 0;
        int cc = -1;
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
                ray_end_point = Point(i, j);
                return true;
            }
            int pxl = inImg.get(j, i);
            if (pxl <= 0)
            {
                ray_end_point = Point(i, j);
                return true;
            }
        }
    }
///------------------------------------------------------------------------------
    void getRaysInv(Point& center, vPoint& points, Img &inImg)
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
///------------------------------------------------------------------------------
    void getRays(Point center, vPoint& points, Img &inImg)
//    void getRays(Point center, vPoint& points, pyImg &inImg)
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
    bool getGrain(Point& p)
    {
        if (checkPointInMask(p)) {
            const bool fill_result = nfill3(p, fill_tol, min_grain_size, max_grain_size);
            if (fill_result)
            {
                getRaysInv(p, grain_pts, np_grain_mask);
                interpolatePolygon(grain_pts, interpolate_polygon_count, true);
                p = getCenter(grain_pts);
                if (moving_average_distance != 0)///---------- Moving Average
                {
                    grain_pts = movingAveragePointVector(grain_pts, moving_average_distance);
                    p = getCenter(grain_pts);
                }
                if (feather_grain_pixels != 0)///---------- Feather Grain
                {
                    featherPoints(grain_pts, p, feather_grain_pixels);
                }
                np_grain_mask.zero();
                fillPolygon(np_grain_mask, grain_pts);
                grain_center = p;
                getRadii(grain_pts, grain_radii, p);
                grain_stats = vStats(grain_radii);
            }
            return fill_result;
        }
        return false;
    }

    void test(Point& p, vPoint &points)
    {
//        mask_floodFill(p,np_mask);
//        pr3(p.x,",",p.y);
//        void getRays(Point center, vPoint& points, Img &inImg)
        getRays(p, points, np_mask);
    }

};
