from timeit import default_timer as timer
import math
from numba import jit
import cv2
import numpy as np
import hcv


def timeit(f, p, n):
    cc = timer()
    for i in range(0, n):
        f(*p)
    t = (- cc + timer()) * 1000
    print("Dial Reading    = [ %.3f \u03BCs ]" % t)


def timeits(f, p, n):
    cc = timer()
    for i in range(0, n):
        f(p)
    t = (- cc + timer()) * 1000
    print("Dial Reading    = [ %.3f \u03BCs ]" % t)


# %% get npArray attributes
def npa(in1):
    print('type   : ', type(in1))
    print('ndim   : ', in1.ndim)
    print('shape  : ', in1.shape)
    print('size   : ', in1.size)
    print('dtype  : ', in1.dtype)
    print('strides: ', in1.strides)
    print('-----------------------')


# %% get npArray Statistics
def nps(in1):
    print('min     : ', in1.min())
    print('max     : ', in1.max())
    print('mean    : ', in1.mean())
    print('std     : ', in1.std())
    print('variance: ', np.var(in1))
    print('m25%    : ', np.percentile(in1, 25))
    print('m50%    : ', np.percentile(in1, 50))
    print('m75%    : ', np.percentile(in1, 75))
    print('-----------------------')


def cv_version():
    print(cv2.__version__)


# %% show image using cv2
def show(img):
    cv2.imshow("img", img);
    cv2.waitKey(0)
    cv2.destroyAllWindows()


# %% crop image
def cropImage(img, x_start, x_end, y_start, y_end):
    width = img.shape[1]
    height = img.shape[0]
    xs = int(width * x_start)
    xe = int(width * x_end)
    ys = int(height * y_start)
    ye = int(height * y_end)
    return img[ys:ye, xs:xe, :]


def createMask(img, border=0):
    w = img.shape[0] + 2 * border
    h = img.shape[1] + 2 * border
    return np.zeros((w, h, 1), dtype=np.uint8)


def invMask(inMask):
    np.bitwise_not(inMask, inMask)


@jit
def pdist(a, b):
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


@jit
def angle(a, b):
    ang = math.atan2(b[0] - a[0], b[1] - a[1]) * 180.0 / math.pi
    if ang < 0.0:
        ang = 360.0 + ang
    return ang


@jit
def getCenter(pts):
    n = pts.shape[0]
    x = 0
    y = 0
    for i in range(0, n):
        x += pts[i, 0]
        y += pts[i, 1]
    return x / n, y / n


@jit
def getCenternp(pts):
    n = pts.shape[0]
    return np.sum(pts[:, 0]) / n, np.sum(pts[:, 1]) / n


@jit
def getAverageRadius(pts, ctr):
    c = pts.shape[0]
    # ctr = getCenternp(pts)
    radius = 0.0
    for i in range(0, c):
        radius += pdist(pts[i], ctr)
    return radius / c


@jit
def getMaxRadius(pts, ctr):
    c = pts.shape[0]
    radius = 0.0
    for i in range(0, c):
        t = pdist(pts[i], ctr)
        if t > radius:
            radius = t
    return radius


@jit
def getMinRadius(pts, ctr):
    c = pts.shape[0]
    radius = 99999999.0
    for i in range(0, c):
        t = pdist(pts[i], ctr)
        if t < radius:
            radius = t
    return radius


@jit
def getRadii(pts, ctr):
    c = pts.shape[0]
    radii = np.zeros(c, dtype=np.double)
    for i in range(0, c):
        radius = pdist(pts[i], ctr)
        radii[i] = radius
    return radii


@jit
def checkPointInMask(mask, p):
    t = mask[p[1], p[0], 0]
    if t == 255:
        return False
    else:
        return True


@jit
def pointInsideImage(img, p):
    if p[0] < 0:
        return False
    if p[0] >= img.shape[1]:
        return False
    if p[1] < 0:
        return False
    if p[1] >= img.shape[1]:
        return False
    return True


@jit
def floodFill(p, img_source, mask, tolerance):
    pointVector = np.zeros((9999999, 2), dtype=np.int32)
    pointer = 0
    pr = p.copy()
    pointVector[pointer, 0] = p[0]
    pointVector[pointer, 1] = p[1]
    flag = pointer + 1
    pixel_count = 0
    while flag:
        p[0] = pointVector[pointer, 0]
        p[1] = pointVector[pointer, 1]
        pointer = pointer - 1
        if pointInsideImage(img_source, p):
            img_source[p[0], p[1], 0] = 255

            if mask[p[0], p[1], 0] < 10:
                pixel_count += 1
                mask[p[0], p[1], 0] = 255

            pr[0] = p[0] - 1
            pr[1] = p[1]
            if pointInsideImage(img_source, pr):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointer = pointer + 1
                    pointVector[pointer, 0] = pr[0]
                    pointVector[pointer, 1] = pr[1]

            pr[0] = p[0] + 1
            pr[1] = p[1]
            if pointInsideImage(img_source, pr):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointer = pointer + 1
                    pointVector[pointer, 0] = pr[0]
                    pointVector[pointer, 1] = pr[1]

            pr[0] = p[0]
            pr[1] = p[1] - 1
            if pointInsideImage(img_source, pr):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointer = pointer + 1
                    pointVector[pointer, 0] = pr[0]
                    pointVector[pointer, 1] = pr[1]

            pr[0] = p[0]
            pr[1] = p[1] + 1
            if pointInsideImage(img_source, pr):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointer = pointer + 1
                    pointVector[pointer, 0] = pr[0]
                    pointVector[pointer, 1] = pr[1]

        flag = pointer + 1
    return pixel_count


@jit
def floodFill_list(p, img_source, mask, tolerance):
    pxl = 0
    pr = [0, 0]
    pointVector = [p.copy()]
    flag = len(pointVector)
    pixel_count = 0
    while flag:
        p = pointVector[flag - 1]
        pointVector.pop()
        if pointInsideImage(img_source, tuple(p)):
            img_source[p[0], p[1], 0] = 255
            if mask[p[0], p[1], 0] < 10:
                pixel_count += 1
                mask[p[0], p[1], 0] = 255
            pr[0] = p[0] - 1
            pr[1] = p[1]
            if pointInsideImage(img_source, tuple(pr)):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointVector.append(pr.copy())
            pr[0] = p[0] + 1
            pr[1] = p[1]
            if pointInsideImage(img_source, tuple(pr)):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointVector.append(pr.copy())
            pr[0] = p[0]
            pr[1] = p[1] + 1
            if pointInsideImage(img_source, tuple(pr)):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointVector.append(pr.copy())
            pr[0] = p[0]
            pr[1] = p[1] - 1
            if pointInsideImage(img_source, tuple(pr)):
                pxl = img_source[pr[0], pr[1], 0]
                if pxl < tolerance:
                    pointVector.append(pr.copy())
        flag = len(pointVector)
    return pixel_count


def fillPolygon(inImg, pts, R=255, G=255, B=255):
    cv2.fillConvexPoly(inImg, pts, int(pts.size()), (B, G, R))


def featherPoints(pts, center, border):
    c = pts.size()
    for i in range(0, c):
        ang = round(angle(center, pts[i]))
        radius = pdist(center, pts[i])
        r = int(round(radius)) + border
        pts[i, 0] = int(center.x + r * math.cos(ang))
        pts[i, 1] = int(center.y + r * math.sin(ang))


def featherGrainMask(grain_mask, pts, center, border=3):
    featherPoints(pts, center, border)
    grain_mask = grain_mask * 0
    fillPolygon(grain_mask, pts)


# def interpolatePolygon(source, result_count, close = True, result=None):
#     if result is None:
#         result = & tmp;
#         flag = true;
#     flag = False
#     tmp = [0,0]
#     p = [0,0]
#     if (close) :
#         source.append(source[0])
#     result->clear();
#
# if (source.size() < 2 | | result_count < 2)
# {
# // degenerate source vector or result_count value
# // for simplicity, this returns an empty result
# // but special cases may be handled when appropriate for the application
# return;
# }
#
# const
# double
# total_length = getLength(source, close);
# const
# double
# segment_length = total_length / (result_count - 1);
#
# // start and finish
# are
# the
# current
# source
# segment
# 's endpoints
# auto
# start = source.begin();
# auto
# finish = start + 1;
#
# double
# src_segment_offset = 0;
# double
# src_segment_length = dist(*start, *finish);
#
# // The
# first
# point in the
# result is the
# same as the
# first
# point
# // in the
# source.
#     result->push_back(*start);
# double
# dd;
#
# for (std::size_t i = 1; i < result_count - 1; ++i)
# {
# const
# double
# next_offset = segment_length * i;
#
# while (src_segment_offset + src_segment_length < next_offset)
#     {
#         src_segment_offset += src_segment_length;
#     start = finish + +;
#     src_segment_length = dist(*start, *finish);
#     }
#
#     const
#     double
#     part_offset = next_offset - src_segment_offset;
#     const
#     double
#     part_ratio = part_offset / src_segment_length;
#
#     dd = finish->x - start->x;
#     p.x = start->x + int(round(part_ratio * dd));
#     dd = finish->y - start->y;
#     p.y = start->y + int(round(part_ratio * dd));
#     result->push_back(p);
#     }
#     result->push_back(source.back());
#
#     if (close)
#     {
#     source.pop_back();
#     result->pop_back();
#     }
#
#     if (flag)
#     {
#     std::vector < cv::Point > out;
#     source = tmp;
#     }
#     }

class hcvRay:
    def __init__(self):
        self.data = []
        self.ray = np.zeros(1000, dtype=np.int32)
        self.grain_pts = np.zeros((9999, 2), dtype=np.int32)
        self.pointVector = np.zeros((999999, 2), dtype=np.int32)
        self.grain_radii = np.zeros(1000, dtype=np.double)
        self.img = np.zeros((1, 1, 1), dtype=np.uint8)
        self.img_tmp = np.zeros((1, 1, 1), dtype=np.uint8)
        self.mask = np.zeros((1, 1, 1), dtype=np.uint8)
        self.grain_mask = np.zeros((1, 1, 1), dtype=np.uint8)
        self.ray_end_point = [0, 0]
        self.grain_center = [0, 0]
        self.im = 0
        self.jm = 0
        self.pixel_count = 0
        self.ray_inc = 1
        self.angle_inc = 5
        self.back_span = 0
        self.tolerance = 25

    def setImage(self, inImg):
        self.img = inImg
        self.im = inImg.shape[0] - 1
        self.jm = inImg.shape[1] - 1
        self.mask = createMask(self.img, 0)
        self.grain_mask = createMask(self.img, 0)
        self.grain_mask = self.grain_mask * 0
        self.mask = self.mask * 0

    def getRay(self, center, angle, inImg):
        c = 0.0
        cc = -1
        di = math.cos(angle)
        dj = math.sin(angle)
        ic = float(center[0])
        jc = float(center[1])
        while 1:
            cc += 1
            c += self.ray_inc
            i = int(round(ic + c * di))
            j = int(round(jc + c * dj))
            if i < 0 or i >= self.im or j < 0 or j >= self.jm:
                if i < 0:
                    i = 0
                if i >= self.im:
                    i = self.im - 1
                if j < 0:
                    j = 0
                if j >= self.jm:
                    j = self.jm - 1
                self.ray_end_point = [i, j]
                return True
            pxl = inImg[j, i]
            if pxl > 0:
                self.ray_end_point = [i, j]
                return True




# =============================================================================
if __name__ == "__main__":
    cv_version()
    fn = r"img\test.jpg"
    simg = cv2.imread(fn)
    i = createMask(simg)
    i = cv2.rectangle(i, (100, 200), (400, 400), (255, 255, 255), 3)
    j = i.copy()
    # print(npa(i))
    show(i)

    # print(angle((11, 11), (9, 9)))
    # a = np.ones((10000, 2), dtype=np.uint8)
    # timeit(getRadii, (a, (0, 0)), 10000)
    # print(getMinRadius(a, (1, 3)))
    # tt = getRadii(a, (1, 5))
    self = np.zeros(2, dtype=np.int32)
    p = [600, 600]
    pp = np.zeros(2, dtype=np.int32)
    pp[0] = 300
    pp[1] = 300
    # timeit(floodFill, (pp, i, j, 10), 1)
    timeit(hcv.floodFill, (110, 210, i, j, 10), 1)
    show(j)

    fnames = [
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/mix6_out.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/5.jpgout.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/7.jpgout.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/2.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/1.jpgout.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/gravels1_out.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/gravels1_out_f_1.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/gravels3_out.jpg",
        "C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/sand.jpg"]

    print(hcv.mainn(fnames[8]))


    im = cv2.imread(fn, cv2.IMREAD_COLOR)
    npa(im)
    im = cv2.resize(im, (0,0), 0.5, 0.5, cv2.INTER_LANCZOS4)
    im2 = im.copy()
    ic = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    npa(im)
    npa(ic)