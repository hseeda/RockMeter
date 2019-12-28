from timeit import default_timer as timer
import time
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
    print("Time            = [ %.3f \u03BCs ]" % t)


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


def drawLine(inImg, p1, p2, i_line_width, color):
    cv2.line(inImg, p1, p2, color, i_line_width, cv2.LINE_AA)


def drawPolygon(inImg, points, i_line_width=1, closed=True,
                R=0, G=255, B=0):
    c = len(points) - 1
    if c > 0:
        for i in range(0, c):
            drawLine(inImg,
                     (points[i].x, points[i].y),
                     (points[i + 1].x, points[i + 1].y),
                     i_line_width, (R, G, B))
        if closed:
            drawLine(inImg,
                     (points[c].x, points[c].y),
                     (points[0].x, points[0].y),
                     i_line_width, (R, G, B))


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


def scaleImg(img, scale):
    width = int(img.shape[1] * scale)
    height = int(img.shape[0] * scale)
    return cv2.resize(img, (width, height), interpolation=cv2.INTER_LANCZOS4)

#
# class Ray:
#     def __init__(self):
#         self.data = []
#         self.ray = np.zeros(1000, dtype=np.int32)
#         self.grain_pts = np.zeros((9999, 2), dtype=np.int32)
#         self.pointVector = np.zeros((999999, 2), dtype=np.int32)
#         self.grain_radii = np.zeros(1000, dtype=np.double)
#         self.ray_end_point = [0, 0]
#         self.grain_center = [0, 0]
#         self.im = 0
#         self.jm = 0
#         self.pixel_count = 0
#         self.ray_inc = 1
#         self.angle_inc = 5
#         self.back_span = 0
#         self.tolerance = 25
#
#     def loadImg(self, file_name):
#         self.img = scaleImg(cv2.imread(file_name, cv2.IMREAD_GRAYSCALE), 0.5)
#         self.img_tmp = self.img.copy()
#         self.mask = self.img[:, :].copy() * 0
#         self.mask_tmp = self.mask.copy()
#         self.grain_mask = self.mask.copy()
#         self.mRay = hcv.Ray(self.img, self.img_tmp, self.mask, self.mask_tmp, self.grain_mask)
#
#     def setImg(self, i_img):
#         self.img = i_img.copy()
#         self.img_tmp = self.img.copy()
#         self.mask = self.img[:, :].copy() * 0
#         self.mask_tmp = self.mask.copy()
#         self.grain_mask = self.mask.copy()
#         self.mRay = hcv.Ray(self.img, self.img_tmp, self.mask, self.mask_tmp, self.grain_mask)
#

# ======================================================================================
def Ray_loadImg(Ray, file_name):
    img = scaleImg(cv2.imread(file_name, cv2.IMREAD_GRAYSCALE), 0.5)
    img_tmp = img.copy()
    mask = img[:, :].copy() * 0
    mask_tmp = mask.copy()
    grain_mask = mask.copy()
    Ray.init(img, img_tmp, mask, mask_tmp, grain_mask)


# ======================================================================================
def Ray_setImg(Ray, i_img):
    img = i_img.copy()
    img_tmp = img.copy()
    mask = img[:, :].copy() * 0
    mask_tmp = mask.copy()
    grain_mask = mask.copy()
    Ray.init(img, img_tmp, mask, mask_tmp, grain_mask)
