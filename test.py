from timeit import default_timer as timer
from numba import jit
import numpy as np
import cv2
import hcv

# hello from vscode

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
    cv2.imshow("img", img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()


fn = "1.jpg"
i = cv2.imread(fn)
# hcv.cvi(i)
# show(i)

a = hcv.Img(i)
a.info()
print(hcv.dist_d(0,0,1,1))