from rock_meter_lib import *
import cv2
import hcv
# =============================================================================
if __name__ == "__main__":
    # fn = r"img\test.jpg"
    # simg = cv2.imread(fn)
    # i = createMask(simg)
    # i = cv2.rectangle(i, (100, 200), (400, 400), (255, 255, 255), 3)
    # j = i.copy()

    r = hcv.Ray()
    Ray_loadImg(r, r"C:\_Hassan\_Dev\_Projects\_HOpenCV\x64\Release\data\Gravelometer\AH\gravels1_out.jpg")

    p = hcv.Point(100, 600)
    pts = hcv.vPoint()

    show(r.img_tmp)

    print(r.getGrain(p))
    drawPolygon(r.img, r.grain_pts, 3, True, 255, 0, 0)
    # show(r.mask_tmp)
    show(r.img)
