import rock_meter_lib as rl
import math  
import hcv

r = hcv.Ray()
cimg = rl.Ray_loadImg(r, r"C:\_Hassan\_Dev\_Projects\_HOpenCV\x64\Release\data\Gravelometer\AH\gravels1_out.jpg")

# 	cv::Mat im, im2, ic;

# std::vector<double> radii;
radii = hcv.vDouble
rl.show(r.img)

r.angle_inc = 5
r.fill_tol    = 0.01
grain_count = 0

r.feather_grain_pixels = 3
r.moving_average_distance = 3

r.interpolate_polygon_count = 45
w = r.im
h = r.jm
cr = 0

mgs = 10000

# rl.show(cimg, 0, "im2")

r.max_grain_size = 10000
r.min_grain_size = 20
stp = int (math.sqrt(r.min_grain_size) / 2)

grain_quality_min = 0.2 # grain quality = stdev(radii)/ mean(radii)


                                

# for (int ix = stp / 2; ix < w; ix += stp)
# for (int iy = stp / 2; iy < h; iy += stp)
for ix in range(stp // 2, w, stp):
   for iy in range(stp // 2, h, stp):
        cr+= 1
        p = hcv.Point(ix,iy)
        # print(p, p.x, p.y)
# %%------------------------------------------
        fill_result = r.getGrain(p)
        if fill_result:
            print(p.x, p.y)
            print (fill_result, r.grain_stats.stdev_ratio)
            if (r.grain_stats.stdev_ratio < grain_quality_min):
                print("mask merged")
                rl.show(r.mask_tmp)
                r.mergeGrainMask()
                rl.show(r.mask)
                
                
                #                 
# 				{
# 					if (r.grain_stats.stdev_ratio < grain_quality_min)
# 					{
# 						r.mergeGrainMask();
# 						r.addGrain();
# 						grain_count++;
# 						pr(grain_count);
# 						hcv::drawPolygon(im2, r.grain_pts, true, 1, 255, 255, 0);
# 						hcv::drawCircle(im2, r.grain_center, 2, 2, 255, 0, 0);
# 					}
# 					if (devisableBy(grain_count, 10))
# 						hcv::showImg("im2", im2, -1);
# 					//}
# 				}





# 			}
# 		//}
# 		r.printGrains();
# 		hcv::showImg("im2", im2, 1);
# 		//writeCSV("c:\\tmp\\1.csv", radii, 6);
# 		hpause;
# 		exit(0);
# 	}
# }