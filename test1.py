int main()
{
	QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

	static const char* names[] = {
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/mix6_out.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/5.jpgout.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/7.jpgout.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/2.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/1.jpgout.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/gravels1_out.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/gravels1_out_f_1.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/Gravelometer/AH/gravels3_out.jpg",
		"C:/_Hassan/_Dev/_Projects/_HOpenCV/x64/Release/data/sand.jpg", nullptr };

	///---------------------------------------------------------
	static const char* fn = names[1];
	///---------------------------------------------------------

	cv::Mat im, im2, ic;

	std::vector<double> radii;

	im = cv::imread(fn, cv::IMREAD_COLOR);

	hcv::checkImg(im);
	hcv::scaleImg(im, im, 0.5);

	im2 = im.clone();
	cv::cvtColor(im, ic, cv::COLOR_BGR2GRAY);
	hcv::printImgSize(im, "im");
	hcv::printImgSize(ic, "ic");
	///----------------------------------------------------
	hcv::hcvRay r;
	///----------------------------------------------------
	r.setImage(&ic);
	r.angle_inc = 5;
	///------------------------------------------------------------ (nfill2 loop)
	if (1) {
		double fill_tol = .01;
		int min_grain_size;
		int max_grain_size;
		int grain_count = 0;

		bool fill_result;
		int feather_grain_pixels = 3;
		int moving_average_distance = 3;
		double grain_quality_min = 0.2; /// grain quality = stdev(radii)/ mean(radii)
		double grain_quality;
		int interpolate_polygon_count = 45;
		int w = im2.size().width;
		int h = im2.size().height;

		int cr = 0;

		int mgs = 10000;

		hcv::showImg("im2", im2, -1);

		max_grain_size = 10000;
		min_grain_size = 20;
		int stp = sqrt(min_grain_size) / 2;

		for (int ix = stp / 2; ix < w; ix += stp)
			for (int iy = stp / 2; iy < h; iy += stp)
			{
				cr++;
				std::cout << cr << "\r";

				cv::Point p = cv::Point(ix, iy);

				///------------------------------------------
				fill_result = r.getGrain(p,
					fill_tol, min_grain_size, max_grain_size,
					feather_grain_pixels,
					moving_average_distance,
					interpolate_polygon_count);
				///------------------------------------------
				if (fill_result)
				{
					if (r.grain_stats.stdev_ratio < grain_quality_min)
					{
						r.mergeGrainMask();
						r.addGrain();
						grain_count++;
						pr(grain_count);
						hcv::drawPolygon(im2, r.grain_pts, true, 1, 255, 255, 0);
						hcv::drawCircle(im2, r.grain_center, 2, 2, 255, 0, 0);
					}
					if (devisableBy(grain_count, 10))
						hcv::showImg("im2", im2, -1);
					//}
				}
			}
		//}
		r.printGrains();
		hcv::showImg("im2", im2, 1);
		//writeCSV("c:\\tmp\\1.csv", radii, 6);
		hpause;
		exit(0);
	}
}git 