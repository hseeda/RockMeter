#pragma once
#define RAND_MAX 0x7fff

static int hround(double v)
{
	return (v - int(v) < 0.5) ? int(v) : int(v) + 1;
}

static double dRand(double fMax=1.0, double fMin=0.0)
{
		double f = (double) rand() / RAND_MAX;
		return fMin + f * (fMax - fMin);
}
static int iRand(int fMax= 255, int fMin = 0)
{
		double f = (double) rand() / RAND_MAX;
		return fMin + int(f * double (fMax - fMin));

}
