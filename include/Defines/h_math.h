#pragma once

///////////////////////////////////////////////////////////////////////////////////////
// destination = source
#define _copy(destination,source,n) memcpy(source,source,n*8);

///////////////////////////////////////////////////////////////////////////////////////
// out = sum( a[i] * b[i])
#define _dot(a,b,out,n)	out=0.0; for(int i=0; i<n;++i) out += a[i] * b[i];

///////////////////////////////////////////////////////////////////////////////////////
// c[n,k] = a[m,n]t x b[m,k]
#define _mult_t1(a,b,c,n,m,k) \
memset(c,0,n*k*8);\
for(int in=0; in<n;++in) \
	for(int ik=0; ik < k; ++ik) \
		for(int im=0; im < m; ++im) \
			c[in + n * ik] += a[im + m * in] * b[im + m * ik];

///////////////////////////////////////////////////////////////////////////////////////
// c[n,k] = a[n,m] x b[m,k]
#define _mult(a,b,c,n,m,k) \
memset(c,0,n*k*8);\
for(int in=0; in < n;++in) \
	for(int ik=0; ik < k; ++ik) \
		for(int im=0; im < m; ++im) \
			c[in + n * ik] += a[in + n * im] * b[im + m* ik];

///////////////////////////////////////////////////////////////////////////////////////
// c[n,k] = a[n,m] x b[k,m]t
#define _mult_t2(a,b,c,n,m,k) \
memset(c,0,n*k*8);\
for(int in=0; in<n;++in) \
	for(int ik=0; ik < k; ++ik) \
		for(int im=0; im < m; ++im) \
			c[in + n * ik] += a[in + n * im] * b[ik + k * im];

///////////////////////////////////////////////////////////////////////////////////////
// c[n,k] = a[m,n]t x b[m,k]
#define _mult_t1_plus(a,b,c,n,m,k) \
for(int in=0; in<n;++in) \
	for(int ik=0; ik < k; ++ik) \
		for(int im=0; im < m; ++im) \
			c[in + n * ik] += a[im + m * in] * b[im + m * ik];

///////////////////////////////////////////////////////////////////////////////////////
// c[n,k] = a[n,m] x b[m,k]
#define _mult_plus(a,b,c,n,m,k) \
for(int in=0; in < n;++in) \
	for(int ik=0; ik < k; ++ik) \
		for(int im=0; im < m; ++im) \
			c[in + n * ik] += a[in + n * im] * b[im + m* ik];

///////////////////////////////////////////////////////////////////////////////////////
// c[n,k] = a[n,m] x b[k,m]t
#define _mult_t2_plus(a,b,c,n,m,k) \
for(int in=0; in<n;++in) \
	for(int ik=0; ik < k; ++ik) \
		for(int im=0; im < m; ++im) \
			c[in + n * ik] += a[in + n * im] * b[ik + k * im];

///////////////////////////////////////////////////////////////////////////////////////
// c = a + b
#define _add(a,b,c,n) for(int i=0;i<n;++i) c[i] = a[i]+b[i];

///////////////////////////////////////////////////////////////////////////////////////
// c = a - b
#define _sub(a,b,c,n) for(int i=0;i<n;++i) c[i] = a[i]-b[i];

///////////////////////////////////////////////////////////////////////////////////////
// c = a + p1 x b
#define _add_p1(a,b,c,p,n) for(int i=0;i<n;++i) c[i] = a[i]+ p * b[i];

///////////////////////////////////////////////////////////////////////////////////////
// a += b
#define _add_plus(a,b,n) for(int i=0;i<n;++i) a[i] += b[i];

///////////////////////////////////////////////////////////////////////////////////////
// a -= b
#define _add_minus(a,b,n) for(int i=0;i<n;++i) a[i] -= b[i];

///////////////////////////////////////////////////////////////////////////////////////
// a += p x b
#define _add_plus_p1(a,b,p,n) for(int i=0;i<n;++i) a[i] += p * b[i];

///////////////////////////////////////////////////////////////////////////////////////
// a[i] = p[i] x a
#define _mult_p(a,p,n) for(int i=0;i<n;++i) a[i] = p * a[i];

///////////////////////////////////////////////////////////////////////////////////////
#define hradians(a) a * 0.017453292519943
#define degrees(a) 	a * 57.29577951308232087679815481410
#define radians(a) 	a * 0.01745329251994329576923690768489
#define dsin(a) 	sin(a*0.01745329251994329576923690768489)
#define dcos(a) 	cos(a*0.01745329251994329576923690768489)
#define dtan(a) 	tan(a*0.01745329251994329576923690768489)
//#define pi 3.141592653589793238462643383279

#define citr(i,n) if(i >= n) i -= n; else if(i < 0) i += n;
//#define citr(i,n)  (i>n)?i-=n:(i<0) i+=n:i
