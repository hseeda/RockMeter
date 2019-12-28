#pragma once

///////////////////////////////////////////////////////////////////////////////////////

#define pr(a) std::cout << a << "\n";

#define pr1_(a) std::cout << "<<<<< " << a << "\n";

#define pr2_(a1,a2) \
 std::cout << "<<<<< " << a1\
 << a2\
 << "\n";

#define pr3_(a1,a2,a3) \
 std::cout << "<<<<< " << a1\
 <<  a2\
 <<  a3\
 << "\n";

#define pr4_(a1,a2,a3,a4) \
 std::cout << "<<<<< " << a1\
 <<  a2\
 <<  a3\
 <<  a4\
 << "\n";

#define pr5_(a1,a2,a3,a4,a5) \
 std::cout << "<<<<< " << a1\
 <<  a2\
 <<  a3\
 <<  a4\
 <<  a5\
 << "\n";

///////////////////////////////////////////////////////////////////////////////////////
// print separator
#define prsc(ichar,a)\
		{int n = 16 * a;\
		for(int i=0;i<n; ++i)\
			std::cout << ichar ;\
		std::cout << "\n";}

///////////////////////////////////////////////////////////////////////////////////////
// print separator
#define prs(a) for(int i=0;i<a; ++i)  std::cout << "---------------|"; std::cout<< "\n";
#define prss(s,a)\
 for(int i=0;i<a; ++i){\
 for(int j=0;j<s-1; ++j)\
 std::cout << "-";\
 std::cout << "|";}\
 std::cout<< "\n";

///////////////////////////////////////////////////////////////////////////////////////
// print end separator --------
#define prs_(a) for(int i=0;i<a; ++i)  std::cout << "----------------"; std::cout<< "\n";

///////////////////////////////////////////////////////////////////////////////////////
// print end separator --------
#define prs_s(a) for(int i=0;i<a; ++i)  std::cout << "****************"; std::cout<< "\n";

///////////////////////////////////////////////////////////////////////////////////////
// print end separator ========
#define prs_e(a) for(int i=0;i<a; ++i)  std::cout << "================"; std::cout<< "\n";

///////////////////////////////////////////////////////////////////////////////////////
// print end separator ====
#define prs_eh(a) for(int i=0;i<a; ++i)  std::cout << "========"; std::cout<< "\n";


///////////////////////////////////////////////////////////////////////////////////////
// print end separator <<<<<<<<
#define prs_l(a) for(int i=0;i<a; ++i)  std::cout << "<<<<<<<<<<<<<<<<"; std::cout<< "\n";

///////////////////////////////////////////////////////////////////////////////////////
//print parameter
#define prp(a1,a2)  std::cout << std::setw(32) << std::left << a1 << " = " << std::setw(13) << std::left << a2 << "\n";

///////////////////////////////////////////////////////////////////////////////////////
//print message
#define prm(a,n)  \
		for(int i=0;i<n; ++i)  \
			std::cout << "---------------|"; \
		std::cout << " [ " << a << " ]" << "\n";

///////////////////////////////////////////////////////////////////////////////////////
// print vector
#define prv(v,n) \
	for(int i=0;i<n; ++i) \
		std::cout << std::showpos << std::setw(16) << v[i];\
		std::cout << std::noshowpos << "\n";
#define prvs(v,n,separator) \
	for(int i=0;i<n-1; ++i) \
		std::cout  << v[i] << separator;\
		std::cout  << v[n-1] <<"\n";
///////////////////////////////////////////////////////////////////////////////////////
// print double array (scientific)
#define prv_s(v,n) \
	for(int i=0;i<n; ++i) \
		std::cout << std::scientific  << std::left << std::showpos << std::setw(16) << v[i]; \
		std::cout << std::fixed << std::noshowpos << "\n";


#define pr1(a1)   std::cout << a1 << "\n";
#define pr2(a1,a2)   std::cout << a1 << a2 << "\n";
#define pr3(a1,a2,a3)   std::cout << a1 << a2 << a3 << "\n";
#define pr4(a1,a2,a3,a4)   std::cout << a1 << a2 << a3 << a4 << "\n";
#define pr5(a1,a2,a3,a4,a5)   std::cout << a1 << a2 <<a3 << a4 << a5 << "\n";
#define pr6(a1,a2,a3,a4,a5,a6)   std::cout << a1 << a2 <<a3 << a4 << a5 << a6 << "\n";
#define pr7(a1,a2,a3,a4,a5,a6,a7)   std::cout << a1 << a2 <<a3 << a4 << a5 << a6 << a7 << "\n";
#define pr8(a1,a2,a3,a4,a5,a6,a7,a8)   std::cout << a1 << a2 <<a3 << a4 << a5 << a6 <<  a7 <<  a8 << "\n";
#define pr9(a1,a2,a3,a4,a5,a6,a7,a8,a9)   std::cout << a1 << a2 <<a3 << a4 << a5 << a6 <<  a7 <<  a8 << a9 << "\n";



#define pr2c(a1,a2)   std::cout << a1 << " , " << a2 << "\n";
#define pr3c(a1,a2,a3)   std::cout << a1 << " , " << a2 << " , " << a3 << "\n";
#define pr4c(a1,a2,a3,a4)   std::cout << a1 << " , " << a2 << " , "<< a3 << " , " << a4 << "\n";
#define pr5c(a1,a2,a3,a4,a5)   std::cout << a1 << " , " << a2 <<  " , " << a3 <<  " , " << a4 <<  " , " << a5 << "\n";
#define pr6c(a1,a2,a3,a4,a5,a6)   std::cout << a1 <<  " , " << a2 <<  " , " << a3 <<  " , " << a4 <<  " , " << a5 <<  " , " << a6 << "\n";
#define pr7c(a1,a2,a3,a4,a5,a6,a7)   std::cout << "<<<<< "<< a1 << " , " << a2<< " , "  << a3<< " , "  << a4 << " , " << a5<< " , "  << a6 << " , " << a7 << "\n";



#define pr2t(a1,a2)   std::cout << a1 << " \t" << a2 << " \t" <<"\n";
#define pr3t(a1,a2,a3)   std::cout << a1 << " \t " << a2 << " \t " << a3 << "\n";
#define pr4t(a1,a2,a3,a4)   std::cout << a1 << " \t" << a2 << " \t" << a3 << " \t" << a4 << "\n";
#define pr5t(a1,a2,a3,a4,a5)   std::cout << a1 << " \t" << a2 <<  " \t" << a3 <<  " \t"  << a4 <<  " \t" << a5 << "\n";
#define pr6t(a1,a2,a3,a4,a5,a6)   std::cout << a1 <<  " \t" << a2 <<  " \t" << a3 <<  " \t" << a4 <<  " \t" << a5 <<  " \t" << a6 << "\n";
#define pr7t(a1,a2,a3,a4,a5,a6,a7)   std::cout << a1 <<  " \t" << a2 <<  " \t" << a3 <<  " \t" << a4 <<  " \t" << a5 <<  " \t" << a6 <<  " \t" << a7 << "\n";

#define pr9f(w,a1,a2,a3,a4,a5,a6,a7,a8,a9)\
	cout_reset(); \
   std::cout \
	<< std::setw(w) << std::left << a1\
	<< std::setw(w) << std::left << a2\
	<< std::setw(w) << std::left << a3\
	<< std::setw(w) << std::left << a4\
	<< std::setw(w) << std::left << a5\
	<< std::setw(w) << std::left << a6\
	<< std::setw(w) << std::left << a7\
	<< std::setw(w) << std::left << a8\
	<< std::setw(w) << std::left << a9\
	<< "\n";



#define set_cout_precision(a) std::cout.precision(a);
#define set_cout_width(a) std::cout.width(a);
#define set_cout_fixed() std::cout.setf(std::ios::fixed);



#define prstr1(name,a1)\
{std::stringstream s; s << a1; name = s.str();}
#define prstr2(name,a1,a2)\
{std::stringstream s; s << a1 << a2; name = s.str();}
#define prstr3(name,a1,a2,a3)\
{std::stringstream s; s << a1 << a2 << a3; name = s.str();}
#define prstr4(name,a1,a2,a3,a4)\
{std::stringstream s; s << a1 << a2 << a3 << a4; name = s.str();}
#define prstr5(name,a1,a2,a3,a4,a5)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5; name = s.str();}
#define prstr6(name,a1,a2,a3,a4,a5,a6)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6; name = s.str();}
#define prstr7(name,a1,a2,a3,a4,a5,a6,a7)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7; name = s.str();}
#define prstr8(name,a1,a2,a3,a4,a5,a6,a7,a8)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8; name = s.str();}
#define prstr9(name,a1,a2,a3,a4,a5,a6,a7,a8,a9)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9; name = s.str();}
#define prstr10(name,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 << a10; name = s.str();}

#define prtxt1(a1)\
{std::stringstream s; s << a1; txt = s.str();}
#define prtxt2(a1,a2)\
{std::stringstream s; s << a1 << a2; txt = s.str();}
#define prtxt3(a1,a2,a3)\
{std::stringstream s; s << a1 << a2 << a3; txt = s.str();}
#define prtxt4(a1,a2,a3,a4)\
{std::stringstream s; s << a1 << a2 << a3 << a4; txt = s.str();}
#define prtxt5(a1,a2,a3,a4,a5)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5; txt = s.str();}
#define prtxt6(a1,a2,a3,a4,a5,a6)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6; txt = s.str();}
#define prtxt7(a1,a2,a3,a4,a5,a6,a7)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7; txt = s.str();}
#define prtxt8(a1,a2,a3,a4,a5,a6,a7,a8)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8; txt = s.str();}
#define prtxt9(a1,a2,a3,a4,a5,a6,a7,a8,a9)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9; txt = s.str();}
#define prtxt10(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)\
{std::stringstream s; s << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 << a10; txt = s.str();}
