#pragma once

///////////////////////////////////////////////////////////////////////////////////////
#define ifd(a)  if (HDEBUG) {{a}}


#define ifd1(a) if (HDEBUG1) {{a}}
#define ifd2(a) if (HDEBUG2) {{a}}

#define hexit getchar();exit(9999)
#define hpause std::cin.get()

#define frewind(f)			f.clear();f.seekg(0);
#define feof(f) 				f.clear();f.seekp(0, std::ios_base::end);
#define fclosef(f) 			f.close();

#define fopenI(f,filename)	f.open(filename, std::ios::in | std::ios::binary);
#define fopenO(f,filename)	f.open(filename, std::ios::out | std::ios::binary);
#define fopenA(f,filename)	f.open(filename, std::ios::app | std::ios::binary);
#define fopenIt(f,filename)	f.open(filename, std::ios::in);
#define fopenOt(f,filename)	f.open(filename, std::ios::out);
#define fopenAt(f,filename)	f.open(filename, std::ios::app);

#define MSGBOX(x) \
{ \
   std::ostringstream oss; \
   oss << x; \
   MessageBox(NULL,oss.str().c_str(), "Msg Title", MB_OK | MB_ICONQUESTION); \
}
