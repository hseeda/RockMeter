#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <thread>
#include <vector>
#include <chrono>
#include <conio.h>

#include <cerrno>
#include <ctime>


#include <locale>
#include <map>
#include <memory>
#include <set>
#include <sstream>


#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <iterator>  // std::begin, std::end

#include "ConsoleColor.h"
#include "windows.h"
#include  <cstdlib>
#include "Defines/h_print.h"
#include "Defines/h_defines.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::chrono::steady_clock::time_point tictoc12345[10];
static std::chrono::steady_clock::time_point tictoc1234567890[10];
static std::ios::fmtflags os_flags(std::cout.flags());
static const std::string WHITESPACE = " \n\r\t\f\v";
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//#                                  MKL
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//static void mklVersion() {
//	MKLVersion Version;
//	mkl_get_version(&Version);
//	printf("======================================================================================\n");
//	printf("Major version: %d\n", Version.MajorVersion);
//	printf("Minor version: %d\n", Version.MinorVersion);
//	printf("Update version: %d\n", Version.UpdateVersion);
//	printf("Product status: %s\n", Version.ProductStatus);
//	printf("Build: %s\n", Version.Build);
//	printf("Platform: %s\n", Version.Platform);
//	printf("Processor optimization: %s\n", Version.Processor);
//	printf("======================================================================================\n");
//	printf("\n");
//}

//static void setMKLThreads(int mklth) {
//	mkl_set_dynamic(false);
//	mkl_set_num_threads_local(mklth);
//}

//static int getMKLThreads() {
//	return mkl_get_max_threads();
//}

static int getMAXThreads() {
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	return concurentThreadsSupported;
}

//static void startMKLmem_usage_monitoring() {
//	mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE);
//}


//static void freeMKL(int printFlag) {
//	if (printFlag != 0) {
//		std::cout << "MKL Peak Memory Usage    (" << dFormatWithCommas((double)mkl_peak_mem_usage(MKL_PEAK_MEM) / 1000) << " KBytes)" << std::endl;
//		//std::cout << "MKL Peak Memory Usage    (" << mkl_peak_mem_usage(MKL_PEAK_MEM) / 1000 << " KBytes)" << std::endl;
//		MKL_INT64	AllocatedBytes;
//		int			N_AllocatedBuffers;
//		AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
//		std::cout << "MKL used    (" << dFormatWithCommas((double)AllocatedBytes / 1000) << "  KBytes) in (" << N_AllocatedBuffers << ") Buffers\n";
//		//std::cout << "MKL used    (" << AllocatedBytes / 1000 << "  KBytes) in (" << N_AllocatedBuffers << ") Buffers\n";
//	}
//
//	mkl_free_buffers();
//
//	if (printFlag != 0) {
//		MKL_INT64	AllocatedBytes;
//		int			N_AllocatedBuffers;
//		AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
//		if (AllocatedBytes > 0) {
//			std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MKL memory leak!\n";
//			std::cout << "MKL used             (" << dFormatWithCommas((double)AllocatedBytes / 1000) << "  KBytes) in (" << N_AllocatedBuffers << ") Buffers\n";
//			//std::cout << "MKL used             (" << AllocatedBytes / 1000 << "  KBytes) in (" << N_AllocatedBuffers << ") Buffers\n";
//			getchar();
//			exit(9999);
//		}
//		std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MKL memory FREED!\n";
//	}
//}

static int addToPath(std::string path) {
	const std::size_t ENV_BUF_SIZE = 10000; // Enough for your PATH?

	char buf[ENV_BUF_SIZE];
	std::size_t bufsize = ENV_BUF_SIZE;
	int e = getenv_s(&bufsize, buf, bufsize, "PATH");
	if (e) {
		std::cerr << "`getenv_s` failed, returned " << e << '\n';
		exit(EXIT_FAILURE);
	}
	std::string env_path = buf;
	//std::cout << "In main process, `PATH`=" << env_path << std::endl;
	env_path += ";" + path;
	e = _putenv_s("PATH", env_path.c_str());
	if (e) {
		std::cerr << "`_putenv_s` failed, returned " << e << std::endl;
		exit(EXIT_FAILURE);
	}
	//e = std::system("echo \"In child process `PATH`=%PATH%\"");
	//if (e) {
	//	std::cerr << "`std::system` failed, returned " << e << std::endl;
	//	exit(EXIT_FAILURE);
	//}
	return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void createDirectory(const std::string& rootDirectory, const std::string& dirName) {
	std::string s = rootDirectory + "\\" + dirName;
	CreateDirectoryA((s.c_str()), NULL);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::string getJobName() {
	char result[MAX_PATH];
	GetModuleFileNameA(NULL, result, MAX_PATH);
	std::string fileName(result);
	size_t position = fileName.find(".");
	return (std::string::npos == position) ? fileName : fileName.substr(0, position);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::string getCurrentDirectory() {
	char result[MAX_PATH];
	GetCurrentDirectoryA(MAX_PATH, result);
	return std::string(result);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int delFile(const std::string& filePath, const std::string& fileNameWOext, const std::string& fileExt) {
	// Set file attributes
	std::string s = filePath + "//" + fileNameWOext + "." + fileExt;
	if (::SetFileAttributesA((s.c_str()), FILE_ATTRIBUTE_NORMAL) == FALSE)
		return ::GetLastError();

	// Delete file
	if (::DeleteFileA(s.c_str()) == FALSE)
		return ::GetLastError();
	else
		return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int delFile(const std::string& filePath, const std::string& fileName) {
	// Set file attributes
	std::string s = filePath + "//" + fileName;
	if (::SetFileAttributesA((s.c_str()), FILE_ATTRIBUTE_NORMAL) == FALSE)
		return ::GetLastError();

	// Delete file
	if (::DeleteFileA(s.c_str()) == FALSE)
		return ::GetLastError();
	else
		return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int delFile(const std::string& filePath) {
	// Set file attributes
	if (::SetFileAttributesA((filePath.c_str()), FILE_ATTRIBUTE_NORMAL) == FALSE)
		return ::GetLastError();

	// Delete file
	if (::DeleteFileA((filePath.c_str())) == FALSE)
		return ::GetLastError();
	else
		return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int delFilesinDirectory(const std::string& rootDirectory)
{
	bool            bSubdirectory = false;       // Flag, indicating whether
												 // subdirectories have been found
	HANDLE          hFile;                       // Handle to directory
	std::string     filePath;                 // Filepath
	std::string     strPattern;                  // Pattern
	WIN32_FIND_DATAA FileInformation;             // File information

	strPattern = rootDirectory + "\\*.*";
	hFile = ::FindFirstFileA(strPattern.c_str(), &FileInformation);
	if (hFile != INVALID_HANDLE_VALUE)
	{
		do
		{
			if (FileInformation.cFileName[0] != '.')
			{
				filePath.erase();
				filePath = rootDirectory + "\\" + FileInformation.cFileName;

				if (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
				{
					bSubdirectory = true;
				}
				else
				{
					// Set file attributes
					if (::SetFileAttributesA(filePath.c_str(),
						FILE_ATTRIBUTE_NORMAL) == FALSE)
						return ::GetLastError();

					// Delete file
					if (::DeleteFileA(filePath.c_str()) == FALSE)
						return ::GetLastError();
				}
			}
		} while (::FindNextFileA(hFile, &FileInformation) == TRUE);

		// Close handle
		::FindClose(hFile);

		DWORD dwError = ::GetLastError();
		if (dwError != ERROR_NO_MORE_FILES)
			return dwError;
		else
		{
			return 0;
		}
	}
	return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int delFilesinDirectory(const std::string& rootDirectory, const std::string& fileExt)
{
	bool            bSubdirectory = false;       // Flag, indicating whether
												 // subdirectories have been found
	HANDLE          hFile;                       // Handle to directory
	std::string     filePath;                 // Filepath
	std::string     strPattern;                  // Pattern
	WIN32_FIND_DATAA FileInformation;             // File information

	strPattern = rootDirectory + "\\*." + fileExt;
	hFile = ::FindFirstFileA(strPattern.c_str(), &FileInformation);
	if (hFile != INVALID_HANDLE_VALUE)
	{
		do
		{
			if (FileInformation.cFileName[0] != '.')
			{
				filePath.erase();
				filePath = rootDirectory + "\\" + FileInformation.cFileName;

				if (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
				{
					bSubdirectory = true;
				}
				else
				{
					// Set file attributes
					if (::SetFileAttributesA((filePath.c_str()),
						FILE_ATTRIBUTE_NORMAL) == FALSE)
						return ::GetLastError();

					// Delete file
					if (::DeleteFileA((filePath.c_str())) == FALSE)
						return ::GetLastError();
				}
			}
		} while (::FindNextFileA(hFile, &FileInformation) == TRUE);

		// Close handle
		::FindClose(hFile);

		DWORD dwError = ::GetLastError();
		if (dwError != ERROR_NO_MORE_FILES)
			return dwError;
		else
		{
			return 0;
		}
	}
	return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int delDirectory(const std::string& rootDirectory,
	bool              bDeleteSubdirectories)
{
	bool            bSubdirectory = false;       // Flag, indicating whether
												 // subdirectories have been found
	HANDLE          hFile;                       // Handle to directory
	std::string     filePath;                 // Filepath
	std::string     strPattern;                  // Pattern
	WIN32_FIND_DATAA FileInformation;             // File information

	strPattern = rootDirectory + "\\*.*";
	hFile = ::FindFirstFileA((strPattern.c_str()), &FileInformation);
	if (hFile != INVALID_HANDLE_VALUE)
	{
		do
		{
			if (FileInformation.cFileName[0] != '.')
			{
				filePath.erase();
				filePath = rootDirectory + "\\" + FileInformation.cFileName;

				if (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
				{
					if (bDeleteSubdirectories)
					{
						// Delete subdirectory
						int iRC = delDirectory(filePath, bDeleteSubdirectories);
						if (iRC)
							return iRC;
					}
					else
						bSubdirectory = true;
				}
				else
				{
					// Set file attributes
					if (::SetFileAttributesA((filePath.c_str()),
						FILE_ATTRIBUTE_NORMAL) == FALSE)
						return ::GetLastError();

					// Delete file
					if (::DeleteFileA((filePath.c_str())) == FALSE)
						return ::GetLastError();
				}
			}
		} while (::FindNextFileA(hFile, &FileInformation) == TRUE);

		// Close handle
		::FindClose(hFile);

		DWORD dwError = ::GetLastError();
		if (dwError != ERROR_NO_MORE_FILES)
			return dwError;
		else
		{
			if (!bSubdirectory)
			{
				// Set directory attributes
				if (::SetFileAttributesA((rootDirectory.c_str()),
					FILE_ATTRIBUTE_NORMAL) == FALSE)
					return ::GetLastError();

				// Delete directory
				if (::RemoveDirectoryA((rootDirectory.c_str())) == FALSE)
					return ::GetLastError();
			}
		}
	}

	return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static int numCPUs() {
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	return concurentThreadsSupported;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void tic_h(int itimer = 0) {
	tictoc12345[itimer] = std::chrono::high_resolution_clock::now();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void tic(int itimer = 0) {
	tictoc12345[itimer] = std::chrono::high_resolution_clock::now();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void ticp(int itimer = 0) {
	tictoc12345[itimer] = std::chrono::high_resolution_clock::now();
	std::cout << "Timer (" << itimer << ") Started." << std::endl;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void tic_l(int itimer = 0) {
	tictoc12345[itimer] = std::chrono::steady_clock::now();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static double toc_l(int itimer = 0) {
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	return std::chrono::duration_cast<std::chrono::milliseconds> (end - tictoc12345[itimer]).count() / 1000.0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static double toc_h(int itimer = 0) {
	auto end = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::microseconds> (end - tictoc12345[itimer]).count() / 1000000.0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static double toc(int itimer = 0) {
	auto end = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::microseconds> (end - tictoc12345[itimer]).count() / 1000000.0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void tocp(std::string s, int itimer = 0) { std::cout.precision(4); std::cout << std::setw(40) << std::left << s << std::setw(8) << std::setprecision(4) << toc(itimer) << " sec" << std::endl; }
static void tocp(int itimer = 0) { std::cout.precision(4); std::cout << std::setw(40) << std::left << "Elapsed Time in seconds = " << std::setw(8) << std::setprecision(4) << toc(itimer) << " sec" << std::endl; }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void tocpr(std::string s, int itimer = 0) { std::cout << std::setw(40) << std::left << s << std::setw(8) << std::setprecision(4) << toc(itimer) << " sec" << std::endl; tic(itimer); }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void hsep(std::string s, char se, int n)
{
	for (int i = 0; i < n; i++) std::cout << se;
	std::cout << s << std::endl;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef Hpprint
#define Hpprint
template <class T>
void mprint(T* h, int nrow, int ncol) {
	int c = 0;
	std::cout.precision(5);;
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			std::cout << "  " << std::left << std::setw(2) << h[c];
			c += 1;
		}
		std::cout << std::endl;
	}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <class T>
void rprint(T* h, int length)									// print vector in row format
{
	std::cout.precision(5);
	for (int i = 0; i < length; i++) {
		std::cout << "  " << std::left << std::setw(2) << *(h + i);
	}
	std::cout << std::endl;
}

template <class T>
void cprint(T* h, int length)									// print vector in col format
{
	std::cout.precision(5);
	for (int i = 0; i < length; i++) {
		std::cout << *(h + i) << std::endl;
	}
}
#endif // !Hpprint
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//--------------------      File handlers                 ---------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::ifstream openRead(std::string file_name) {
	std::ifstream f(file_name, std::ios::in);
	return f;
}
int parseLine(std::ifstream* f, char* delim, std::vector<std::string> &v, std::string &sline) {
	int c = 0;
	char* token;
	char* tmp;
	//--------------------------------------------------------- Read line from file
	
	std::getline(*f, sline);
	if (sline == "") v[0] = "";

//#ifdef _DEBUG
//	std::cout << sline << std::endl;
//	std::cout << "-------------------------- " << c << std::endl;
//#endif
	token = strtok_s((char*)sline.c_str(), delim, &tmp);
	while (token)
	{
		v[c] = token;
		c = c + 1;
		//std::cout << token << " ; ";
		token = strtok_s(NULL, delim, &tmp);
	}
	//std::cout << std::endl;
	return c;
}
static std::string iFormatWithCommas(int value)
{
	std::stringstream ss;
	ss.imbue(std::locale(""));
	ss << std::fixed << value;
	return ss.str();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::string dFormatWithCommas(double value)
{
	std::stringstream ss;
	ss.imbue(std::locale(""));
	ss << std::fixed << value;
	return ss.str();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static void printSeparator() { std::cout << "---------|---------|---------|---------|---------|---------|---------|---------|\n"; }


static void geotechnics() {
	pr("================================================================");
	std::cout << red;
	std::cout << "                         __            __          _          " << std::endl;
	std::cout << white << "       ____ " << red;
	std::cout << "____  ____  / /____  _____/ /_  ____  (_)_________" << std::endl;
	std::cout << white << "      / __ `" << red;
	std::cout << "/ _ \\/ __ \\/ __/ _ \\/ ___/ __ \\/ __ \\/ / ___/ ___/" << std::endl;
	std::cout << white << "     / /_/ " << red;
	std::cout << "/  __/ /_/ / /_/  __/ /__/ / / / / / / / /__(__  ) " << std::endl;
	std::cout << white << "     \\__, /" << red;
	std::cout << "\\___/\\____/\\__/\\___/\\___/_/ /_/_/ /_/_/\\___/____/  " << std::endl;
	std::cout << white;
	std::cout << "    /____/                                                    " << std::endl;
}

static void finish() {
	pr("================================================================");

	pr("                ____ _         _        __  ");
	pr("               / __/(_)____   (_)_____ / /_ ");
	pr("              / /_ / // __ \\ / // ___// __ \\");
	pr("             / __// // / / // /(__  )/ / / /");
	pr("            /_/  /_//_/ /_//_//____//_/ /_/ ");
	std::cout << cyan;
	pr("================================================================")
}

static void debug() {
	pr("    ____   ______ ____   __  __ ______   __  ___ ____   ____   ______");
	pr("   / __ \\ / ____// __ ) / / / // ____/  /  |/  // __ \\ / __ \\ / ____/");
	pr("  / / / // __/  / __  |/ / / // / __   / /|_/ // / / // / / // __/   ");
	pr(" / /_/ // /___ / /_/ // /_/ // /_/ /  / /  / // /_/ // /_/ // /___  ");
	pr("/_____//_____//_____/ \\____/ \\____/  /_/  /_/ \\____//_____//_____/  ");
	pr("");
}
static void fdeser(std::string file_name,
	std::vector<std::string>& name_list,
	std::vector<int>* i1_list,
	std::vector<int>* i2_list,
	std::vector<int>* col_list,
	std::vector<int>* row_list)
{
	int found = 0;
	std::fstream f;
	fopenI(f, file_name);
	if (f.is_open()) {
		std::streampos pos;
		size_t idsize;
		int irow, icol, varsize, nextPos;
		int i1, i2;
		std::string iname;
		while (found == 0) {
			std::string iid(80, '\0');
			f.read((char*)&idsize, sizeof(size_t));
			f.read((char*)&iid[0], idsize);
			f.read((char*)&irow, sizeof(int));
			f.read((char*)&icol, sizeof(int));
			f.read((char*)&varsize, sizeof(int));
			f.read((char*)&i1, sizeof(int));
			f.read((char*)&i2, sizeof(int));
			f.read((char*)&nextPos, sizeof(int));

			if (f.eof()) {
				fclosef(f);
				return;
			}
			iname = iid.substr(0, idsize);
			//std::cout << "[" << iname << "]\n";
			name_list.push_back(iname);
			if (i1_list != nullptr) i1_list->push_back(i1);
			if (i2_list != nullptr) i2_list->push_back(i2);
			if (col_list != nullptr) col_list->push_back(icol);
			if (row_list != nullptr) row_list->push_back(irow);

			pos = f.tellg();
			pos += nextPos;
			f.seekg(pos);
		}
		fclosef(f);
	}
	fclosef(f);
}
/////////////////////////////////////////////////////////////////////////////////////////////
static void fdeser(std::string file_name)
{
	std::vector<std::string> name_list;
	std::vector<int> i1_list;
	std::vector<int> i2_list;
	std::vector<int> col_list;
	std::vector<int> row_list;

	fdeser(file_name, name_list, &i1_list, &i2_list, &col_list, &row_list);

	pr5_("File [", file_name, "] contains [", name_list.size(), "] variables");
	//for (auto i = name_list.begin(); i != name_list.end(); ++i) {
	//	pr(*i);
	//}
	std::cout << "---------------------------------------------------------------" << std::endl;
	std::cout << std::left << std::setw(32) << "  Variable Name";
	std::cout << std::left << std::setw(13) << " [size]" << " \ti1  \ti2" << std::endl;
	std::cout << "---------------------------------------------------------------" << std::endl;

	int c = (int)name_list.size();
	for (int i = 0; i < c; i++) {
		std::cout << "  " << std::left << std::setw(30) << name_list[i];
		std::cout << " [" << std::left << std::setw(6) << row_list[i] << "," << std::left << std::setw(2) << col_list[i];
		std::cout << "] \t" << i1_list[i] << " \t" << i2_list[i] << std::endl;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////
static int fdeserN(std::string file_name)
{
	std::fstream f;
	int c = 0;
	fopenI(f, file_name);
	if (f.is_open()) {
		std::streampos pos;
		size_t idsize;
		int irow, icol, varsize, nextPos;
		int i1, i2;
		std::string iname;
		while (1) {
			std::string iid(80, '\0');
			f.read((char*)&idsize, sizeof(size_t));
			f.read((char*)&iid[0], idsize);
			f.read((char*)&irow, sizeof(int));
			f.read((char*)&icol, sizeof(int));
			f.read((char*)&varsize, sizeof(int));
			f.read((char*)&i1, sizeof(int));
			f.read((char*)&i2, sizeof(int));
			f.read((char*)&nextPos, sizeof(int));

			if (f.eof()) {
				fclosef(f);
				return c;
			}
			c++;
			pos = f.tellg();
			pos += nextPos;
			f.seekg(pos);
		}
		fclosef(f);
	}
	fclosef(f);
	return c;
}

///////////////////////////////////////////////////////////////////////////////////////////////
static double  _norm(double* v, int n) {
	double out = 0.0;
	//#pragma loop( ivdep )
	//#pragma loop( hint_parallel(8) )
	//#pragma loop( no_vector )
	for (int i = 0; i < n; ++i)
		out += pow(v[i], 2);
	return pow(out, 0.5);
}
//
/////////////////////////////////////////////////////////////////////////////////////////////
static std::string ltrim(const std::string& s)
{
	size_t start = s.find_first_not_of(WHITESPACE);
	return (start == std::string::npos) ? "" : s.substr(start);
}

/////////////////////////////////////////////////////////////////////////////////////////////
static std::string rtrim(const std::string& s)
{
	size_t end = s.find_last_not_of(WHITESPACE);
	return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

/////////////////////////////////////////////////////////////////////////////////////////////
static std::string trim(const std::string& s)
{
	return rtrim(ltrim(s));
}

/////////////////////////////////////////////////////////////////////////////////////////////
static int ftimeM(std::string filename)
{
	struct stat fileInfo;

	if (stat(filename.c_str(), &fileInfo) != 0) {  // Use stat( ) to get the info
		std::cerr << "Error: " << "ftimeM" << '\n';
		return(EXIT_FAILURE);
	}

	//std::cout << "Type:         : ";
	//if ((fileInfo.st_mode & S_IFMT) == S_IFDIR) { // From sys/types.h
	//	std::cout << "Directory\n";
	//}
	//else {
	//	std::cout << "File\n";
	//}

	//std::cout << "Size          : " <<
	//	fileInfo.st_size << '\n';               // Size in bytes

	//std::cout << "Device        : " <<
	//	(char)(fileInfo.st_dev + 'A') << '\n';  // Device number

	//std::cout << "Created       : " <<
	//	std::ctime(&fileInfo.st_ctime);         // Creation time

	//std::cout << "Modified      : " <<
	//	std::ctime(&fileInfo.st_mtime);         // Last mod time

	return (int)fileInfo.st_mtime;
}

/////////////////////////////////////////////////////////////////////////////////////////////
static int ftimeC(std::string filename)
{
	struct stat fileInfo;

	if (stat(filename.c_str(), &fileInfo) != 0) {  // Use stat( ) to get the info
		std::cerr << "Error: " << "ftimeC" << '\n';
		return(EXIT_FAILURE);
	}

	return (int)fileInfo.st_ctime;
}
/////////////////////////////////////////////////////////////////////////////////////////////
static void cout_reset() {
	std::cout.flags(os_flags);
}
///////////////////////////////////////////////////////////////////////////////////////////////
static void showConsoleCursor(bool showFlag)
{
	HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO     cursorInfo;
	GetConsoleCursorInfo(out, &cursorInfo);
	cursorInfo.bVisible = showFlag; // set the cursor visibility
	SetConsoleCursorInfo(out, &cursorInfo);
}

template<typename _T>
void writeCSV(std::string file_name, std::vector<_T>& v, int stride = 1) {
	std::fstream f;
	fopenOt(f, file_name);
	size_t n = v.size() / stride;
	int c = 0;
	for (size_t j = 0; j < n; j++) {
		for (int i = 0; i < stride - 1; i++) {
			f << v[c] << ",";
			c++;
		}
		f << v[c] << "\n";
		c++;
	}
	fclosef(f);
}


static bool devisableBy(int i, int n)
{
	if(i%n) return false;
	else return true;
}

static bool devisableBy(double i, DOUBLE n)
{
	if( abs(i/n - int(i/n)) > 1e-6 ) return false;
	else return true;
}