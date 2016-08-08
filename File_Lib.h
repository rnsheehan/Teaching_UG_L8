#ifndef FILE_LIB_H
#define FILE_LIB_H

//#ifndef ATTACH_H
//#include "Attach.h"
//#endif

// The functions in this file will be used to read data from files, and also to write data to files
// R. Sheehan 31 - 1 - 2013

namespace file_funcs{

	double *read_vector_from_file(string filename, int &size); // Use this function to read a vector out of a file and into memory

	void write_vector_to_file(string filename, double *vec, int size); // Use this function to write a vector to a file

}

#endif