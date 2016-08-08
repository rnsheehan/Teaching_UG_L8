#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of functions declared in the namespace file_funcs
// R. Sheehan 31 - 1 - 2013

double *file_funcs::read_vector_from_file(string filename, int &size)
{
	// Use this function to read a single column of data out of a file and into memory
	// R. Sheehan 31 - 1 - 2013

	// Read the numeric contents of a file into a vector
	// It is assumed that a single column of numbers is stored in the file
	// R. Sheehan 22 - 2 - 2012
	
	ifstream read;
	read.open(filename.c_str(),ios_base::in); // open the file for reading
	
	double *arr_ptr = nullptr; // declare a pointer to an array and assign it to nullptr

	if(read.is_open()){
		// Since you are assuming a single column of numerical data you can use the stream extraction operator
				
		// First item is to count the number of data points in the file, again using stream operators
		size = -1;

		// Count the number of lines in the file
		while(read.ignore(1280,'\n')){
			size++;
		}
		
		read.clear(); // empty the buffer
		read.seekg(0,ios::beg); // move to the start of the file

		arr_ptr = new (double [size+1] ); // create array to hold the data from the file
		
		// loop over the lines and read the data into memory
		for(int i=1; i <= size; i++){
			read>>arr_ptr[i];
		}
		
		read.close(); // close the file holding the data
	}
	else{
		cout<<"Error: could not open "<<filename<<endl;
	}

	return arr_ptr; // return the data as an array of doubles
}

void file_funcs::write_vector_to_file(string filename, double *vec, int size)
{
	// Use this function to write a single column of data to a file
	// R. Sheehan 31 - 1 - 2013

	ofstream write; // create the object used to write the data
	write.open(filename.c_str(), ios_base::out|ios_base::trunc); 

	if(write.is_open()){
		
		for(int i=1; i<=size; i++){
			write<<setprecision(15)<<vec[i]<<endl;
		}

	}
	else{
		cout<<"Could not open "<<filename<<" for writing\n";
	}
}