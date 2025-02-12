// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file hdf2tiff.c
 * 
 * \a hdf2tiff.exe - Extract the data from a HDF4 file and save as 
 * a tiff (grey-scale, 16 bit file). 
 * 
 * Usage: hdf2tiff.exe \<input hdf file\> \<output tiff file\> \<optional: data block name\> \par 
 *
 * \par Example:
 * \verbatim  hdf2tiff.exe my_data.hdf my_data.tiff  \endverbatim
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include <Double_2D.h>
#include <io.h>

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 2 or 3 arguements
  if(argc<3 || argc>4 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "hdf2tiff <input hdf file> <output tiff file>" << endl;
    return 1;
  }

  //read the data block in the file
  Double_2D data;
  int status;

  //read the data into an array
  if(argc==4)
    status = read_hdf4(argv[1], data, argv[3]);
  else
    status = read_hdf4(argv[1], data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_tiff(argv[2], data);
      
  return 0;
}
