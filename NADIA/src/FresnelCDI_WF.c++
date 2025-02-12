#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <cstdlib> 
#include <Complex_2D.h>
#include <FresnelCDI_WF.h>
#include <io.h> 
#include <sstream>

using namespace std;

#define FORWARD  -1
#define BACKWARD +1

FresnelCDI_WF::FresnelCDI_WF(Complex_2D & initial_guess,
			     double beam_wavelength,
			     double zone_focal_length,
			     double focal_detector_length,
			     double pixel_size,
			     int n_best)
  :BaseCDI(initial_guess,n_best),
   wavelength(beam_wavelength),
   zone_to_focal_length(zone_focal_length),
   focal_to_detector_length(focal_detector_length),
   pixel_length(pixel_size),
   coefficient(nx/2,ny/2)
   //   backward_coefficient(nx,ny)
{


  //set-up the coefficients
  //it's easier to do it once and reuse the matrix.

  double x_mid = (nx-1)/2.0;
  double y_mid = (ny-1)/2.0;

  double z12 = zone_to_focal_length;
  double z23 = focal_to_detector_length;

  double scaling_x = beam_wavelength*z23/(pixel_length*nx);
  double scaling_y = beam_wavelength*z23/(pixel_length*ny);

  double norm = 1/(sqrt(nx*ny));

  for(int i=0; i<nx/2; i++){
    for(int j=0; j<ny/2; j++){

      double rho_2 = pow(scaling_x*(x_mid-i),2) + pow(scaling_y*(y_mid-j),2);
      //double phi = 2*z12 + 2*z23 + rho_2/z12 + rho_2/z23 ;

      double phi = rho_2/z12 + rho_2/z23 ;

      phi*= M_PI/beam_wavelength;

      coefficient.set_real(i,j,sin(phi)*norm);
      coefficient.set_imag(i,j,-cos(phi)*norm);

    }
  }
  
  /**  Double_2D result(nx,ny);
  forward_coefficients_const->get_2d(PHASE,result);
  write_ppm("forward_coefficients_phase.ppm",result);**/

}


void FresnelCDI_WF::multiply_factors(Complex_2D & c, int direction){
  
  double x_mid = (nx-1)/2.0;
  double y_mid = (ny-1)/2.0;

  double old_real, old_imag;
  double coef_real, coef_imag;
  int i_, j_;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      old_real = c.get_real(i,j);
      old_imag = c.get_imag(i,j);

      int i_ = x_mid - fabs(x_mid - i);
      int j_ = y_mid - fabs(y_mid - j);

      coef_real = coefficient.get_real(i_,j_);
      coef_imag = direction*coefficient.get_imag(i_,j_);
      
      c.set_real(i,j,old_real*coef_real - old_imag*coef_imag);
      c.set_imag(i,j,old_imag*coef_real + old_real*coef_imag);
     
    }
  }
}


void FresnelCDI_WF::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(!support.get(i,j)){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	//	double random = (65000.0*rand()/(double) RAND_MAX);
	double mag = intensity_sqrt.get(i,j); 
	complex.set_value(i,j,REAL,mag); 
	complex.set_value(i,j,IMAG,0);
      }
    }
  }
}




int FresnelCDI_WF::iterate(){

  //Double_2D result(nx,ny);
  
  // complex.get_2d(PHASE,result);
  // write_ppm("b0.ppm",result);

  propagate_from_detector(complex);

  // complex.get_2d(PHASE,result);
  // write_ppm("b1.ppm",result);

  apply_support(complex);

  //  complex.get_2d(PHASE,result);
  // write_ppm("b3.ppm",result);
  
  propagate_to_detector(complex);

  //  complex.get_2d(PHASE,result);
  // write_ppm("b4.ppm",result);
  
  scale_intensity(complex);

  //  complex.get_2d(PHASE,result);
  // write_ppm("b5.ppm",result);

  return SUCCESS;
}

void FresnelCDI_WF::propagate_from_detector(Complex_2D & c){
  //go to the focal plane
  c.perform_backward_fft();
  c.invert(true);  
  multiply_factors(c, BACKWARD);
  //  c.multiply(backward_coefficient);

  //go back to zone plate plane. 
  c.perform_backward_fft();

}

void FresnelCDI_WF::propagate_to_detector(Complex_2D & c){

  //go to the focal plane again.
  c.perform_forward_fft();
  //  c.multiply(forward_coefficient);
  multiply_factors(c, FORWARD);

  //and back to the detector plane
  c.invert(true);
  c.perform_forward_fft();

}

void FresnelCDI_WF::set_support(double z_diameter, double size){

  double center_pixel = (nx-1)/2;
  double delta_z = (zone_to_focal_length/focal_to_detector_length)*pixel_length;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){ 
      //calculate the distance from the center
      double i_ = i-center_pixel;
      double j_ = j-center_pixel;
      double r = delta_z*sqrt(i_*i_+j_*j_);
      if(r < size*z_diameter/2.0){
	support.set(i,j,1);
//	support_init.set(i,j,1);

      }
    }
  }
  convolve(support,3,5);
}
