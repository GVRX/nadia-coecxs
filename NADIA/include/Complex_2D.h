// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file Complex_2D.h
 * @class Complex_2D
 * @author Nadia Davidson 
 * @date Last Modified on 21/2/2011
 *
 * @brief A 2-dimensional array of complex numbers 
 *
 * This class represents a 2D complex field. Setter and getter methods
 * are provided along with some other useful functions (add,
 * multiplying and fast fourier transforming). Complex_2D objects are
 * used in the CDI reconstruction to represent the ESW in a single
 * plane.
 */

#ifndef COMPLEX_2D_H
#define COMPLEX_2D_H

#include <math.h>
#include <fftw3.h>
#include <Double_2D.h>
#include <types.h>

#ifndef NUM_THREADS
#define NUM_THREADS 1
#endif


template<class T>
class ComplexR_2D{

  /** nx/ny are number of samplings in x/y */
  int nx, ny;

  /* A flag which is passed to fftw when plans are created */
  int fftw_type;

  /** "array" holds the data */
  FFTW_COMPLEX *array;
  /* A fftw plan for forward fourier transforms */
  FFTW_PLAN f_forward;
  /* A fftw plan for backward fourier transforms */
  FFTW_PLAN f_backward;

  int malloc_size;

  
 public:

  /**
   * Constructor that creates a 2D object with the given dimensions.
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   *
   */
  ComplexR_2D(int x_size, int y_size);

  /**
   * Destructor
   *
   */
  ~ComplexR_2D();

  /**
   *Copy Constructor
   */
  ComplexR_2D(const ComplexR_2D& object);

  /**
   *Copy Constructor from a Double_2D
   */
  ComplexR_2D(const Double_2D& object);


  /**
   *Assignment Operator
   */

  ComplexR_2D& operator=(const ComplexR_2D& rhs){

    nx = rhs.get_size_x();
    ny = rhs.get_size_y();

    malloc_size=sizeof(FFTW_COMPLEX)*rhs.get_size_x()*rhs.get_size_y();

    FFTW_FREE(array);
    if(f_forward)
          FFTW_DESTROY_PLAN(f_forward);
    if(f_backward)
          FFTW_DESTROY_PLAN(f_backward);

    array = (FFTW_COMPLEX*) FFTW_MALLOC(malloc_size);

    memcpy(array, rhs.array, malloc_size);


    //initalise the fftw plans to null (not created yet. We will
    //create them when needed to avoid unnecessary time overhead).
    f_forward = 0;
    f_backward = 0;
    fftw_type = FFTW_MEASURE;

/*    for(int i=0; i <rhs.get_size_x(); i++){
      for(int j=0; j < rhs.get_size_y(); j++){
	set_real(i, j, rhs.get_real(i, j));
	set_imag(i, j, rhs.get_imag(i, j));
      }
    }
*/
  };

  /**
   *Assignment Operator to make a Complex_2D with zero imaginary 
   *component from a Double_2D
   */

  ComplexR_2D& operator=(const Double_2D& rhs){

    nx = rhs.get_size_x();
    ny = rhs.get_size_y();

    FFTW_FREE(array);
    if(f_forward)
      FFTW_DESTROY_PLAN(f_forward);
    if(f_backward)
      FFTW_DESTROY_PLAN(f_backward);
 
    malloc_size=sizeof(FFTW_COMPLEX)*rhs.get_size_x()*rhs.get_size_y();
    array=(FFTW_COMPLEX*) FFTW_MALLOC(malloc_size);

    //initalise the fftw plans to null (not created yet. We will
    //create them when needed to avoid unnecessary time overhead).
    f_forward = 0;
    f_backward = 0;
    fftw_type = FFTW_MEASURE;

    for(int i=0; i <rhs.get_size_x(); i++){
      for(int j=0; j < rhs.get_size_y(); j++){
	set_real(i, j, rhs.get(i, j));
	set_imag(i, j, 0);
      }
    }

  };



  /**
   * Set the value at point x,y. Note that this is
   * the slow method. For fast, but unsafe, methods 
   * use set_real or set_imag.
   * 
   * @param x The x position
   * @param y The y position
   * @param type Which component to set. The options are either: 
   * "REAL" or "IMAG"
   * @param value The value which it will be set to
   *  
   */
  void set_value(int x, int y, int type, T value);


  /**
   * Set the real components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @param value The value which it will be set to
   *  
   */
  inline void set_real(int x, int y, T value){
    array[x*ny+y][REAL] = value;
  };

  /**
   * Set the imaginary components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @param value The value which it will be set to
   *  
   */
  inline void set_imag(int x, int y, T value){
    array[x*ny+y][IMAG] = value;
  };

  /**
   * Set the magnitude at point x,y, while preserving the phase. Note
   * that this is an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @param value The value which it will be set to
   *  
   */
  inline void set_mag(int x, int y, T value){
    double mag = get_mag(x,y);
    if(mag==0){
      set_real(x,y,value);
      set_imag(x,y,0);
    }
    else{
      array[x*ny+y][REAL]*=value/mag;
      array[x*ny+y][IMAG]*=value/mag;
    }
  };

  /**
   * Set the phase at point x,y, while preserving the magnitude. Note
   * that this is an unsafe method as no bounds checking is performed.
   * 
   * @param x The x position
   * @param y The y position
   * @param value The value which it will be set to
   *  
   */
  inline void set_phase(int x, int y, T value){
    double mag = get_mag(x,y);
    set_real(x,y,mag*cos(value));
    set_imag(x,y,mag*sin(value));
  };


  /**
   * Get the real components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */
  inline T get_real(int x, int y) const{
    return array[x*ny+y][REAL];
  };

  /**
   * Get the imaginary components at point x,y. Note that this is
   * an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */
  inline T get_imag(int x, int y) const{
    return array[x*ny+y][IMAG];
  };

  /**
   * Get the magnitude at point x,y, @f$ \sqrt{\mathrm{real}^2 + \mathrm{imag}^2} @f$
   * Note that this is an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */
  inline T get_mag(int x, int y) const{
    return sqrt(array[x*ny+y][REAL]*array[x*ny+y][REAL]+
	array[x*ny+y][IMAG]*array[x*ny+y][IMAG]);
  };

  /**
   * Get the phase at point x,y. It goes between -pi and pi.  Please
   * note that this method is ineffectual when the magnitude is
   * 0. Hence, if you are setting values in a Complex_2D array for the
   * first time, you must use set_mag and then set_phase.  Also not,
   * this is an unsafe method as no bounds checking is performed.
   * 
   * @param x The horizontal position
   * @param y The vertical position
   * @return The value at (x,y)  
   */

  inline T get_phase(int x, int y) const{
    T phase = atan2(get_imag(x,y),get_real(x,y));

    if( phase > M_PI )
      return phase - 2*M_PI;

    return phase;

  };


  /**
   * Get the value at point x,y. Note that this is
   * the slow method. For fast, but unsafe, methods 
   * use get_real, get_imag or get_mag.
   * 
   * @param x The x position
   * @param y The y position
   * @param type Which component to set. The options are either: 
   * "REAL","IMAG","MAG","MAG_SQ" or "PHASE"
   * @return The value
   *  
   */
  T get_value(int x, int y, int type) const; 


  /**
   * Get the size in x;
   * 
   * @return The number of horizontal points.
   *  
   */
  int get_size_x() const{
    return nx;
  };

  /**
   * Get the size in y;
   * 
   * @return The number of vertical points.
   *  
   */
  int get_size_y() const{
    return ny;
  };

  /**
   * Get a 2D array of real numbers. 
   * 
   * @param type Which type of value is wanted. The options are either: 
   * "REAL","IMAG","MAG","MAG_SQ" or "PHASE"
   * @param result A 2D array. The array will be filled with the
   * result. Note: This method does not allocated memory for the array,
   * so this should be done before making the call.
   */
  void get_2d(int type, Double_2D & result) const;

  /**
   * Scale the real and imaginary components of the array by a factor. 
   * 
   * @param scale_factor Number which is multiplied by all components
   * of the field.
   */
  void scale(T scale_factor);

  /**
   * Add another Complex_2D to this Complex_2D. The values in this object
   * will be modified.
   * 
   * @param c2 The Complex_2D to add.
   * @param scale If this value is non-empty, or not 1, c2 will be
   * scaled before being added to the Complex_2D,
   * @f $\mathrm{this} = \mathrm{this} +
   * \mathrm{scale} \times \mathrm{c2} $ @f . 
   * Using this function is more efficient than
   * calling Complex_2D::scale() followed by Complex_2D:add()
   * separately.
   */
  void add(ComplexR_2D & c2, T scale=1);


  /**
   * Multiple another Complex_2D to this Complex_2D. The values in this object
   * will be modified.
   * 
   * @param c2 The Complex_2D to add.
   * @param scale If this value is non-empty, or not 1, c2 will be
   * scaled before being added to the Complex_2D,
   * @f $\mathrm{this} = \mathrm{this} \times
   * \mathrm{scale} \times \mathrm{c2} $ @f . 
   * Using this function is more efficient than
   * calling Complex_2D::scale() followed by Complex_2D:multiply()
   * separately.
   */

  void multiply(ComplexR_2D & c2, T scale=1);

  /**
   * Multiple a Double_2D to this Complex_2D. The values in this object
   * will be modified.
   * 
   * @param d2 The Double_2D to multiply.
   * @param scale If this value is non-empty, or not 1, d2 will be
   * scaled before being multiplied with this Complex_2D,
   * @f $\mathrm{this} = \mathrm{this} \times
   * \mathrm{scale} \times \mathrm{c2} $ @f . 
   * Using this function is more efficient than
   * calling Complex_2D::scale() followed by Complex_2D:multiply()
   * separately.
   */
  void multiply(Double_2D & d2, T scale=1);




  /**
   * Get the norm of this Complex_2D: @f $     $ @f.
   * 
   * @return @f $ \sqrt{ \sum_{x,y}{ \mathrm{|C(x,y)|^2} }  } $ @f.
   */
  T get_norm() const;

  /**
   * Create a new Complex_2D with the same values as this one.
   * 
   * @return The new Complex_2D
   */
  ComplexR_2D * clone() const;

  /**
   * Copy the values from another Complex_2D to this one.
   * 
   * @param c The Complex_2D which will be copied from.
   */
  void copy(const ComplexR_2D & c);


  /**
   * Rearrange the matrix so that the corners are placed in the
   * center. This is used after fourier transforming. i.e.
   * if the array is represented by as 4 quadrants: 
   * <br> 1 2
   * <br> 3 4
   * <p> It becomes:
   * <br> 4 3
   * <br> 2 1
   * 
   * This function may also scale the array while it inverts it,
   * required. This is aimed at reducing the CPU time needed, as both
   * tasks are often performed after fourier transforming.
   *
   * @param scale Turn the scaling on-true or off-false. If on, the
   * matrix will be scaled to 1/sqrt(nx*ny) while it is being inverted.
   * By default the array is not scaled.
   */
  void invert(bool scale=false);

  void conjugate();


  /**
   * Forward fourier transform the Complex_2D object. The
   * Complex_2D is not scaled to give the same normalisation as 
   * before the transform. The invert() function can be used to
   * acheive this.
   *
   */
  void perform_forward_fft();

  /**
   * Backward fourier transform the Complex_2D object. The
   * Complex_2D is not scaled to give the same normalisation as 
   * before the transform. The invert() function can be used to
   * acheive this.
   *
   */
  void perform_backward_fft();

  void mirror();


  void perform_backward_fft_real(Double_2D & result);

  void perform_forward_fft_real(Double_2D & input);

  /**
   * A flag which is passed to fftw when plans are created. It maybe
   * either FFTW_ESTIMATE, FFTW_MEASURE or FFTW_PATIENT. See the fftw
   * documentation for a description of each. By default we use
   * FFTW_MEASURE. FFTW_ESTIMATE is used for testing purposes.
   *
   * @param type - FFTW_ESTIMATE, FFTW_MEASURE or FFTW_PATIENT
   */
  void set_fftw_type(int type){
    fftw_type = type;
  };

  /**
   * Create the same complex with some padding. The padding is filled with 0s
   **/

  ComplexR_2D get_padded(int x_add, int y_add);

  /**
   * Create the same complex without the some padding. 
   */
  ComplexR_2D get_unpadded(int x_add, int y_add);



 private:

  /**
   * Check that an (x,y) position is within the bounds of the array
   * 
   * @param x The horizontal position to check
   * @param y The vertical position to check
   */
  int check_bounds(int x, int y) const;


  /**
   * Create the fftw plans for forward and backward fast fourier transforms
   */
  void initialise_fft();


};
//////////////////////////////////////////
#ifndef DOUBLE_PRECISION
typedef ComplexR_2D<float> Complex_2D;
#else
typedef ComplexR_2D<double> Complex_2D;
#endif //DOUBLE_PRECISION
//////////////////////////////////////////

#endif
