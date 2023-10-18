# Programming-Phase-field-in-Fortran

The repository shows the use of Fortran programming language for the phase-field methods &mdash; model A, model B, and model C. The codes are 2D. Model A ( Allen-Cahn equation ), and Model B ( Cahn-Hilliard ) examples are given; solidification of a single component alloy is also provided.

![Output](images/modelAB.png)
___
## **Significance** 

 * <span style="color:blue"> **Sequential code**</span>.

 * Self-explanatory <span style="color:blue"> **comments**</span> for each section.

 * <span style="color:blue"> **Stand-alone**</span> codes.

 * Code description in <span style="color:blue"> **README**</span> files.

 * <span style="color:blue"> **Examples:**</span> FeCr alloy for model B; grain evolution for model A. 

 * Some of the <span style="color:blue"> **best programming practices**</span> .

      * *if statement* rather than *if then* construct

 * Demonstration of <span style="color:blue"> **concurrent programming technique**</span> .

      * Model B example

 * File savings with <span style="color:blue"> **logical unit number**</span>.

      * Model A example

 * <span style="color:blue"> **Screen shots**</span> for each output code.

 * <span style="color:blue"> **Data availability**</span> for each code.

 * <span style="color:blue"> **Compiler comparison**</span>  for further studies.

      * Model B README

 * Instructions for <span style="color:blue"> **different operating systems**</span>.

 * <span style="color:blue"> **Optimization**</span>  options.

 * Integrated <span style="color:blue"> **dislin graphical software**</span> in each code:

    * Single statement <span style="color:blue"> **quick plots**</span>.

    * Routines for <span style="color:blue"> **continuous animation**</span>. 

    * <span style="color:blue"> **Multiplot**</span> techniques.

 * <span style="color:blue"> **gnuplot script**</span> for output files:
    * Customized <span style="color:blue"> **color plots**</span>.

    * commands for <span style="color:blue"> **continuous animation**</span>. 

    * <span style="color:blue"> **Multiplot**</span> approach.

## **Computational Tools**

The simulations were performed on the system with the following details:

|                OS      |      Compilers and versions               |  Integrated graphics library  |  Output graphics library   |
| -----------------------| ----------------------------------------- |------------------------------ |----------------------------|
| Linux (ubuntu 20.04)   |     gfortran (12.1.0)                     |  Dislin ( 11.5 )              |     gnuplot ( 5.4 )        |
| Windows (10, 64 bit)   |     gfortran (12.1.0), intel (2021.6.0)   |  Dislin ( 11.5 )              |     gnuplot ( 5.4 )        |


## **Future work**

* [Whole array programming](https://github.com/Shahid718/Phase-field-Fortran-codes-using-whole-array)
  
*  Procedural programming

      * [internal procedures](https://github.com/Shahid718/Fortran-Phase-field-codes-using-Internal-Procedures)
      * external procedures
               
* Modular programming

* Object-oriented programming

* Parallel programming with 

  * Co-arrays
  * OpenMP
  * MPI
    
* GPU programming with 

  * CUDA Fortran
  * OpenACC
    
## **Contact**

shahid.maqbool@rwth-aachen.de

shahid@postech.ac.kr
