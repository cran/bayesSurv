// Template functions for G-spline smoothing
//
// 03/11/2004: 'writeToFile_1'
//             'writeThreeToFile_1'
//             'writeTwoToFile_1'
//             'writeAddToFile_1'
//
#ifndef TEMPLATE_GS_FUN_CPP
#define TEMPLATE_GS_FUN_CPP

#include "templatefun_GS.h"

using namespace std;

// =============================================================================================
// ***** writeToFile_1: Function to write a numeric array as one row to a file             *****
// =============================================================================================
 template <typename dd>
void
writeToFile_1(const dd* array,  const int& nC,     std::ofstream& ofile,
              const int& prec,  const int& width)
{
  for (int j = 0; j < nC; j++){
    if (array[j] >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array[j] < 1 && array[j] > -1 && array[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array[j];
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array[j];
        ofile << "   ";
      }
    }
  }
  ofile << std::endl;

  return;
}


// =================================================================================================
// ***** writeFiveToFile_1: Function to write five numeric arrays as one row to a file         *****
//     this is actually 5x copy of writeToFile_1 function
// =================================================================================================
 template <typename dd>
void
writeFiveToFile_1(const dd* array1,      const dd* array2,  const dd* array3,  const dd* array4,  const dd* array5,
                  const int& nC1,        const int& nC2,    const int& nC3,    const int& nC4,    const int& nC5,
                  std::ofstream& ofile,  const int& prec,   const int& width)
{
  int j;

  for (j = 0; j < nC1; j++){
    if (array1[j] >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array1[j] < 1 && array1[j] > -1 && array1[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array1[j];
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array1[j];
        ofile << "   ";
      }
    }
  }

  for (j = 0; j < nC2; j++){
    if (array2[j] >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array2[j] < 1 && array2[j] > -1 && array2[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array2[j];
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array2[j];
        ofile << "   ";
      }
    }
  }

  for (j = 0; j < nC3; j++){
    if (array3[j] >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array3[j] < 1 && array3[j] > -1 && array3[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array3[j];
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array3[j];
        ofile << "   ";
      }
    }
  }

  for (j = 0; j < nC4; j++){
    if (array4[j] >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array4[j] < 1 && array4[j] > -1 && array4[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array4[j];
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array4[j];
        ofile << "   ";
      }
    }
  }

  for (j = 0; j < nC5; j++){
    if (array5[j] >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array5[j] < 1 && array5[j] > -1 && array5[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array5[j];
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array5[j];
        ofile << "   ";
      }
    }
  }

  ofile << std::endl;

  return;
}


// ========================================================================================================
// ***** writeTwoToFile_1: Function to write part of one array and whole second array as one row in a file
// ========================================================================================================
//
// col1 ......... index of the column from the first array that is to be written
// nC2 .......... number of columns into which the second array should be splitted 
//
 template <typename dd1, typename dd2>
void
writeTwoToFile_1(const dd1* array1,     const dd2* array2,     
                 const int& col1,       const int& nC2,
                 std::ofstream& ofile,  const int& prec,  const int& width)
{
      if (array1[col1] >= FLT_MAX){                                                 
        ofile << std::setw(width) << "1e50";
        ofile << "   ";
      }
      else{
        if (array1[col1] < 1 && array1[col1] > -1 && array1[col1] != 0){
          ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array1[col1];
          ofile << "   ";
        }
        else{
          ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array1[col1];
          ofile << "   ";
        }
      }

      for (int j = 0; j < nC2; j++){
        if (array2[j] >= FLT_MAX){                                                 
          ofile << std::setw(width) << "1e50";
          ofile << "   ";
	}
        else{
          if (array2[j] < 1 && array2[j] > -1 && array2[j] != 0){
            ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array2[j];
            ofile << "   ";
          }
          else{
            ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array2[j];
            ofile << "   ";
          }
        }
      }
      ofile << std::endl;


    return;
}   /** end of function writeTwoToFile_1 **/

// =============================================================================================
// ***** writeAddToFile_1: Function to write a numeric array as one row in a file          *****
//         - some value is added to all elements of the array before
//           writing them to the file
//         - useful when recording C++/R indeces
// =============================================================================================
//
// nR .......... number of rows into which the array should be splitted
//
 template <typename dd>
void
writeAddToFile_1(const dd* array,       const int& nC,    const dd& add,
                 std::ofstream& ofile,  const int& prec,  const int& width)
{
  for (int j = 0; j < nC; j++){
    if (array[j] + add >= FLT_MAX){
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (array[j] + add < 1 && array[j] + add > -1 && array[j] != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << array[j] + add;
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << array[j] + add;
        ofile << "   ";
      }
    }
  }
  ofile << std::endl;

  return;
}   /** end of function writeAddToFile_1 **/


#endif
