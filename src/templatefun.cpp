// Template functions
//
// 01/02/2004: 'printArray'
//             'changePointers'
//             'writeToFile'
//             'writeAddToFile'
//             'writeTwoToFile'
//             'writeToFile2'
//             'writeRaggedToFile'
//             'readFromFile'
// 08/03/2004: 'cumsum'
// 19/06/2012: std::cout replaced by Rprintf
//
#ifndef TEMPLATE_FUN_CPP
#define TEMPLATE_FUN_CPP

#include "templatefun.h"

// using namespace std;

// ====================================================
// Function to compute cumulative sums of an array
// ====================================================
//
// csum ...... an array of cumulative sums
// vals ...... original values
// kP ........ length of arrays csum and vals
//
 template <typename dd>
void
cumsum(dd* csum,   const dd* vals,   const int* kP)
{
  int i;
  csum[0] = vals[0];
  for (i = 1; i < *kP; i++){
    csum[i] = csum[i - 1] + vals[i];
  }
  
  return;
}

// ================================================================
// Function to print an array on a screen
//  -> useful especially for debugging
// ================================================================
 template <typename pp>
void
printArray(const pp* array, const int* length)
{
  //for (int i = 0; i < *length; i++) std::cout << array[i] << ",  ";
  //std::cout << std::endl;
  for (int i = 0; i < *length; i++) Rprintf("%g,  ", array[i]);
  Rprintf("\n");
  return;
}


// ================================================================================
// Function to change two pointers
// ================================================================================
 template <typename pp>
void
changePointers(pp** aP, pp** bP)
{
  pp* helpP;
  helpP = *aP;
  *aP = *bP;
  *bP = helpP;
  return;
}


// ================================================================================
// ***** writeToFile: Function to write a numeric array to a file             *****
// ================================================================================
// The content of the array is written to the file in ROW major order.
//
// nR .......... number of rows into which the array should be splitted
// nC .......... number of columns into which the array should be splitted 
// flag ........ it can be 'a': append
//                         'o': owerwrite
//                         'n': don't replace
//
 template <typename dd>
void
writeToFile(const dd* array,          const int& nR,                const int& nC,
            const std::string& dir,   const std::string& filename,  const char& flag,
            const int& prec,          const int& width)
{
  try{
    std::string path = dir + filename;
    std::ofstream out;
    openFile(out, path, flag);

    std::ostringstream s;
    unsigned int mlen = width;

    /* Passes up to 5 rows to get things to line up nicely */
    for (int i = 0; i < nR && i < 5; i++) {
      for (int j = 0; j < nC; j++) {
	s.str("");        
        if (array[i*nC + j] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (array[i*nC + j] < 1 && array[i*nC + j] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << array[i*nC + j] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << array[i*nC + j] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

//    s.str("");
    for (int i = 0; i < nR; i++) {
      for (int j = 0; j < nC; j++){
        if (array[i*nC + j] >= FLT_MAX){                                                 
          out << std::setw(mlen) << "1e50";
          out << "   ";
	}
        else{
          if (array[i*nC + j] < 1 && array[i*nC + j] > -1){
            out << std::scientific << std::setw(mlen) << std::setprecision(prec) << array[i*nC + j];
            out << "   ";
          }
          else{
            out << std::fixed << std::setw(mlen) << std::setprecision(prec) << array[i*nC + j];
            out << "   ";
          }
        }
      }
      out << std::endl;
    }
//    out << s.str();
    out.close();
    return;
  }  // end of try
  catch(returnR){
    throw;
  }  
}   // end of function writeToFile


// ================================================================================
// ***** writeAddToFile: Function to write a numeric array to a file          *****
//         - some value is added to all elements of the array before
//           writting them to the file
//         - useful when recording C++/R indeces
// ================================================================================
// The content of the array is written to the file in ROW major order.
//
// nR .......... number of rows into which the array should be splitted
// nC .......... number of columns into which the array should be splitted 
// flag ........ it can be 'a': append
//                         'o': owerwrite
//                         'n': don't replace
//
 template <typename dd>
void
writeAddToFile(const dd* array,          const int& nR,                const int& nC,     const dd& add,
               const std::string& dir,   const std::string& filename,  const char& flag,
               const int& prec,          const int& width)
{
  try{
    std::string path = dir + filename;
    std::ofstream out;
    openFile(out, path, flag);

    std::ostringstream s;
    unsigned int mlen = width;

    /* Passes up to 5 rows to get things to line up nicely */
    for (int i = 0; i < nR && i < 5; i++) {
      for (int j = 0; j < nC; j++) {
	s.str("");        
        if (array[i*nC + j] + add >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (array[i*nC + j] + add < 1 && array[i*nC + j] + add > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << array[i*nC + j] + add << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << array[i*nC + j] + add << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

//    s.str("");
    for (int i = 0; i < nR; i++) {
      for (int j = 0; j < nC; j++){
        if (array[i*nC + j] + add >= FLT_MAX){
          out << std::setw(mlen) << "1e50";
          out << "   ";
	}
        else{
          if (array[i*nC + j] + add < 1 && array[i*nC + j] + add > -1){
            out << std::scientific << std::setw(mlen) << std::setprecision(prec) << array[i*nC + j] + add;
            out << "   ";
          }
          else{
            out << std::fixed << std::setw(mlen) << std::setprecision(prec) << array[i*nC + j] + add;
            out << "   ";
          }
        }
      }
      out << std::endl;
    }
//    out << s.str();
    out.close();
    return;
  }  // end of try
  catch(returnR){
    throw;
  }  
}   // end of function writeAddToFile


// ================================================================================
// Function to write part of one array and whole second array to a file
// ================================================================================
// * used mainly in function 'writeToFiles' to write number of mixture components,
//   mixture intercept and mixture scale into one file
//
// nR1 .......... number of rows into which the first array should be splitted
// nC1 .......... number of columns into which the first array should be splitted 
// col1 ......... index of the column from the first array that is to be written
// nR2 .......... number of rows into which the second array should be splitted
// nC2 .......... number of columns into which the second array should be splitted 
// flag ......... it can be 'a': append
//                          'o': owerwrite
//                          'n': don't replace
//
 template <typename dd1, typename dd2>
void
writeTwoToFile(const dd1* array1,        const int& nR1,                 const int& nC1,      const int& col1,
               const dd2* array2,        const int& nR2,                 const int& nC2,
               const std::string& dir,   const std::string& filename,    const char &flag,
               const int& prec,          const int& width)
{
  try{
    if (nR1 != nR2) throw returnR("C++ programming error: contact the author", 99);
    std::string path = dir + filename;
    std::ofstream out;
    openFile(out, path, flag); 

    std::ostringstream s;
    unsigned int mlen = width;

    /* Passes up to 5 rows of the second array to get things to line up nicely */
    for (int i = 0; i < nR2 && i < 5; i++) {
      for (int j = 0; j < nC2; j++) {
	s.str("");        
        if (array2[i*nC2 + j] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (array2[i*nC2 + j] < 1 && array2[i*nC2 + j] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << array2[i*nC2 + j] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << array2[i*nC2 + j] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }


    /* Write to files */
//    s.str("");
    for (int i = 0; i < nR1; i++) {
      if (array1[i*nC1 + col1] >= FLT_MAX){                                                 
        out << std::setw(mlen) << "1e50";
        out << "   ";
      }
      else{
        if (array1[i*nC1 + col1] < 1 && array1[i*nC1 + col1] > -1){
          out << std::scientific << std::setw(mlen) << std::setprecision(prec) << array1[i*nC1 + col1];
          out << "   ";
        }
        else{
          out << std::fixed << std::setw(mlen) << std::setprecision(prec) << array1[i*nC1 + col1];
          out << "   ";
        }
      }
      for (int j = 0; j < nC2; j++){
        if (array2[i*nC2 + j] >= FLT_MAX){                                                 
          out << std::setw(mlen) << "1e50";
          out << "   ";
	}
        else{
          if (array2[i*nC2 + j] < 1 && array2[i*nC2 + j] > -1){
            out << std::scientific << std::setw(mlen) << std::setprecision(prec) << array2[i*nC2 + j];
            out << "   ";
          }
          else{
            out << std::fixed << std::setw(mlen) << std::setprecision(prec) << array2[i*nC2 + j];
            out << "   ";
          }
        }
      }
      out << std::endl;
    }
//    out << s.str();
    out.close();
    return;
  }  // end of try
  catch(returnR){
    throw;
  }  
}   // end of function writeTwoToFile


// ================================================================================
// Function to write a numeric array (double indexed) to a file
// ================================================================================
// The first index gives columns, the second index rows in the resulting file
//
// dd .......... double indexed array
// n1 .......... range of the first index -> columns in the resulting file
// n2 .......... range of the second index -> rows in the resulting file
// flag ........ it can be 'a': append
//                         'o': owerwrite
//                         'n': don't replace
//
 template <typename dd>
void
writeToFile2(dd** array,               const int n1,                 const int n2,
             const std::string& dir,   const std::string& filename,  const char &flag,
             const int& prec,          const int& width)
{
  try{
    std::string path = dir + filename;
    std::ofstream out;
    openFile(out, path, flag); 

    std::ostringstream s;
    unsigned int mlen = width;

    /* Passes up to 5 rows to get things to line up nicely */
    for (int i = 0; i < n2 && i < 5; i++) {
      for (int j = 0; j < n1; j++) {
	s.str("");        
        if (array[j][i] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (array[j][i] < 1 && array[j][i] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << array[j][i] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << array[j][i] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

//    s.str("");
    for (int i = 0; i < n2; i++) {
      for (int j = 0; j < n1; j++){
        if (array[j][i] >= FLT_MAX){                                                 
          out << std::setw(mlen) << "1e50";
          out << "   ";
	}
        else{
          if (array[j][i] < 1 && array[j][i] > -1){
            out << std::scientific << std::setw(mlen) << std::setprecision(prec) << array[j][i];
            out << "   ";
          }
          else{
            out << std::fixed << std::setw(mlen) << std::setprecision(prec) << array[j][i];
            out << "   ";
          }
        }
      }
      out << std::endl;
    }
//    out << s.str();
    out.close();
    return;
  }  // end of try
  catch(returnR){
    throw;
  }  
}   // end of function writeToFile2


// ===========================================
// Function to write a ragged array to file
// ===========================================
// 
// array[nR*maxnC] ..... array to be written to a file
//                       * values stored in ROW MAJOR order
//                       * this array is assumed to be regular (not ragged),
//                         however, maxnC - multnC * nC[i] from the i-th row will be ignored for writting
// nR .................. number of rows into which the array is to be splitted
// maxnC ............... maximal number of columns (i.e. nC[j]*multnC <= maxnC for all j)
// nC[nR]............... number of columns on each row that should be after multiplying by multnC written to the file
// multnC .............. each nC is multiplied by multnC to get number of columns that should be written to the file
//                       (useful for writing bivariate (multnc = 2) means, e.g. in bayesHistogram)
//
 template <typename dd>
void
writeRaggedToFile(const dd* array,         const int& nR,                const int& maxnC,  
                  const int* nC,           const int& multnC,
                  const std::string& dir,  const std::string& filename,  const char &flag,
                  const int& prec,         const int& width)
{
  try{
    int i, j;
    std::string path = dir + filename;
    std::ofstream out;
    openFile(out, path, flag);

    /*** Write to the file ***/
    std::ostringstream s;
    unsigned int mlen = width;
   
    /* Passes up to 5 rows to get things to line up nicely */
    for (i = 0; i < nR && i < 5; i++) {
      if (multnC * nC[i] > maxnC) throw returnR("C++ Error: multnC * nC must be <= maxnC in writeRaggedToFile", 1);
      for (j = 0; j < multnC * nC[i]; j++) {
	s.str("");        
        if (array[i*maxnC + j] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (array[i*maxnC + j] < 1 && array[i*maxnC + j] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << array[i*maxnC + j] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << array[i*maxnC + j] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

    /* Write */
//    s.str("");
    for (i = 0; i < nR; i++) {
      if (multnC * nC[i] > maxnC) throw returnR("C++ Error: multnC * nC must be <= maxnC in writeRaggedToFile", 1);
      for (j = 0; j < multnC * nC[i]; j++){
        if (array[i*maxnC + j] >= FLT_MAX){
          out << std::setw(mlen) << "1e50";
          out << "   ";
	}
        else{
          if (array[i*maxnC + j] < 1 && array[i*maxnC + j] > -1){
            out << std::scientific << std::setw(mlen) << std::setprecision(prec) << array[i*maxnC + j];
            out << "   ";
          }
          else{
            out << std::fixed << std::setw(mlen) << std::setprecision(prec) << array[i*maxnC + j];
            out << "   ";
          }
        }
      }
      out << std::endl;
    }
//    out << s.str();
    out.close();
    return;  
  }
  catch(returnR){
    throw;
  }  
}  /** end of function writeRaggedToFile **/



// =========================================
// Function to read an array from a file
// =========================================
// The file is read by ROWS and written to an array
//
// array........ array where to read it        (length = nread * nC)
// nread ....... on OUTPUT: number of read rows
// nR .......... number of rows in file THAT ARE AVAILABLE for reading after skipping possible header
//               (this will usually be MCMC sample size)
// nC .......... number of columns in file THAT ARE TO BE READ
// header ...... 0/1 is there header?
// skip ........ how many rows at the beginning of the file are to be skipped (do not count header)
// by .......... only every by-th mixture will be read
// dir ......... directory
// filename .... name of the file
// skipOnRow ... number of values that are to be skipped at the beginning of each row
// skip ........ number of rows thar are to be skipped at the beginnig of the file
//
 template <typename dd>
void
readFromFile(dd* array,               int* nread,       
             const int& nR,           const int& nC,               const int& header,
             const int& skip,         const int& by,
             const std::string& dir,  const std::string& filename,  
             const int& skipOnRow)
{
  try{
    int i, j, ii;
    int size = nR * nC;
    if (size <= 0) throw returnR("C++ Error: File of null size is to be read.", 99);
    if (skip < 0) throw returnR("C++ Error: 'skip' parameter must be >= 0 in 'readFromFile'", 1);
    if (by <= 0) throw returnR("C++ Error: 'by' parameter must be > 0 in 'readFromFile'", 1);
    if (skip >= nR) throw returnR("C++ Error: too many rows are to be skipped by 'readFromFile'", 1);

    std::string path = dir + filename;

    std::string errmes, mess;
    char cmess[200];
    std::ifstream file(path.c_str(), std::ios::in);    

    dd temp;
    if (!file){
      errmes = std::string("C++ Error: Could not open ") + path;
      throw returnR(errmes, 99);
    } 
    else{
      mess = std::string("Reading ") + path + "\n";
      strcpy(cmess, mess.c_str());
      Rprintf(cmess);

      /*** Skip what is to be skipped (header included) ***/
      char ch;
      for (i = 0; i < skip + header; i++){
        file.get(ch);        
        while (ch != '\n') file.get(ch);
      }

      /*** Read the first row to be read ***/
      *nread = 1;
      double* veld = array;
      if (file.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + path + std::string(" before ") 
                 + char(*nread) + std::string(" rows were read.");
        throw returnR(errmes, 99);
      }
      for (j = 0; j < skipOnRow; j++){
        if (file.eof()){
          errmes = std::string("C++ Error: Reached end of file ") + path + std::string(" before ") 
                   + char(*nread) + std::string(" rows were read.");
          throw returnR(errmes, 99);
        }
        file >> temp;
      }
      for (j = skipOnRow; j < nC + skipOnRow; j++){
        if (file.eof()){
          errmes = std::string("C++ Error: Reached end of file ") + path + " before "
                   + char(*nread) + std::string(" rows were read.");
          throw returnR(errmes, 99);
        }
        file >> (*veld);
        veld++;
      }

      /*** Read remaining rows to be read ***/
      for (i = skip + 1 + by; i <= nR; i += by){

        /** Skip by-1 rows **/
        for (ii = 0; ii < by - 1; ii++){
          file.get(ch);        
          while (ch != '\n') file.get(ch);
        }

        /** Read the values **/
        (*nread)++;
        for (j = 0; j < skipOnRow; j++){
          if (file.eof()){
            errmes = std::string("C++ Error: Reached end of file ") + path + std::string(" before ") 
                     + char(*nread) + std::string(" rows were read.");
            throw returnR(errmes, 99);
          }
          file >> temp;
        }
        for (j = skipOnRow; j < nC + skipOnRow; j++){
          if (file.eof()){
            errmes = std::string("C++ Error: Reached end of file ") + path + " before "
                     + char(*nread) + std::string(" rows were read.");
            throw returnR(errmes, 99);
          }
          file >> (*veld);
          veld++;
        }
        file.get(ch);                 
        while (ch != '\n') file.get(ch);
      }
    }
    file.close();
    return;
  }  // end of try
  catch(returnR){
    throw;
  }  
}    // end of the function readFromFile


#endif
