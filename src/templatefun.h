#ifndef TEMPLATE_FUN_H
#define TEMPLATE_FUN_H

#include <R.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <climits>
#include <string>
#include "AK_Error.h"

template <typename dd>
void
cumsum(dd* csum,   const dd* vals,   const int* kP);

template <typename pp>
void
printArray(const pp* array, const int* length);

template <typename pp>
void
changePointers(pp** aP, pp** bP);

template <typename dd>
void
writeToFile(const dd* array,          const int nR,                 const int nC,
            const std::string& dir,   const std::string& filename,  const char &flag = 'n',
            const int& prec = 6,      const int& width = 0);

template <typename dd>
void
writeTwoToFile(const dd* array1,         const int nR1,                 const int nC1,      const int col1,
               const dd* array2,         const int nR2,                 const int nC2,
               const std::string& dir,   const std::string& filename,   const char &flag,
               const int& prec = 6,      const int& width = 0);

template <typename dd>
void
writeToFile2(dd** array,               const int n1,                 const int n2,
             const std::string& dir,   const std::string& filename,  const char &flag = 'n',
             const int& prec = 6,      const int& width = 0);

template <typename dd>
void
readFromFile(dd* array,               const int nR,                 const int nC,
             const std::string& dir,  const std::string& filename,  
             const int skipOnRow,     const int skip);


#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
// Necessary for template instantiation with some compilers.
# include "templatefun.cpp"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
