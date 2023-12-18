#ifndef TEMPLATE_FUN_H
#define TEMPLATE_FUN_H

#include <R.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AK_Error.h"
#include "openFile.h"

template <typename dd>
void
cumsum(dd* csum,   const dd* vals,   const int* kP);

//template <typename pp>                                    // commented on 14/12/2023
//void
//printArray(const pp* array, const int* length);
//
// /** explicit printArrayI() and printArrayD() moved to printArray.[h,cpp], 14/12/2023 **/
  
template <typename pp>
void
changePointers(pp** aP, pp** bP);

template <typename dd>
void
writeToFile(const dd* array,          const int& nR,                const int& nC,
            const std::string& dir,   const std::string& filename,  const char &flag = 'n',
            const int& prec = 6,      const int& width = 0);

template <typename dd>
void
writeAddToFile(const dd* array,          const int& nR,                const int& nC,           const dd& add,
               const std::string& dir,   const std::string& filename,  const char& flag = 'n',
               const int& prec = 6,      const int& width = 0);

template <typename dd1, typename dd2>
void
writeTwoToFile(const dd1* array1,       const int& nR1,               const int& nC1,    const int& col1,
               const dd2* array2,       const int& nR2,               const int& nC2,
               const std::string& dir,  const std::string& filename,  const char &flag,  
               const int& prec = 6,     const int& width = 0);

template <typename dd>
void
writeToFile2(dd** array,               const int n1,                 const int n2,
             const std::string& dir,   const std::string& filename,  const char &flag = 'n',
             const int& prec = 6,      const int& width = 0);

template <typename dd>
void
writeRaggedToFile(const dd* array,         const int& nR,                const int& maxnC,  
                  const int* nC,           const int& multnC,
                  const std::string& dir,  const std::string& filename,  const char &flag = 'n',
                  const int& prec = 6,     const int& width = 0);

template <typename dd>
void
readFromFile(dd* array,               int* nread,       
             const int& nR,           const int& nC,               const int& header,
             const int& skip,         const int& by,
             const std::string& dir,  const std::string& filename,  
             const int& skipOnRow);

#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
// Necessary for template instantiation with some compilers.
# include "templatefun.cpp"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
