#ifndef TEMPLATE_GS_FUN_H
#define TEMPLATE_GS_FUN_H

#include <R.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "AK_Error.h"

template <typename dd>
void
writeToFile_1(const dd* array,      const int& nC,        std::ofstream& ofile,
              const int& prec = 6,  const int& width = 1);

template <typename dd>
void
writeFiveToFile_1(const dd* array1,      const dd* array2,     const dd* array3,     const dd* array4,  const dd* array5,
                  const int& nC1,        const int& nC2,       const int& nC3,       const int& nC4,    const int& nC5,
                  std::ofstream& ofile,  const int& prec = 6,  const int& width = 1);

template <typename dd1, typename dd2>
void
writeTwoToFile_1(const dd1* array1,     const dd2* array2,     
                 const int& col1,       const int& nC2,
                 std::ofstream& ofile,  const int& prec = 6,  const int& width = 1);

template <typename dd>
void
writeAddToFile_1(const dd* array,       const int& nC,        const dd& add,
                 std::ofstream& ofile,  const int& prec = 6,  const int& width = 1);


#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
// Necessary for template instantiation with some compilers.
# include "templatefun_GS.cpp"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
