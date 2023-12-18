#ifndef PRINT_ARRAY_CPP
#define PRINT_ARRAY_CPP

#include "printArray.h"

void                                                                          // added on 14/12/2023
printArrayI(const int* array, const int* length)
{
  for (int i = 0; i < *length; i++) Rprintf("%d,  ", array[i]);
  Rprintf("\n");
  return;
}

void                                                                         // added on 14/12/2023
printArrayD(const double* array, const int* length)
{
  for (int i = 0; i < *length; i++) Rprintf("%g,  ", array[i]);
  Rprintf("\n");
  return;
}

#endif
