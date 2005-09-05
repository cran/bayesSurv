#ifndef _OPEN_FILE_H
#define _OPEN_FILE_H

#include <R.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "AK_Error.h"

void
openFile(std::ofstream& ofile,  const std::string& path,  const char& write_flag);

#endif
