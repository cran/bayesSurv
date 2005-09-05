#include "openFile.h"


// ====================================================================================
// ***** openFile: open file for output                                          *****/
// ====================================================================================
void
openFile(std::ofstream& ofile,  const std::string& path,  const char& write_flag)
{
  try{
    bool err = false;

    std::string errmess;
    if (write_flag == 'n') {
      std::fstream temp(path.c_str(), std::ios::in);
      if (!temp) {
	ofile.open(path.c_str(), std::ios::out);
      } 
      else {
	temp.close();
	err = true;
      }
    }
    else if (write_flag == 'o')
           ofile.open(path.c_str(), std::ios::out | std::ios::trunc);
         else if (write_flag == 'a')
                ofile.open(path.c_str(), std::ios::out | std::ios::app);
              else {
                errmess = std::string("C++ Error: Incorrect flag for writing to ") + path + ". ";
		returnR error(errmess, 99);
                throw error;
              }

    if (!ofile || err) {
      errmess = std::string("C++ Error: Could not open ") + path + " for writing. ";
      returnR error(errmess, 99);
      throw error; 
    }
    return;
  }
  catch(returnR){
    throw;
  }  
}

