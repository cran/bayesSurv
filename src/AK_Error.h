/* AK_Error.h
 *
 * This header provedes class definition for errors thrown by my functions
 * such that thay can be catched and solved without aborting R.
 */


#ifndef AK_ERROR_H
#define AK_ERROR_H

#include <exception>
#include <string>
#include <sstream>
#include <iostream>


  /**** Class which in combination with catch in user's functions provides a safe return to R.
   *    Its constructor outputs the string given to it as its parameter 
   *    thereby notifying the user of why the program crashed.            
  /* ======================================================================================= */
  class returnR
  {
    public:
      returnR (std::string & sterr, const int err)
               : fout_ (err)
      {
        std::cerr << sterr << std::endl;
        std::cerr << std::endl;
      }

      returnR (const char * sterr, const int err)
	: fout_ (err)
      { 
	std::cerr << sterr << std::endl;
	std::cerr << std::endl;
      }

      returnR (const int err)
               : fout_ (err)
      {
        fout_ = 99;
      }

      returnR & operator= (const returnR &rr) 
      {
        fout_ = rr.fout_;
        return *this;
      }

      ~returnR () 
      {
      }

      inline int errflag ()
      {
        return(fout_);
      }

    private:
      int fout_;
  };  // end of class returnR

#endif



