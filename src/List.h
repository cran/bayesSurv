// List.h
//
// Class to store a dynamic list of values
//
#ifndef LIST_H
#define LIST_H

#include <iostream>
#include "AK_Error.h"

template <typename T>
class List
{
  struct Node{
    T item;
    Node* next;
  };               // structure to store elements of the list

  private:
    Node* first; 
    Node* last;
    int size;

  public:
    List();

    List(const List<T>& ls);

    ~List();

    List<T>& 
    operator=(const List<T>& ls);

    T& 
    operator[] (const int i) const;

    int 
    length() const;

    bool 
    isEmpty() const;
    
    void
    addNode(const T& element);

    template <typename S>
    friend std::ostream&
    operator<<(std::ostream& os, const List<S>& ls);
    

};    // end of the class List


#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
// Necessary for template instantiation with some compilers.
# include "List.cpp"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
