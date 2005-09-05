// Methods for the class List (a dynamic list of some values)
//
// 23/01/2004: start working on it
//
#ifndef LIST_CPP
#define LIST_CPP

#include "List.h"

using namespace std;

// ***** Implicit constructor
// ============================
 template<typename T>
List<T>::List() 
  : first (NULL), 
    last (NULL), 
    size(0)
{
}


// ***** Copy constructor
// ========================
 template<typename T>
List<T>::List (const List<T>& ls)
{
  size = 0;
  first = NULL;
  last = NULL;
  
  Node* temp = ls.first;
  while (temp != NULL){
    Node* added = new Node;
    if (added == NULL){
      returnR error("C++ Error: Unable to allocate a memory.", 99);
      throw error; 
    }
    added->item = temp->item;
    added->next = NULL;
    size++;

    if (first == NULL) first = added;
    else               last->next = added;
    last = added;
    temp = temp->next;
  }
}


// ***** Destructor
// =================
 template<typename T>
List<T>::~List()
{
  Node* temp;
  while (first != NULL){
    temp = first;
    first = first->next;
    delete temp;
  }
}


// ***** Assignment operator
// ==========================
 template<typename T>
List<T>&
List<T>::operator=(const List<T>& ls)
{
    // Protect against itself assignment
  if (this == &ls) return *this;

    // Delete possible old stuff from the left hand side
  Node* temp;
  while (first != NULL){
    temp = first;
    first = first->next;
    delete temp;
  }
  
    // Copy the right hand side to the left hand side
  size = 0;
  first = NULL;
  last = NULL;
  
  temp = ls.first;
  while (temp != NULL){
    Node* added = new Node;
    if (added == NULL){
      returnR error("C++ Error: Unable to allocate a memory.", 99);
      throw error; 
    }
    added->item = temp->item;
    added->next = NULL;
    size++;

    if (first == NULL) first = added;
    else               last->next = added;
    last = added;
    temp = temp->next;
  }

  return *this;
}


// ***** Retrieve the ith element of the List
// ==========================================
 template<typename T>
T& 
List<T>::operator[] (const int i) const
{
  int j;

  if (i >= length() || i < 0){
    returnR error("C++ List Error: Index out of range.", 99);
    throw error; 
  }

  Node* temp = first;
  for (j = 0; j < i; j++){
    temp = temp->next;
  }
  return temp->item;
}


// ***** Return the length of the List
// ====================================
 template<typename T>
int
List<T>::length() const
{
  return size;
}


// ***** Find out whether the List is empty
// =========================================
 template<typename T>
bool
List<T>::isEmpty() const
{
  if (size == 0) return true;
  else           return false;
}


// ***** Add a new node to the end of the List
// =============================================
 template<typename T>
void
List<T>::addNode(const T& element)
{
  Node* added = new Node;
  if (added == NULL){
    returnR error("C++ Error: Unable to allocate a memory.", 99);
    throw error; 
  }
  added->item = element;
  added->next = NULL;
  size++;
  if (first == NULL) first = added;
  else               last->next = added;
  last = added;

  return;
}


// ***** operator << 
// ===================
 template<typename S>
std::ostream&
operator<<(std::ostream& os, const List<S>& ls)
{
  if (ls.size == 0){
    os << "NULL";
    return os;
  }

  typename List<S>::Node* temp = ls.first;
  while (temp != NULL){
    os << temp->item << ",   ";
    temp = temp->next;
  }
  return os;

}


#endif
