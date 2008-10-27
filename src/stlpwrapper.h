// -*- C++ -*-
// $Id: stlpwrapper.h,v 1.6 2006/10/23 09:45:34 hs Exp $


/** @file stlpwrapper.h
 *  @brief simple pointer wrapper
 *
 *  Contains a pointer wrapper derived from Mhumit Khan stl.html
 *  @todo implement a proper smart pointer
 */


// the user still has to delete the objects from the heap

#ifndef STLPWRAPPER_H
#define STLPWRAPPER_H

/** @brief A simple pointer wrapper
 *
 *  Wraps a pointer, so that it can be stored in an STL container.
 *  The current implementation has however little advantages. 
 *  Care has to be taken, for the destructor does not delete dynamic allocated objects!
 *  Users still have to delete these objects from the heap.
 */
template <class TYPE>
class PtrWrapper {
  private:
    /// The pointer
    TYPE* m_pObj;   
  public:
    /// Constructor with a pointer to the TYPE, default is NULL pointer
    PtrWrapper(TYPE* pObj = 0) : m_pObj(pObj) {};
    /// Copy constructor
    PtrWrapper(const PtrWrapper<TYPE>& wrapper) : m_pObj(wrapper.m_pObj) {}
    
    /// Assignment operator
    PtrWrapper<TYPE>& operator= (const PtrWrapper<TYPE>& wrapper) {
      m_pObj = wrapper.m_pObj;
      return *this;
    }
    
    /// Destructor doesn't free pointer
    ~PtrWrapper() {};
    
    /// Typecast to the original pointer type
    operator const TYPE* () const { return m_pObj; }
    /// Typecast to the original pointer type
    operator TYPE* () { return m_pObj; }
    /// Typecast to the original pointer type
    TYPE* operator->() { return m_pObj; }	
    /// Explicit method for the typecast
    TYPE* pObj () const { return m_pObj; };

    /// Compares the two objects by address
    friend bool operator== (const PtrWrapper<TYPE>& left, const PtrWrapper<TYPE>& right) {
      return (left.m_pObj == right.m_pObj);
    };

    /// Compares the two objects
    friend bool operator!= (const PtrWrapper<TYPE>& left, const PtrWrapper<TYPE>& right) {
      return (left.m_pObj != right.m_pObj);
    };

};

#endif // STLPWRAPPER_H
