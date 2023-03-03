#pragma once

#include <memory>
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>


// Design problem:
// - Once Rcpp creates an object in an XPtr, it permanently takes
//   ownership and *will* eventually call the destructor
// - Since MatrixLoaders and FragmentLoaders now take ownership of their inputs via unique_ptr,
//   they will also eventually call the destructor for their inputs
// - This could result in double-free, where both R and C++ try to free an input loader
// Solution:
// - Rather than using XPtr<T>, use XPtr<unique_ptr<T>>. This works because unique_ptr is capable of
//   giving up ownership of an object. So then when R calls the destructor on the unique_ptr, the unique
//   ptr will only call the destructor on contents that it owns.
// - This pair of functions allows creating an object and giving ownership to R, then taking back ownership
//   later when passed the object from R.

// Analog of make_unique, but wrapped in an XPtr
template<class T, class... Args>
SEXP make_unique_xptr(Args&&... args) {
  return Rcpp::wrap(Rcpp::XPtr<std::unique_ptr<T>>(new std::unique_ptr<T>(new T(std::forward<Args>(args)...))));
}

// Take ownership of the unique_ptr from the XPtr
template<class T>
std::unique_ptr<T> take_unique_xptr(SEXP &sexp) {
  Rcpp::XPtr<std::unique_ptr<T>> ptr(sexp);
  return std::unique_ptr<T>(ptr->release());
}

// Get a pointer to the object within the XPtr without assuming ownership
template<class T>
T* peek_unique_xptr(SEXP &sexp) {
  Rcpp::XPtr<std::unique_ptr<T>> ptr(sexp);
  return ptr->get();
}