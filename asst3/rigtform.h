#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS175_EPS2);
  }
  // creating the RigTForm when given a Cvec3 and Quaternion
  RigTForm(const Cvec3& t, const Quat& r) {
	//TODO
     t_ = t;
	 r_ = r;
  }
  // what is it asking us to do??? so confused!
  explicit RigTForm(const Cvec3& t) {
	 // TODO
	  t_ = t;
	  r_ = Quat();
  }
  //what are we supposed to do?? so confused!!
  explicit RigTForm(const Quat& r) {
    // TODO
	  t_ = Cvec3(0, 0, 0);
	  r_ = r;
  }

  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }
  // multiplying a Cvec4
  Cvec4 operator * (const Cvec4& a) const {
	  // A.r * c  + A.t if A is the RigTForm and c is the cvec4 (according to the book)
	  return (r_* a + Cvec4(t_, 0) * a[3]);
  }
  // does it matter which one we use as one and two?
  RigTForm operator * (const RigTForm& a) const {
	  //translation  t1 + r1t2, rotation r1r2
	  return RigTForm(t_ + Cvec3((r_* Cvec4(a.getTranslation(), 0))), r_ * a.getRotation());
  }
};
// syntax here right?
inline RigTForm inv(const RigTForm& tform) {
	return RigTForm(Cvec3(inv(tform.getRotation()) * Cvec4(tform.getTranslation(), 0) * -1), inv(tform.getRotation()));
	// TODO
	// answer.t = -r-1(t)
	// anaswer.r = inv(Rbt.r)
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}
// convert rigTform to a matrix, from the book (is getTranslation the same as makeTranslation)??
inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
	// TODO
	Matrix4 T = Matrix4::makeTranslation(tform.getTranslation());
	Matrix4 R = quatToMatrix(tform.getRotation());
	Matrix4 m = T * R;
	return m;
}

#endif