#ifndef DRAWER_H
#define DRAWER_H

#include <vector>
#include <cstddef>
#include <string>
#include <memory>
#include <list>
#include <stdexcept>

#include "scenegraph.h"
#include "asstcommon.h"

class Drawer : public SgNodeVisitor {
protected:
  std::vector<RigTForm> rbtStack_;
  const ShaderState& curSS_;
public:
  Drawer(const RigTForm& initialRbt, const ShaderState& curSS)
    : rbtStack_(1, initialRbt)
    , curSS_(curSS) {}

  virtual bool visit(SgTransformNode& node) {
    rbtStack_.push_back(rbtStack_.back() * node.getRbt());
    return true;
  }

  virtual bool postVisit(SgTransformNode& node) {
    rbtStack_.pop_back();
    return true;
  }

  virtual bool visit(SgShapeNode& shapeNode) {
    const Matrix4 MVM = rigTFormToMatrix(rbtStack_.back()) * shapeNode.getAffineMatrix();
	/*for (int i = 0; i < 4; i++) {
	  for (int j = 0; j < 4; j++) {
		std::cout << (rigTFormToMatrix(rbtStack_.back()))(i,j) << " ";
	  }
	  std::cout << std::endl;
	}*/
    sendModelViewNormalMatrix(curSS_, MVM, normalMatrix(MVM));
    shapeNode.draw(curSS_);
    return true;
  }

  virtual bool postVisit(SgShapeNode& shapeNode) {
    return true;
  }

  const ShaderState& getCurSS() const {
    return curSS_;
  }
};

#endif



