#include "./shape.hpp"

class Character {
  Simulator *simulator = nullptr;

  vector<Cylinder *> bones;
public:
  
  void useSimulator(Simulator *s) {
    this->simulator = s;
    for (auto bone : bones) {
      s->addShape(bone);
    }
  }

  void addBone(Cylinder *bone) {
    bones.push_back(bone);
    if (simulator != nullptr) {
      simulator->addShape(bone);
    }
  }

  void update(float t) {
    
  }
}
