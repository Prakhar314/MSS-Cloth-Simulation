#include "./shape.hpp"
#include "./simulator.hpp"
#include <glm/fwd.hpp>
#include <vector>
#include <utility>

class ShapeContainer {
private:
  std::vector<Shape *> shapes;
  std::vector<ShapeContainer *> children;
  std::vector<std::pair<float, float>> timeSteps;

  glm::vec3 recenter, position, rotation_axis, scale;
  float angle = 0.0f;

public:
  ShapeContainer();

  void setRecenter(glm::vec3 recenter);
  void setPosition(glm::vec3 position);
  void setRotation(glm::vec3 rotation_axis);
  void setScale(glm::vec3 scale);

  void setTimeSteps(std::vector<std::pair<float, float>> timeSteps);

  void addShape(Shape *s, Simulator *sim);
  void addChild(ShapeContainer *c);
  void update(float t, glm::mat4 p = glm::mat4(1.0f));

  ~ShapeContainer() {
    for (Shape *s : shapes) {
      delete s;
    }
    for (ShapeContainer *c : children) {
      delete c;
    }
  }
};
