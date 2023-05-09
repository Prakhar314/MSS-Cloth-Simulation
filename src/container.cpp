#include "./container.hpp"
#include "glm/gtx/transform.hpp"

ShapeContainer::ShapeContainer() {
  recenter = glm::vec3(0.0f);
  position = glm::vec3(0.0f);
  scale = glm::vec3(1.0f);
  rotation_axis = glm::vec3(0.0f, 1.0f, 0.0f);
  angle = 0.0f;
}

void ShapeContainer::addShape(Shape *s, Simulator *sim) {
  shapes.push_back(s);
  sim->addShape(s);
}

void ShapeContainer::addChild(ShapeContainer *c) { children.push_back(c); }

void ShapeContainer::update(float t, glm::mat4 p) {
  glm::mat4 m = glm::translate(p, position);
  m = glm::rotate(m, glm::radians(angle), rotation_axis);
  m = glm::translate(m, recenter);
  m = glm::scale(m, scale);
  for (Shape *s : shapes) {
    s->setTransform(m);
  }
  for (ShapeContainer *c : children) {
    c->update(t, m);
  }
}

void ShapeContainer::setRecenter(glm::vec3 recenter) {
  this->recenter = recenter;
}

void ShapeContainer::setPosition(glm::vec3 position) {
  this->position = position;
}

void ShapeContainer::setRotation(float angle, glm::vec3 rotation_axis) {
  this->angle = angle;
  this->rotation_axis = rotation_axis;
}

void ShapeContainer::setScale(glm::vec3 scale) { this->scale = scale; }
