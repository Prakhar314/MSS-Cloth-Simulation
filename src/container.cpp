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
  if (timeSteps.size() > 0) {
    int lastTime = ceil(timeSteps[timeSteps.size() - 1].first + 1.0f);
    int currentTime = ceil(t);
    int cycle = currentTime / lastTime;
    t -= cycle * lastTime;
    for (size_t i = 0; i < timeSteps.size(); i++) {
      if (t < timeSteps[i].first) {
        float t0 = 0.0f;
        float t1 = 0.0f;
        float t2 = timeSteps[i].first;
        float t3 = 0.0f;
        float a0 = 0.0f;
        float a1 = 0.0f;
        float a2 = timeSteps[i].second;
        float a3 = 0.0f;

        if (i == 0) {
          angle = timeSteps[i].second * t / t2;
          break;
        }
        if (i > 0) {
          t1 = timeSteps[i - 1].first;
          a1 = timeSteps[i - 1].second;
        }
        if (i > 1) {
          t0 = timeSteps[i - 2].first;
          a0 = timeSteps[i - 2].second;
        }
        if (i < timeSteps.size() - 1) {
          t3 = timeSteps[i + 1].first;
          a3 = timeSteps[i + 1].second;
        } else {
          t3 = t2 + (t2 - t1);
          a3 = a2;
        }

        float dq_dt_1 = (a2 - a1) / (t2 - t1) * (t1 - t0) / (t2 - t0) +
                        (a1 - a0) / (t1 - t0) * (t2 - t1) / (t2 - t0);
        float dq_dt_2 = (a3 - a2) / (t3 - t2) * (t2 - t1) / (t3 - t1) +
                        (a2 - a1) / (t2 - t1) * (t3 - t2) / (t3 - t1);
        float dt = (t - t1) / (t2 - t1);
        float dt_2 = dt * dt;
        float dt_3 = dt_2 * dt;
        angle = (2 * dt_3 - 3 * dt_2 + 1) * a1 +
                (dt_3 - 2 * dt_2 + dt) * dq_dt_1 + (-2 * dt_3 + 3 * dt_2) * a2 +
                (dt_3 - dt_2) * dq_dt_2;

        break;
      }
    }
  }
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

void ShapeContainer::setRotation(glm::vec3 rotation_axis) {
  this->rotation_axis = rotation_axis;
}

void ShapeContainer::setScale(glm::vec3 scale) { this->scale = scale; }

void ShapeContainer::setTimeSteps(
    std::vector<std::pair<float, float>> timeSteps) {
  this->timeSteps = timeSteps;
}
