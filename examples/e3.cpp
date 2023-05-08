#include "../src/shape.hpp"
#include "../src/simulator.hpp"
#include <glm/gtx/transform.hpp>
#include <iostream>

int main() {
  int width = 640, height = 480;

  Simulator s;
  s.init("Animation", width, height);

  glm::mat4 m = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, 0.0f, 0.0f));

  Sheet *sheet = new Sheet();
  sheet->setDimensions(10, 10, 0.1f);
  sheet->init();
  sheet->fixParticle(0, 10);
  sheet->fixParticle(10, 10);
  sheet->initTransform(m);
  s.addShape(sheet);

  m = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -0.9f, 0.0f));
  Plane *p = new Plane();
  p->setDimensions(glm::vec3(0, 1, 0));
  p->init();
  p->initTransform(m);
  p->color = glm::vec3(0.5f, 0.5f, 0.5f);
  s.addShape(p);

  m = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, -2.8f, -0.25f));
  Cylinder *c = new Cylinder();
  c->setDimensions(0.1f, 2.0f, 20, 10);
  c->init();
  c->initTransform(m);
  c->initConstantVelocity(glm::vec3(0.0f, 0.2f, 0.0f));
  c->color = glm::vec3(0.0f, 0.0f, 0.5f);
  s.addShape(c);

  m = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, -0.7f, -0.0f));
  Sphere *sphere = new Sphere();
  sphere->setDimensions(0.2f, 20, 20);
  sphere->color = glm::vec3(1.0f, 0.0f, 0.0f);
  sphere->init();
  sphere->initTransform(m);
  s.addShape(sphere);

  float t = SDL_GetTicks64() / 1e3;
  while (!s.shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    s.update((tnew - t) * 1.0f);
    s.render();
    t = tnew;
  }
}
