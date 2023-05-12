
#include "../src/shape.hpp"
#include "../src/simulator.hpp"
#include <glm/gtx/transform.hpp>
#include <iostream>

// collision with sphere

int main() {
  int width = 640, height = 480;

  Simulator s;
  s.init("Animation", width, height);

  glm::mat4 m = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, -0.0f, 0.0f));

  Sheet *sheet = new Sheet();
  sheet->setDimensions(10, 10, 0.1f);
  sheet->init();
  sheet->fixParticle(0, 10);
  sheet->fixParticle(10, 10);
  sheet->setTransform(m);
  s.addShape(sheet);

  m = glm::translate(m, glm::vec3(0.0f, -0.9f, -0.1f));

  Sphere *sphere = new Sphere();
  sphere->setDimensions(0.2f, 20, 20);
  sphere->init();
  sphere->setTransform(m);
  sphere->setMaterial(0.5f, 0.5f);
  sphere->color = glm::vec3(1.0f, 0.0f, 0.0f);
  s.addShape(sphere);

  float t = SDL_GetTicks64() / 1e3;
  while (!s.shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    s.update((tnew - t) * 1.0f);
    s.render();
    t = tnew;
  }
}
