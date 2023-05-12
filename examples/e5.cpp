#include "../src/shape.hpp"
#include "../src/simulator.hpp"
#include <glm/gtx/transform.hpp>
#include <iostream>

// hung cloth, with or without self collision

int main(int argc, char *argv[]) {
  int width = 640, height = 480;

  glm::mat4 m = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, 0.0f, 0.0f));
  m = glm::rotate(m, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));

  Simulator s;
  s.init("Animation", width, height);

  Sheet *sheet = new Sheet();
  float mass = 0.001f;
  float ksStr = 0.12f, kdStr = 0.0012f;
  float ksBend = 0.01f, kdBend = 0.0001f;
  float ksShear = 0.06f, kdShear = 0.0006f;
  float g = -0.98f;
  sheet->setSelfCollisions(argc > 1);
  sheet->setDimensions(20, 20, 0.05f);
  sheet->setMass(mass / 4);
  sheet->setGravity(g / 4);
  sheet->setSpringConstants(ksStr / 8, kdStr / 8, ksBend / 8, kdBend / 8,
                            ksShear / 8, kdShear / 8);
  sheet->init();
  sheet->setTransform(m);
  s.addShape(sheet);
  // Sheet *sheet = new Sheet(argc > 1);
  // sheet->setDimensions(10, 10, 0.1f);
  // sheet->init();
  // sheet->setTransform(m);
  // s.addShape(sheet);

  m = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -1.0f, 0.0f));
  Plane *p = new Plane();
  p->setDimensions(glm::vec3(0, 1, 0));
  p->init();
  p->setTransform(m);
  p->color = glm::vec3(0.5f, 0.5f, 0.5f);
  s.addShape(p);

  float t = SDL_GetTicks64() / 1e3;
  while (!s.shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    s.update((tnew - t) * 1.0f);
    s.render();
    t = tnew;
  }
}
