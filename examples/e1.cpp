#include "../src/shape.hpp"
#include "../src/simulator.hpp"
#include <glm/gtx/transform.hpp>
#include <iostream>

// hung cloth, with or without pbd

int main(int argc, char **argv) {
  int width = 640, height = 480;

  glm::mat4 m = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, 0.0f, 0.0f));

  Simulator s;
  s.init("Animation", width, height);
  Sheet *sheet = new Sheet();
  sheet->setUsePBD(argc > 1);
  sheet->setDimensions(10, 10, 0.1f);
  sheet->init();
  sheet->fixParticle(0, 10);
  sheet->fixParticle(10, 10);
  sheet->setTransform(m);
  s.addShape(sheet);

  float t = SDL_GetTicks64() / 1e3;
  while (!s.shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    s.update((tnew - t) * 1.0f);
    s.render();
    t = tnew;
  }
}
