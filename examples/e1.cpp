#include "../src/simulator.hpp"
#include "../src/shape.hpp"
#include <iostream>

int main() {
  int width = 640, height = 480;

  Simulator s;
  s.init("Animation", width, height);
  Sheet *sheet = new Sheet();
  sheet->setDimensions(10, 10, 0.1f);
  sheet->init();
  s.addShape(sheet);

  float t = SDL_GetTicks64() / 1e3;
  while (!s.shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    s.update((tnew - t) * 1.0f);
    s.render();
    t = tnew;
  }
}
