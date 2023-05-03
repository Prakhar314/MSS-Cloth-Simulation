#include "./simulator.hpp"
#include <iostream>

int main() {
  int width = 640, height = 480;

  Simulator s;
  s.init("Animation", width, height);
  s.addSheet(10, 10, 0.1);

  float t = SDL_GetTicks64() / 1e3;
  while (!s.shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    std::cout << tnew - t << std::endl;
    s.update((tnew - t) * 1.0f);
    s.render();
    t = tnew;
  }
}
