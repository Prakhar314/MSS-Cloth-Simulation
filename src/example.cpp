#include "./simulator.hpp"

int main() {
  int width = 640, height = 480;

  Simulator s;
  s.init("Animation", width, height);
  s.createScene();

  while (!s.shouldQuit()) {
    float t = SDL_GetTicks64() / 1e3;
    s.update(t);
    s.render();
  }
}
