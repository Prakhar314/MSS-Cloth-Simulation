#pragma once
#include "./camera.hpp"
#include "./hw.hpp"
#include "./shape.hpp"
#include <glm/glm.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace GL = COL781::OpenGL;

class Simulator {
  GL::Rasterizer r;
  GL::ShaderProgram program;
  COL781::CameraControl camCtl;
  std::vector<Shape *> shapes;

public:
  ~Simulator() {
  }

  void init(const std::string &filename, const int width, const int height);
  void addShape(Shape* shape);
  void update(float t);
  void render();
  void setCameraZ(float z);
  bool shouldQuit() { return r.shouldQuit(); }
};
