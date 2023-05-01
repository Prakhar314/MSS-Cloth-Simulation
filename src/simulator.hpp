#pragma once
#include "./camera.hpp"
#include "./hw.hpp"
#include <glm/glm.hpp>
#include <string>

namespace GL = COL781::OpenGL;

class Simulator {
  GL::Rasterizer r;
  GL::ShaderProgram program;
  GL::Object object;
  GL::AttribBuf vertexBuf, normalBuf;
  COL781::CameraControl camCtl;

  const int nv = 4;
  const int nt = 2;
  glm::vec3 *vertices;
  glm::vec3 *normals;
  glm::ivec3 *triangles;

public:
  ~Simulator() {
    delete[] vertices;
    delete[] normals;
    delete[] triangles;
  }

  void init(const std::string &filename, const int width, const int height);
  void createScene();
  void update(float t);
  void render();

  bool shouldQuit() { return r.shouldQuit(); }
};
