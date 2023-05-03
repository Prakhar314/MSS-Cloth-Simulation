#pragma once
#include "./camera.hpp"
#include "./hw.hpp"
#include <glm/glm.hpp>
#include <string>
#include <tuple>

namespace GL = COL781::OpenGL;

class Shape {
protected:
  uint32_t nv, nn, nt;
  glm::vec3 *vertices = nullptr;
  glm::vec3 *normals = nullptr;
  glm::ivec3 *triangles = nullptr;

  void recomputeNormals();

public:
  virtual void update(float t) {}
  void getObject(uint32_t &nv, glm::vec3 *&vertices, uint32_t &nn,
                 glm::vec3 *&normals, uint32_t &nt, glm::ivec3 *&triangles) {
    nv = this->nv;
    nn = this->nn;
    nt = this->nt;
    vertices = this->vertices;
    normals = this->normals;
    triangles = this->triangles;
  }
  ~Shape() {
    if (vertices == nullptr)
      return;
    delete[] vertices;
    delete[] normals;
    delete[] triangles;
  }
};

struct Spring {
  uint32_t i, j;
  float restLength, ks, kd;

public:
  Spring(uint32_t i, uint32_t j, float restLength, float ks, float kd)
      : i(i), j(j), restLength(restLength), ks(ks), kd(kd) {}
  Spring() {}
};

class Sheet : public Shape {
  Spring *springs = nullptr;
  glm::vec3 *velocities = nullptr;
  glm::vec3 *acc = nullptr;
  uint32_t nSprings, width, height, spacing;
  float mass = 0.001f;
  float ksStr = 0.12f, kdStr = 0.0012f;
  float ksBend = 0.01f, kdBend = 0.0001f;
  float ksShear = 0.06f, kdShear = 0.0006f;
  float g = -0.98f;

public:
  Sheet(uint32_t width, uint32_t height, float spacing);
  ~Sheet() {
    if (springs == nullptr)
      return;
    delete[] springs;
  }

  void update(float t) override;
};

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
  Sheet *sheet;

public:
  ~Simulator() {
    delete[] vertices;
    delete[] normals;
    delete[] triangles;
    delete sheet;
  }

  void init(const std::string &filename, const int width, const int height);
  void addSheet(uint32_t width, uint32_t height, float spacing);
  void update(float t);
  void render();

  bool shouldQuit() { return r.shouldQuit(); }
};
