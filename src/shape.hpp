#pragma once
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
  GL::Object object;
  GL::AttribBuf vertexBuf, normalBuf;
  glm::vec3 color = glm::vec3(1.0, 0.5, 0.0);

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
  uint32_t nSprings = 0, width = 0, height = 0;
  float spacing = 0;
  float mass = 0.001f;
  float ksStr = 0.12f, kdStr = 0.0012f;
  float ksBend = 0.01f, kdBend = 0.0001f;
  float ksShear = 0.06f, kdShear = 0.0006f;
  float g = -0.98f;

public:
  Sheet() {}
  void setDimensions(uint32_t width, uint32_t height, float spacing);
  void setMass(float mass) { this->mass = mass; }
  void setSpringConstants(float ksStr, float kdStr, float ksBend, float kdBend,
                          float ksShear, float kdShear) {
    this->ksStr = ksStr;
    this->kdStr = kdStr;
    this->ksBend = ksBend;
    this->kdBend = kdBend;
    this->ksShear = ksShear;
    this->kdShear = kdShear;
  }
  void init();
  ~Sheet() {
    if (springs == nullptr)
      return;
    delete[] springs;
  }

  void update(float t) override;
};
