#pragma once
#include "./hw.hpp"
#include <cstdint>
#include <glm/fwd.hpp>
#include <glm/glm.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace GL = COL781::OpenGL;

#define PI 3.14159265358979323846f
#define COLLISION_ERR 0.02f

class Shape {
protected:
  uint32_t nv, nn, nt;
  glm::vec3 *vertices = nullptr;
  glm::vec3 *normals = nullptr;
  glm::ivec3 *triangles = nullptr;
  glm::vec3 *velocities = nullptr;
  glm::vec3 *acc = nullptr;

  glm::vec3 constantVelocity = glm::vec3(0.0f);
  glm::vec3 constantOmega = glm::vec3(0.0f);

  glm::mat4 transform = glm::mat4(1.0f);

  float e = 0.5f, mu = 0.5f;

  void initBuffers();
  void recomputeNormals();
  virtual glm::vec3 getCenter() = 0;

public:
  GL::Object object;
  GL::AttribBuf vertexBuf, normalBuf;
  glm::vec3 color = glm::vec3(1.0, 0.5, 0.0);

  virtual void update(float t, const std::vector<Shape *> &shapes);

  void setConstantVelocity(const glm::vec3 &v);
  void setConstantOmega(const glm::vec3 &omega);
  void setTransform(const glm::mat4 &M);

  void setMaterial(float e, float mu) {
    this->e = e;
    this->mu = mu;
  }

  glm::mat4 getTransform() { return transform; }

  virtual bool isSheet() { return false; }

  virtual bool checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                              glm::vec3 *dp, glm::vec3 *dv) {
    return false;
  }

  void getDV(const glm::vec3 &v, const glm::vec3 &poc, const glm::vec3 &normal,
             float mass, glm::vec3 *dv);

  void getObject(uint32_t &nv, glm::vec3 *&vertices, uint32_t &nn,
                 glm::vec3 *&normals, uint32_t &nt, glm::ivec3 *&triangles) {
    nv = this->nv;
    nn = this->nn;
    nt = this->nt;
    vertices = this->vertices;
    normals = this->normals;
    triangles = this->triangles;
  }

  virtual ~Shape() {
    if (vertices == nullptr)
      return;
    delete[] vertices;
    delete[] normals;
    delete[] triangles;
    if (isSheet()) {
      delete[] velocities;
      delete[] acc;
    }
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
  uint32_t nSprings = 0, width = 0, height = 0;
  std::vector<uint32_t> fixedParticles;
  float spacing = 0;
  float mass = 0.001f;
  float ksStr = 0.12f, kdStr = 0.0012f;
  float ksBend = 0.01f, kdBend = 0.0001f;
  float ksShear = 0.06f, kdShear = 0.0006f;
  float g = -0.98f;

  bool usePBD = true;
  bool selfCollisions = true;
  glm::vec3 *oldPositions = nullptr;
  std::vector<uint16_t> *bins = nullptr;

public:
  Sheet() {}

  bool isSheet() override { return true; }

  void selfCollide(bool updateVel = true);
  void collide(Shape *s, bool updateVel = true);

  void setSelfCollisions(bool selfCollisions);
  void setUsePBD(bool usePBD);
  void setDimensions(uint32_t width, uint32_t height, float spacing);
  void setGravity(float g) { this->g = g; }
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
  void fixParticle(uint32_t i, uint32_t j);
  void init();
  ~Sheet() {
    if (vertices == nullptr)
      return;
    delete[] vertices;
    delete[] normals;
    delete[] triangles;
    delete[] velocities;
    delete[] acc;
    delete[] springs;
    if (bins != nullptr) {
      delete[] bins;
    }
    if (oldPositions != nullptr) {
      delete[] oldPositions;
    }
  }

  glm::vec3 getCenter() override;

  void update(float t, const std::vector<Shape *> &shapes) override;
};

class Sphere : public Shape {
  float radius = 0.0f;
  uint32_t nLat = 0;
  uint32_t nLong = 0;

public:
  void setDimensions(float radius, uint32_t nLat, uint32_t nLong) {
    this->radius = radius;
    this->nLat = nLat;
    this->nLong = nLong;
  }

  bool checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                      glm::vec3 *dp, glm::vec3 *dv) override;

  glm::vec3 getCenter() override;
  void init();
};

class Plane : public Shape {
  glm::vec3 normal = glm::vec3(0, 1, 0);

public:
  void setDimensions(glm::vec3 normal) { this->normal = normal; }

  bool checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                      glm::vec3 *dp, glm::vec3 *dv) override;

  glm::vec3 getCenter() override;
  void init();
};

class Cylinder : public Shape {
  float radius = 0.0f;
  float height = 0.0f;
  uint32_t nLong = 0;
  uint32_t nLat = 0;

public:
  void setDimensions(float radius, float height, uint32_t nLong,
                     uint32_t nLat) {
    this->radius = radius;
    this->height = height;
    this->nLong = nLong;
    this->nLat = nLat;
  }

  bool checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                      glm::vec3 *dp, glm::vec3 *dv) override;

  glm::vec3 getCenter() override;
  void init();
};
