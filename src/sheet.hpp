#pragma once
#include "./shape.hpp"

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

