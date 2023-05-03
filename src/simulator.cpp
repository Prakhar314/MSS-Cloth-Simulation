#include "simulator.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;
using namespace std;

void Simulator::init(const std::string &filename, const int width,
                     const int height) {
  r.initialize(filename, width, height);
  camCtl.initialize(width, height);
  camCtl.camera.setCameraView(vec3(0.5, -0.5, 1.5), vec3(0.5, -0.5, 0.0),
                              vec3(0.0, 1.0, 0.0));
  program = r.createShaderProgram(r.vsPhongShading(), r.fsPhongShading());
}

void Simulator::addSheet(uint32_t width, uint32_t height, float spacing) {

  object = r.createObject();

  sheet = new Sheet(width, height, spacing);
  uint32_t nv, nn, nt;
  sheet->getObject(nv, vertices, nn, normals, nt, triangles);

  vertexBuf = r.createVertexAttribs(object, 0, nv, vertices);
  normalBuf = r.createVertexAttribs(object, 1, nn, normals);
  r.createTriangleIndices(object, nt, triangles);
}

void Simulator::update(float t) {

  sheet->update(t);
  uint32_t nv, nn, nt;
  sheet->getObject(nv, vertices, nn, normals, nt, triangles);

  r.updateVertexAttribs(vertexBuf, nv, vertices);
  r.updateVertexAttribs(normalBuf, nv, normals);
}

void Simulator::render() {
  camCtl.update();
  Camera &camera = camCtl.camera;

  r.clear(vec4(1.0, 1.0, 1.0, 1.0));
  r.enableDepthTest();
  r.useShaderProgram(program);

  r.setUniform(program, "model", glm::mat4(1.0));
  r.setUniform(program, "view", camera.getViewMatrix());
  r.setUniform(program, "projection", camera.getProjectionMatrix());
  r.setUniform(program, "lightPos", camera.position);
  r.setUniform(program, "viewPos", camera.position);
  r.setUniform(program, "lightColor", vec3(1.0f, 1.0f, 1.0f));

  r.setupFilledFaces();
  r.setUniform(program, "objectColor", vec3(1.0f, 0.5f, 0.0f));
  r.drawObject(object);

  r.setupWireFrame();
  r.setUniform(program, "objectColor", vec3(0.0f, 0.0f, 0.0f));
  r.drawObject(object);

  r.show();
}

Sheet::Sheet(uint32_t width, uint32_t height, float spacing) {
  this->width = width;
  this->height = height;
  this->spacing = spacing;

  uint32_t m = width;
  uint32_t n = height;
  nv = (m + 1) * (n + 1);
  nn = nv;
  nt = 2 * m * n;
  vertices = new glm::vec3[nv];
  velocities = new glm::vec3[nv];
  normals = new glm::vec3[nn];
  triangles = new glm::ivec3[nt];
  acc = new glm::vec3[nv];
  // Create the vertices, normals and velocities
  for (size_t i = 0; i < m + 1; i++) {
    for (size_t j = 0; j < n + 1; j++) {
      vertices[i * (n + 1) + j] =
          glm::vec3(1.0f * i * spacing, 0, -1.0f * j * spacing) -
          glm::vec3(0.0, 0, 0);
      normals[i * (n + 1) + j] = glm::vec3(0, 0, 1);
      velocities[i * (n + 1) + j] = glm::vec3(0, 0, 0);
    }
  }
  // Create the triangles
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      triangles[2 * (i * n + j)] = glm::ivec3(
          i * (n + 1) + j, i * (n + 1) + j + 1, (i + 1) * (n + 1) + j);
      triangles[2 * (i * n + j) + 1] =
          glm::ivec3(i * (n + 1) + j + 1, (i + 1) * (n + 1) + j + 1,
                     (i + 1) * (n + 1) + j);
    }
  }

  // Create the springs
  uint32_t nStructuralSprings = 2 * m * n + m + n;
  uint32_t nShearSprings = 2 * m * n;
  uint32_t nBendSprings = (m + 1) * (n - 1) + (m - 1) * (n + 1);
  nSprings = nStructuralSprings + nShearSprings + nBendSprings;
  springs = new Spring[nStructuralSprings + nShearSprings + nBendSprings];
  uint32_t springIndex = 0;
  // Structural springs
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n + 1; j++) {
      springs[springIndex++] =
          Spring(i * (n + 1) + j, (i + 1) * (n + 1) + j, spacing, ksStr, kdStr);
    }
  }
  for (size_t i = 0; i < m + 1; i++) {
    for (size_t j = 0; j < n; j++) {
      springs[springIndex++] =
          Spring(i * (n + 1) + j, i * (n + 1) + j + 1, spacing, ksStr, kdStr);
    }
  }
  // Shear springs
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      springs[springIndex++] =
          Spring(i * (n + 1) + j, (i + 1) * (n + 1) + j + 1, spacing * sqrt(2),
                 ksShear, kdShear);
      springs[springIndex++] =
          Spring(i * (n + 1) + j + 1, (i + 1) * (n + 1) + j, spacing * sqrt(2),
                 ksShear, kdShear);
    }
  }
  // Bend springs
  for (size_t i = 0; i < m + 1; i++) {
    for (size_t j = 0; j < n - 1; j++) {
      springs[springIndex++] = Spring(i * (n + 1) + j, i * (n + 1) + j + 2,
                                      2 * spacing, ksBend, kdBend);
    }
  }
  for (size_t i = 0; i < m - 1; i++) {
    for (size_t j = 0; j < n + 1; j++) {
      springs[springIndex++] = Spring(i * (n + 1) + j, (i + 2) * (n + 1) + j,
                                      2 * spacing, ksBend, kdBend);
    }
  }
  // print first spring
  assert(springIndex == nSprings);
}

void Shape::recomputeNormals() {
  for (size_t i = 0; i < nv; i++) {
    normals[i] = glm::vec3(0, 0, 0);
  }
  for (size_t i = 0; i < nt; i++) {
    glm::vec3 e1 = vertices[triangles[i].y] - vertices[triangles[i].x];
    glm::vec3 e2 = vertices[triangles[i].z] - vertices[triangles[i].x];
    glm::vec3 normal =
        glm::normalize(glm::cross(e2, e1)) / glm::length(e1) / glm::length(e2);
    normals[triangles[i].x] += normal;
    normals[triangles[i].y] += normal;
    normals[triangles[i].z] += normal;
  }
  for (size_t i = 0; i < nv; i++) {
    normals[i] = glm::normalize(normals[i]);
  }
}

void Sheet::update(float t) {
  // Compute the forces
  glm::vec3 gravity = glm::vec3(0, g, 0);
  fill_n(acc, nv, gravity);
  for (size_t i = 0; i < nSprings; i++) {
    Spring s = springs[i];
    uint32_t i1 = s.i;
    uint32_t i2 = s.j;
    glm::vec3 p1 = vertices[i1];
    glm::vec3 p2 = vertices[i2];
    glm::vec3 v1 = velocities[i1];
    glm::vec3 v2 = velocities[i2];
    glm::vec3 dp = p1 - p2;
    glm::vec3 dv = v1 - v2;
    float e = glm::length(dp) / s.restLength - 1;
    float e_dot = glm::dot(dv, dp) / glm::length(dp) / s.restLength;
    glm::vec3 f =
        (s.ks * e + s.kd * e_dot) * dp / glm::length(dp);
    acc[i1] -= f / mass;
    acc[i2] += f / mass;
  }
  // keep the top left and right corners fixed
  acc[height] = glm::vec3(0, 0, 0);
  acc[width * height + width + height] = glm::vec3(0, 0, 0);
  for (size_t i = 0; i < nv; i++) {
    velocities[i] += acc[i] * t;
    vertices[i] += velocities[i] * t;
    assert(glm::length(velocities[i]) < 1000);
  }
  recomputeNormals();
}
