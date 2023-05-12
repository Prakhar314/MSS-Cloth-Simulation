#include "./shape.hpp"
#include <cassert>
#include <cstdint>
#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

using namespace std;

void Shape::initBuffers() {
  vertices = new glm::vec3[nv];
  normals = new glm::vec3[nn];
  triangles = new glm::ivec3[nt];

  if (isSheet()) {
    velocities = new glm::vec3[nv];
    acc = new glm::vec3[nv];

    std::fill_n(velocities, nv, glm::vec3(0, 0, 0));
    std::fill_n(acc, nv, glm::vec3(0, 0, 0));
  }
}

void Shape::setTransform(const glm::mat4 &M) { transform = M; }

void Shape::update(float t, const vector<Shape *> &shapes) {
  if (!isSheet()) {
    glm::vec3 c = getCenter();
    for (int i = 0; i < nv; ++i) {
      glm::vec3 relV = glm::cross(constantOmega, vertices[i] - c);
      vertices[i] += (constantVelocity + relV) * t;
    }
  }
}

void Shape::setConstantVelocity(const glm::vec3 &v) { constantVelocity = v; }
void Shape::setConstantOmega(const glm::vec3 &omega) { constantOmega = omega; }

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

void Shape::getDV(const glm::vec3 &v, const glm::vec3 &poc,
                  const glm::vec3 &normal, float mass, glm::vec3 *dv) {
  if (glm::dot(v, normal) > 0) {
    *dv = glm::vec3(0, 0, 0);
    return;
  }

  glm::vec3 velocityOfContact =
      constantVelocity + glm::cross(constantOmega, poc - getCenter());
  glm::vec3 v_rel = v - velocityOfContact;

  glm::vec3 v_n = glm::dot(v_rel, normal) * normal;
  glm::vec3 v_t = v_rel - v_n;
  glm::vec3 j_n = -(1 + e) * mass * v_n;
  float v_t_mag = glm::length(v_t);
  float j_t_mag = min(mu * glm::length(j_n), mass * v_t_mag);
  glm::vec3 j_t = -j_t_mag * v_t;
  // v_t may be 0
  if (v_t_mag < 0.001f) {
    j_t = glm::vec3(0, 0, 0);
  } else {
    j_t /= v_t_mag;
  }

  *dv = (j_n + j_t) / mass;
  assert(glm::length(*dv) < 1000.0f);
}

void Sphere::init() {
  uint32_t m = nLong;
  uint32_t n = nLat;
  assert(radius > 0 && nLat > 0 && nLong > 0);
  nv = m * (n - 1) + 2;
  nt = 2 * m * (n - 1);
  nn = nv;

  initBuffers();

  float pi = 3.1415926535;
  float long_angle = 2 * pi / m;
  float lat_angle = pi / n;
  float lat_off = -pi / 2 + lat_angle;
  // Create the vertices over a sphere
  for (size_t i = 0; i < n - 1; i++) {
    for (size_t j = 0; j < m; j++) {
      float x = radius * cos(long_angle * j) * cos(lat_angle * i + lat_off);
      float y = radius * sin(long_angle * j) * cos(lat_angle * i + lat_off);
      float z = radius * sin(lat_angle * i + lat_off);
      vertices[i * m + j] = glm::vec3(x, y, z);
      normals[i * m + j] = glm::vec3(x, y, z) / radius;
    }
  }
  vertices[m * (n - 1)] = glm::vec3(0, 0, radius);
  normals[m * (n - 1)] = glm::vec3(0, 0, 1);
  vertices[m * (n - 1) + 1] = glm::vec3(0, 0, -radius);
  normals[m * (n - 1) + 1] = glm::vec3(0, 0, -1);
  // Create the triangles
  for (size_t i = 0; i < n - 2; i++) {
    for (size_t j = 0; j < m; j++) {
      triangles[2 * (i * m + j)] =
          glm::ivec3(i * m + j, i * m + (j + 1) % m, (i + 1) * m + j);
      triangles[2 * (i * m + j) + 1] = glm::ivec3(
          i * m + (j + 1) % m, (i + 1) * m + (j + 1) % m, (i + 1) * m + j);
    }
  }
  for (size_t j = 0; j < m; j++) {
    triangles[2 * (m * (n - 2) + j)] =
        glm::ivec3(m * (n - 2) + j, m * (n - 2) + (j + 1) % m, m * (n - 1));
    triangles[2 * (m * (n - 2) + j) + 1] =
        glm::ivec3(j, m * (n - 1) + 1, (j + 1) % m);
  }
}

glm::vec3 Sphere::getCenter() {
  return (vertices[nv - 1] + vertices[nv - 2]) / 2.0f;
}

bool Sphere::checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                            glm::vec3 *dp, glm::vec3 *dv) {
  glm::vec3 origin = getCenter();
  float dist = glm::length(p - origin);
  float c_radius = radius + COLLISION_ERR;
  if (dist > c_radius) {
    return false;
  }

  glm::vec3 poc = glm::normalize(p - origin) * c_radius + origin;

  glm::vec3 normal = glm::normalize(p - origin);
  getDV(v, poc, normal, m, dv);
  *dp = poc - p;
  return true;
}

void Plane::init() {
  assert(glm::length(normal) == 1);
  nv = 4;
  nn = 4;
  nt = 2;

  initBuffers();

  triangles[0] = glm::ivec3(0, 1, 2);
  triangles[1] = glm::ivec3(1, 3, 2);

  glm::vec3 init_normal = glm::vec3(0, 1, 0);
  float angle = acos(glm::dot(init_normal, normal));

  glm::mat4 rot = glm::mat4(1.0f);
  if (abs(angle) > 0.01f) {
    glm::vec3 perp = glm::normalize(glm::cross(init_normal, normal));
    glm::mat4 rot = glm::rotate(glm::mat4(1.0f), angle, perp);
  }

  vertices[0] = rot * glm::vec4(-10, 0, -10, 1);
  vertices[1] = rot * glm::vec4(10, 0, -10, 1);
  vertices[2] = rot * glm::vec4(-10, 0, 10, 1);
  vertices[3] = rot * glm::vec4(10, 0, 10, 1);

  for (size_t i = 0; i < 4; i++) {
    normals[i] = rot * glm::vec4(0, 1, 0, 0);
  }
}

glm::vec3 Plane::getCenter() { return (vertices[0] + vertices[3]) / 2.0f; }

bool Plane::checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                           glm::vec3 *dp, glm::vec3 *dv) {
  glm::vec3 normal = glm::normalize(
      glm::cross(vertices[2] - vertices[0], vertices[1] - vertices[0]));
  float dist = glm::dot(p - vertices[0], normal);
  if (dist > COLLISION_ERR) {
    return false;
  }

  glm::vec3 poc = p - normal * dist;
  getDV(v, poc, normal, m, dv);
  *dp = (COLLISION_ERR - dist) * normal;
  return true;
}

void Cylinder::init() {
  nv = (nLat + 1) * nLong + 2;
  nn = nv;
  nt = nLat * nLong * 2 + nLong * 2;

  initBuffers();

  float lat_angle = 2 * M_PI / nLong;
  float spacing = height / nLat;

  for (size_t i = 0; i < nLat + 1; i++) {
    for (size_t j = 0; j < nLong; j++) {
      float x = spacing * i;
      float y = radius * cos(lat_angle * j);
      float z = radius * sin(lat_angle * j);
      vertices[i * nLong + j] =
          glm::vec3(x, y, z) - glm::vec3(height / 2, 0, 0);
      normals[i * nLong + j] = glm::vec3(0, y, z) / radius;
    }
  }

  vertices[nv - 2] = glm::vec3(-height / 2, 0, 0);
  vertices[nv - 1] = glm::vec3(height / 2, 0, 0);
  normals[nv - 2] = glm::vec3(-1, 0, 0);
  normals[nv - 1] = glm::vec3(1, 0, 0);

  for (size_t i = 0; i < nLat; i++) {
    for (size_t j = 0; j < nLong; j++) {
      triangles[2 * (i * nLong + j)] = glm::ivec3(
          i * nLong + j, i * nLong + (j + 1) % nLong, (i + 1) * nLong + j);
      triangles[2 * (i * nLong + j) + 1] =
          glm::ivec3(i * nLong + (j + 1) % nLong,
                     (i + 1) * nLong + (j + 1) % nLong, (i + 1) * nLong + j);
    }
  }

  for (size_t j = 0; j < nLong; j++) {
    triangles[2 * (nLat * nLong + j)] =
        glm::ivec3(nLat * nLong + j, nLat * nLong + (j + 1) % nLong, nv - 1);
    triangles[2 * (nLat * nLong + j) + 1] =
        glm::ivec3(j, nv - 2, (j + 1) % nLong);
  }
}

glm::vec3 Cylinder::getCenter() {
  return (vertices[nv - 1] + vertices[nv - 2]) / 2.0f;
}

bool Cylinder::checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                              glm::vec3 *dp, glm::vec3 *dv) {
  glm::vec3 bottom = vertices[nv - 2];
  glm::vec3 top = vertices[nv - 1];
  glm::vec3 normal = glm::normalize(top - bottom);

  float distAlongAxis = glm::dot(p - bottom, normal);
  float distFromAxis = glm::length(p - (bottom + distAlongAxis * normal));

  if (distAlongAxis < -COLLISION_ERR ||
      distAlongAxis > height + COLLISION_ERR ||
      distFromAxis > radius + COLLISION_ERR) {
    return false;
  }

  // check if close to faces or curve
  if (distAlongAxis < COLLISION_ERR) {
    glm::vec3 poc = p - normal * (distAlongAxis + COLLISION_ERR);
    getDV(v, poc, -normal, m, dv);
    *dp = poc - p;
  } else if (distAlongAxis > height - COLLISION_ERR) {
    glm::vec3 poc =
        p - normal * distAlongAxis + top - bottom + COLLISION_ERR * normal;
    getDV(v, poc, normal, m, dv);
    *dp = poc - p;
  } else {
    glm::vec3 poc = bottom + distAlongAxis * normal +
                    (radius + COLLISION_ERR) *
                        glm::normalize(p - bottom - distAlongAxis * normal);
    normal = glm::normalize(p - bottom - distAlongAxis * normal);
    getDV(v, poc, normal, m, dv);
    *dp = poc - p;
  }
  return true;
}
