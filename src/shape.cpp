#include "./shape.hpp"
#include <cstdint>
#include <iostream>

using namespace std;

void Shape::initBuffers() {
  vertices = new glm::vec3[nv];
  normals = new glm::vec3[nn];
  triangles = new glm::ivec3[nt];
  velocities = new glm::vec3[nv];
  acc = new glm::vec3[nv];

  std::fill_n(velocities, nv, glm::vec3(0, 0, 0));
  std::fill_n(acc, nv, glm::vec3(0, 0, 0));
}

void Shape::initTransform(const glm::mat4 &M) {
  for (size_t i = 0; i < nv; i++) {
    vertices[i] = glm::vec3(M * glm::vec4(vertices[i], 1.0));
    normals[i] = glm::vec3(M * glm::vec4(normals[i], 0.0));
    velocities[i] = glm::vec3(M * glm::vec4(velocities[i], 0.0));
  }
  origin = glm::vec3(M * glm::vec4(origin, 1.0));
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

void Sheet::setDimensions(uint32_t width, uint32_t height, float spacing) {
  this->width = width;
  this->height = height;
  this->spacing = spacing;
}

void Sheet::fixParticle(uint32_t i, uint32_t j) {
  fixedParticles.push_back(i * (height + 1) + j);
}

void Sheet::init() {
  assert(width > 0 && height > 0 && spacing > 0.0f);
  uint32_t m = width;
  uint32_t n = height;
  nv = (m + 1) * (n + 1);
  nn = nv;
  nt = 2 * m * n;

  initBuffers();

  // Create the vertices, normals and velocities
  for (size_t i = 0; i < m + 1; i++) {
    for (size_t j = 0; j < n + 1; j++) {
      vertices[i * (n + 1) + j] =
          glm::vec3(1.0f * i * spacing, 0, -1.0f * j * spacing) +
          glm::vec3(-spacing / 2 * (m + 1), 0, spacing / 2 * (n + 1));
      normals[i * (n + 1) + j] = glm::vec3(0, 0, 1);
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
    glm::vec3 f = (s.ks * e + s.kd * e_dot) * dp / glm::length(dp);
    acc[i1] -= f / mass;
    acc[i2] += f / mass;
  }
  // keep the top left and right corners fixed
  for (auto &i : fixedParticles) {
    acc[i] = glm::vec3(0, 0, 0);
  }
  // acc[height] = glm::vec3(0, 0, 0);
  // acc[width * height + width + height] = glm::vec3(0, 0, 0);
  for (size_t i = 0; i < nv; i++) {
    velocities[i] += acc[i] * t;
    vertices[i] += velocities[i] * t;
    assert(glm::length(velocities[i]) < 1000);
  }
  recomputeNormals();
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

bool Sphere::checkCollision(const glm::vec3 p, const glm::vec3 v, float m,
                            glm::vec3 *dp, glm::vec3 *dv) {
  float dist = glm::length(p - origin);
  float c_radius = radius * 1.1f;
  if (dist > c_radius) {
    return false;
  }

  glm::vec3 v_n =
      glm::dot(v, glm::normalize(p - origin)) * glm::normalize(p - origin);
  glm::vec3 v_t = v - v_n;
  glm::vec3 j_n = -(1 + e) * m * v_n;
  float v_t_mag = glm::length(v_t);
  float j_t_mag = min(mu * glm::length(j_n), m * v_t_mag);
  glm::vec3 j_t = -j_t_mag * glm::normalize(v_t);

  *dv = (j_n + j_t) / m;
  *dp = glm::normalize(p - origin) * (c_radius) + origin - p;
  return true;
}

void Sheet::collide(Shape *s) {
  for (size_t i = 0; i < nv; i++) {
    glm::vec3 p = vertices[i];
    glm::vec3 v = velocities[i];
    glm::vec3 dp;
    glm::vec3 dv;
    if (s->checkCollision(p, v, mass, &dp, &dv)) {
      std::cout << "collision" << std::endl;
      std::cout << "p: " << p.x << " " << p.y << " " << p.z << std::endl;
      std::cout << "v: " << v.x << " " << v.y << " " << v.z << std::endl;
      std::cout << "dp: " << dp.x << " " << dp.y << " " << dp.z << std::endl;
      std::cout << "dv: " << dv.x << " " << dv.y << " " << dv.z << std::endl;
      std::cout << std::endl;

      vertices[i] += dp;
      velocities[i] += dv;
    }
  }
}
