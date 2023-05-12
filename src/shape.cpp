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

void Sheet::setDimensions(uint32_t width, uint32_t height, float spacing) {
  this->width = width;
  this->height = height;
  this->spacing = spacing;
}

void Sheet::fixParticle(uint32_t i, uint32_t j) {
  fixedParticles.push_back(i * (height + 1) + j);
}

glm::vec3 Sheet::getCenter() { return glm::vec3(0, 0, 0); }

void Sheet::init() {
  assert(width > 0 && height > 0 && spacing > 0.0f);
  uint32_t m = width;
  uint32_t n = height;
  nv = (m + 1) * (n + 1);
  nn = nv;
  nt = 2 * m * n;

  initBuffers();
  if (selfCollisions) {
    bins = new vector<uint16_t>[m * n * max(m, n)];
  }

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

void Sheet::selfCollide() {

  glm::vec3 max_p = vertices[0];
  glm::vec3 min_p = vertices[0];

  // find number of bins
  for (int i = 0; i < nv; i++) {
    max_p = glm::max(max_p, vertices[i]);
    min_p = glm::min(min_p, vertices[i]);
  }

  int bin_x = ceil((max_p.x - min_p.x) / spacing);
  int bin_y = ceil((max_p.y - min_p.y) / spacing);
  int bin_z = ceil((max_p.z - min_p.z) / spacing);

  // clear bins
  for (int i = 0; i < bin_x * bin_y * bin_z; i++) {
    bins[i].clear();
  }

  // fill bins
  for (uint16_t i = 0; i < nv; i++) {
    uint16_t x = floor((vertices[i].x - min_p.x) / spacing);
    uint16_t y = floor((vertices[i].y - min_p.y) / spacing);
    uint16_t z = floor((vertices[i].z - min_p.z) / spacing);
    bins[x * bin_y * bin_z + y * bin_z + z].push_back(i);
  }

  for (int i = 0; i < nv; ++i) {
    // get bin
    int x = floor((vertices[i].x - min_p.x) / spacing);
    int y = floor((vertices[i].y - min_p.y) / spacing);
    int z = floor((vertices[i].z - min_p.z) / spacing);

    // check all 3x3x3 bins
    for (int x_d = x - 1; x_d < x + 2; x_d++) {
      for (int y_d = y - 1; y_d < y + 2; y_d++) {
        for (int z_d = z - 1; z_d < z + 2; z_d++) {
          if (x_d < 0 || x_d >= bin_x || y_d < 0 || y_d >= bin_y || z_d < 0 ||
              z_d >= bin_z) {
            continue;
          }
          int bin_idx = x_d * bin_y * bin_z + y_d * bin_z + z_d;
          for (int j : bins[bin_idx]) {
            if (i == j) {
              continue;
            }
            float l = glm::length(vertices[i] - vertices[j]);
            if (l < spacing) {
              bool i_fixed = false;
              bool j_fixed = false;

              for (auto &k : fixedParticles) {
                if (i == k) {
                  i_fixed = true;
                }
                if (j == k) {
                  j_fixed = true;
                }
              }

              glm::vec3 x_hat = (vertices[i] - vertices[j]) / l;
              glm::vec3 dp = (spacing - l) * x_hat;

              if (i_fixed) {
                vertices[j] -= dp;
              } else if (j_fixed) {
                vertices[i] += dp;
              } else {
                vertices[i] += dp / 2.0f;
                vertices[j] -= dp / 2.0f;
              }
              glm::vec3 rv = velocities[i] - velocities[j];
              glm::vec3 v_n = glm::dot(rv, x_hat) * x_hat;
              glm::vec3 v_t = rv - v_n;
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

              if (i_fixed) {
                velocities[j] -= (j_n + j_t) / mass;
              } else if (j_fixed) {
                velocities[i] += (j_n + j_t) / mass;
              } else {
                velocities[i] += (j_n + j_t) / mass / 2.0f;
                velocities[j] -= (j_n + j_t) / mass / 2.0f;
              }
            }
          }
        }
      }
    }
  }
}

void Sheet::update(float t, const vector<Shape *> &shapes) {
  // Compute the forces
  glm::vec3 gravity = glm::vec3(0, g, 0);
  fill_n(acc, nv, glm::inverse(this->transform) * glm::vec4(gravity, 0));
  for (size_t i = 0; i < nSprings; i++) {
    Spring s = springs[i];
    if (s.ks < 0)
      continue; // skip PBD constraints
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
  // fixed particles have zero acceleration
  for (auto &i : fixedParticles) {
    acc[i] = glm::vec3(0, 0, 0);
  }
  for (size_t i = 0; i < nv; i++) {
    velocities[i] += acc[i] * t;
    vertices[i] += velocities[i] * t;
    assert(glm::length(velocities[i]) < 1000);
  }

  // constraint resolution
  if (selfCollisions) {
    selfCollide();
  }
  for (int j = 0; j < shapes.size(); ++j) {
    Shape *other = shapes[j];
    if (!other->isSheet()) {
      collide(other);
    }
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

void Sheet::collide(Shape *s, bool updatePos, bool updateVel) {
  glm::mat4 sheet_t = getTransform();
  glm::mat4 sheet_inv_t = glm::inverse(sheet_t);
  glm::mat4 shape_t = s->getTransform();
  glm::mat4 shape_inv_t = glm::inverse(shape_t);
  glm::mat4 sheet_to_shape_t = shape_inv_t * sheet_t;
  glm::mat4 shape_to_sheet_t = sheet_inv_t * shape_t;
  for (size_t i = 0; i < nv; i++) {
    glm::vec3 p = sheet_to_shape_t * glm::vec4(vertices[i], 1.0f);
    glm::vec3 v = sheet_to_shape_t * glm::vec4(velocities[i], 0.0f);
    glm::vec3 dp;
    glm::vec3 dv;
    if (s->checkCollision(p, v, mass, &dp, &dv)) {
      if (updatePos)
        vertices[i] += glm::vec3(shape_to_sheet_t * glm::vec4(dp, 0.0f));
      if (updateVel)
        velocities[i] += glm::vec3(shape_to_sheet_t * glm::vec4(dv, 0.0f));
    }
  }
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
