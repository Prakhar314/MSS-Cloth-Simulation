#include "./sheet.hpp"

using namespace std;

void Sheet::setUsePBD(bool usePBD) { this->usePBD = usePBD; }
void Sheet::setSelfCollisions(bool selfCollisions) {
  this->selfCollisions = selfCollisions;
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
  if (usePBD) {
    oldPositions = new glm::vec3[nv];
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
      springs[springIndex++] = Spring(i * (n + 1) + j, (i + 1) * (n + 1) + j,
                                      spacing, usePBD ? -1 : ksStr, kdStr);
    }
  }
  for (size_t i = 0; i < m + 1; i++) {
    for (size_t j = 0; j < n; j++) {
      springs[springIndex++] = Spring(i * (n + 1) + j, i * (n + 1) + j + 1,
                                      spacing, usePBD ? -1 : ksStr, kdStr);
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

void Sheet::selfCollide(bool updateVel) {

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
              if (updateVel) {
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
}

void Sheet::update(float t, const vector<Shape *> &shapes) {
  if (t < 0.00001f)
    return;
  if (usePBD) {
    std::copy(vertices, vertices + nv, oldPositions);
  }
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
  if (usePBD) {
    for (int i = 0; i < 20; i++) {
      // structural constraints
      for (int j = 0; j < nSprings; j++) {
        if (springs[j].ks > 0) {
          continue;
        }
        Spring s = springs[j];
        glm::vec3 p1 = vertices[s.i];
        glm::vec3 p2 = vertices[s.j];
        glm::vec3 dp = p1 - p2;
        float l = glm::length(dp);
        glm::vec3 x_hat = dp / l;
        glm::vec3 imp = (s.restLength - l) * x_hat;

        bool i_fixed = false;
        bool j_fixed = false;
        for (auto &k : fixedParticles) {
          if (s.i == k) {
            i_fixed = true;
          }
          if (s.j == k) {
            j_fixed = true;
          }
        }
        if (i_fixed) {
          vertices[s.j] -= imp;
        } else if (j_fixed) {
          vertices[s.i] += imp;
        } else {
          vertices[s.i] += imp / 2.0f;
          vertices[s.j] -= imp / 2.0f;
        }
      }
      // self collisions
      if (selfCollisions)
        selfCollide(false);
      // other collisions
      // for (int j = 0; j < shapes.size(); ++j) {
      //   Shape *other = shapes[j];
      //   if (!other->isSheet()) {
      //     collide(other, false);
      //   }
      // }
    }
    // update velocities
    for (size_t i = 0; i < nv; i++) {
      velocities[i] = (vertices[i] - oldPositions[i]) / t;
    }
  }
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

void Sheet::collide(Shape *s, bool updateVel) {
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
      vertices[i] += glm::vec3(shape_to_sheet_t * glm::vec4(dp, 0.0f));
      if (updateVel)
        velocities[i] += glm::vec3(shape_to_sheet_t * glm::vec4(dv, 0.0f));
    }
  }
}

