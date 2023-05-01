#include "simulator.hpp"

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

void Simulator::init(const std::string &filename, const int width,
                     const int height) {
  r.initialize(filename, width, height);
  camCtl.initialize(width, height);
  camCtl.camera.setCameraView(vec3(0.5, -0.5, 1.5), vec3(0.5, -0.5, 0.0),
                              vec3(0.0, 1.0, 0.0));
  program = r.createShaderProgram(r.vsPhongShading(), r.fsPhongShading());
}

void Simulator::createScene() {
  vertices = new vec3[nv];
  normals = new vec3[nv];
  triangles = new ivec3[nt];

  object = r.createObject();
  vertices[0] = vec3(0, 0, 1);
  vertices[1] = vec3(1, 0, 1);
  vertices[2] = vec3(1, 0, 0);
  vertices[3] = vec3(0, 0, 0);
  vertexBuf = r.createVertexAttribs(object, 0, nv, vertices);
  normals[0] = vec3(0, 0, 1);
  normals[1] = vec3(0, 0, 1);
  normals[2] = vec3(0, 0, 1);
  normals[3] = vec3(0, 0, 1);
  normalBuf = r.createVertexAttribs(object, 1, nv, normals);
  triangles[0] = ivec3(0, 1, 2);
  triangles[1] = ivec3(0, 2, 3);
  r.createTriangleIndices(object, nt, triangles);
}

void Simulator::update(float t) {
  float freq = 2, amp = 1;
  float phase0 = 0, phase1 = 1;
  float theta0 = amp * cos(freq * t + phase0),
        theta1 = amp * cos(freq * t + phase1);
  vertices[0] = vec3(0, -cos(theta0), sin(theta0));
  vertices[1] = vec3(1, -cos(theta1), sin(theta1));
  r.updateVertexAttribs(vertexBuf, nv, vertices);
  normals[0] = glm::normalize(
      glm::cross(vertices[1] - vertices[0], vertices[3] - vertices[0]));
  normals[1] = glm::normalize(
      glm::cross(vertices[2] - vertices[1], vertices[0] - vertices[1]));
  normals[2] = glm::normalize(
      glm::cross(vertices[3] - vertices[2], vertices[1] - vertices[2]));
  normals[3] = glm::normalize(
      glm::cross(vertices[0] - vertices[3], vertices[2] - vertices[3]));
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
