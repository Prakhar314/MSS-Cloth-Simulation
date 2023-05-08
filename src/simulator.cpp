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

void Simulator::addShape(Shape *shape) {
  shape->object = r.createObject();

  uint32_t nv, nn, nt;
  glm::vec3 *vertices, *normals;
  glm::ivec3 *triangles;
  shape->getObject(nv, vertices, nn, normals, nt, triangles);

  shape->vertexBuf = r.createVertexAttribs(shape->object, 0, nv, vertices);
  shape->normalBuf = r.createVertexAttribs(shape->object, 1, nn, normals);
  r.createTriangleIndices(shape->object, nt, triangles);

  shapes.push_back(shape);
}

void Simulator::update(float t) {
  for (int i = 0; i < shapes.size(); ++i) {
    Shape *shape = shapes[i];
    shape->update(t);
    uint32_t nv, nn, nt;
    glm::vec3 *vertices, *normals;
    glm::ivec3 *triangles;
    shape->getObject(nv, vertices, nn, normals, nt, triangles);
    r.updateVertexAttribs(shape->vertexBuf, nv, vertices);
    r.updateVertexAttribs(shape->normalBuf, nv, normals);
  }
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

  for (auto shape : shapes) {
    r.setupFilledFaces();
    r.setUniform(program, "objectColor", shape->color);
    r.drawObject(shape->object);

    r.setupWireFrame();
    r.setUniform(program, "objectColor", vec3(0.0f, 0.0f, 0.0f));
    r.drawObject(shape->object);
  }

  r.show();
}
