#include "../src/container.hpp"
#include <glm/fwd.hpp>
#include <glm/gtx/transform.hpp>

int main() {
  int width = 640, height = 480;
  Simulator *sim = new Simulator();
  sim->init("Animation", width, height);

  glm::mat4 m = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -0.0f, 0.0f));
  Sheet *sheet = new Sheet();
  sheet->setDimensions(10, 10, 0.1f);
  sheet->init();
  sheet->setTransform(m);
  sim->addShape(sheet);

  m = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -0.9f, 0.0f));
  Plane *p = new Plane();
  p->setDimensions(glm::vec3(0, 1, 0));
  p->init();
  p->setTransform(m);
  p->color = glm::vec3(0.5f, 0.5f, 0.5f);
  sim->addShape(p);

  ShapeContainer *chest_c = new ShapeContainer();
  // chest_c->setScale(glm::vec3(0.5f, 0.5f, 0.5f));
  chest_c->setPosition(glm::vec3(0.5f, -0.5f, 0.0f));
  Cylinder *chest = new Cylinder();
  chest->setDimensions(0.3f, 0.8f, 20, 10);
  chest->init();
  chest->color = glm::vec3(0.0f, 0.0f, 1.0f);
  chest_c->addShape(chest, sim);

  ShapeContainer *head_c = new ShapeContainer();
  head_c->setPosition(glm::vec3(0.0f, 0.7f, 0.0f));
  Sphere *head = new Sphere();
  head->setDimensions(0.35f, 20, 20);
  head->init();
  head->color = glm::vec3(0.0f, 0.0f, 1.0f);
  head_c->addShape(head, sim);
  chest_c->addChild(head_c);

  ShapeContainer *arm_l_c = new ShapeContainer();
  arm_l_c->setPosition(glm::vec3(0.45f, 0.0f, 0.0f));
  arm_l_c->setRecenter(glm::vec3(0.25f, 0.0f, 0.0f));
  Cylinder *arm_l = new Cylinder();
  arm_l->setDimensions(0.2f, 0.5f, 20, 10);
  arm_l->init();
  arm_l->color = glm::vec3(0.2f, 0.0f, 0.8f);
  arm_l_c->addShape(arm_l, sim);

  ShapeContainer *larm_l_c = new ShapeContainer();
  larm_l_c->setPosition(glm::vec3(0.3f, 0.0f, 0.0f));
  larm_l_c->setRecenter(glm::vec3(0.25f, 0.0f, 0.0f));
  Cylinder *larm_l = new Cylinder();
  larm_l->setDimensions(0.2f, 0.5f, 20, 10);
  larm_l->init();
  larm_l->color = glm::vec3(0.4f, 0.0f, 0.6f);
  larm_l_c->addShape(larm_l, sim);
  arm_l_c->addChild(larm_l_c);

  ShapeContainer *arm_r_c = new ShapeContainer();
  arm_r_c->setPosition(glm::vec3(-0.45f, 0.0f, 0.0f));
  arm_r_c->setRecenter(glm::vec3(-0.25f, 0.0f, 0.0f));
  arm_r_c->setRotation(15.0f, glm::vec3(0.0f, 0.0f, -1.0f));
  Cylinder *arm_r = new Cylinder();
  arm_r->setDimensions(0.2f, 0.5f, 20, 10);
  arm_r->init();
  arm_r->color = glm::vec3(0.2f, 0.0f, 0.8f);
  arm_r_c->addShape(arm_r, sim);

  ShapeContainer *larm_r_c = new ShapeContainer();
  larm_r_c->setPosition(glm::vec3(-0.3f, 0.0f, 0.0f));
  larm_r_c->setRecenter(glm::vec3(-0.25f, 0.0f, 0.0f));
  larm_r_c->setRotation(15.0f, glm::vec3(0.0f, 0.0f, -1.0f));
  Cylinder *larm_r = new Cylinder();
  larm_r->setDimensions(0.2f, 0.5f, 20, 10);
  larm_r->init();
  larm_r->color = glm::vec3(0.4f, 0.0f, 0.6f);
  larm_r_c->addShape(larm_r, sim);
  arm_r_c->addChild(larm_r_c);

  chest_c->addChild(arm_l_c);
  chest_c->addChild(arm_r_c);

  float ts = SDL_GetTicks64() / 1e3;
  while (!sim->shouldQuit()) {
    float tnew = SDL_GetTicks64() / 1e3;
    chest_c->update((tnew - ts) * 1.0f);
    sim->update((tnew - ts) * 1.0f);
    sim->render();
    ts = tnew;
  }

  delete sim;
}
