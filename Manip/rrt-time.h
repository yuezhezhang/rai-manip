#pragma once

#include "spdlog/spdlog.h"

#include <algorithm>

#include <Gui/opengl.h>
#include <GL/gl.h>

#include <PlanningSubroutines/ConfigurationProblem.h>
#include "timedPath.h"

struct Node{
  arr q;
  arr x;
  double t;
  double cost = -1;

  Node *parent = nullptr;
  std::vector<Node*> children;

  Node(arr _q, double _t): q(_q), t(_t) {};

  void addChild(Node* node){
    children.push_back(node);
  }

  void print() const{
    std::cout << "t: " << t << std::endl;
    std::cout << "q: " << q << std::endl;
  }

  void updateCost(const double newCost){
    const double diff = cost - newCost;

    for (auto c: children){
      c->updateCost(c->cost - diff);
    }
    cost = newCost;
  }
};


struct Tree: GLDrawer {
  double comp_time_us{0.};

  std::vector<Node*> nodes;
  bool reverse;

  std::function<double (const Node&, const Node&)> distance;

  Tree(std::function<double (const Node&, const Node&)> _distance, bool _reverse=false):reverse(_reverse), distance(_distance){};
  ~Tree(){
    //std::cout << "deleting stuff" << std::endl;
    for (auto n: nodes){

      delete n;
    }
  }

  Node * getNearest(const Node &target, const bool print = false){
    const auto start_time = std::chrono::high_resolution_clock::now();

    Node *close = nullptr;
    double min_dist = 1e7;
    double d;
    for (Node* n: nodes){
      if (reverse) {d = distance(target, *n);}
      else {d = distance(*n, target);}

      if (print){
        target.print();
        n->print();
        std::cout << "dist: " << d << std::endl << std::endl;
      }

      if (d < min_dist){
        close = n;
        min_dist = d;
      }
    }

    if (print){
      std::cout << std::endl;
    }

    const auto end_time = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count();

    comp_time_us += duration;

    return close;
  };

  std::vector<Node*> getNear(const Node &start, const double radius){
    const auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<Node*> closeNodes;
    double d;
    for (Node* n: nodes){
      if (n == &start){
        continue;
      }

      if (!reverse) {d = distance(start, *n);}
      else {d = distance(*n, start);}

      if (d < radius){
        closeNodes.push_back(n);
      }
    }

    const auto end_time = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count();

    comp_time_us += duration;

    return closeNodes;
  }

  bool rewire(Node *nNode, TimedConfigurationProblem &TP){
    bool rewired = false;
    auto closeNodes = getNear(*nNode, 1);
    for (Node *n: closeNodes){
      //std::cout << reverse << std::endl;
      //std::cout << nNode->q << " " << nNode->t << std::endl;
      //std::cout << n->q << " " << n->t << std::endl;

      if (n->parent == nNode){
        continue;
      }

      double d;
      if (!reverse) {d = distance(*nNode, *n);}
      else {d = distance(*n, *nNode);}

      //std::cout << d << std::endl;

      if (n->cost > nNode->cost + d){
        //std::cout << n->cost << " " << nNode->cost + d << std::endl;
        if (!TP.checkEdge(nNode->q, n->t, n->q, n->t)) {continue;}

        n->updateCost(nNode->cost + d);

        std::vector<Node*> &childVec = n->parent->children;
        childVec.erase(std::remove(childVec.begin(), childVec.end(), n), childVec.end());
        n->parent = nNode;
        nNode->addChild(n);

        rewired = true;
      }
    }
    return rewired;
  }

  void addNode(Node* n){
    nodes.push_back(n);
  }

  void glDraw(OpenGL &gl){
    if (reverse){
      glColor(0., 0., 0.);
    }
    else{
      glColor(1., 1., 1.);
    }
    glLineWidth(2.f);
    glBegin(GL_LINES);
    for(uint i=0;i<nodes.size();i++){
      if (nodes[i]->parent){
        glVertex3dv(&(nodes[i]->parent->x(0)));
        glVertex3dv(&(nodes[i]->x(0)));
      }
    }
    glEnd();
    glLineWidth(1.f);
  }
};

using TimedGoalSampler = std::function<void (const double, arr&)>;
struct PathFinder_RRT_Time{
  TimedConfigurationProblem &TP;

  double vmax = .1;
  double lambda = .9;
  double step_time = 1;

  double goalSampleProbability = 0.9;//0.9

  uint maxInitialSamples = 50;

  uint maxIter = 1000;
  bool verbose = false;
  bool disp = false;

  double tol = 1e-3;
  bool connect = true;
  bool optimize = false;
  uint conv_iter = 500;
  bool sampleAdditionalNodes = false;

  bool informed_sampling = false;

  std::vector<bool> periodicDimensions;

  double edge_checking_time_us{0.};
  double nn_time_us{0.};

  FrameL prePlannedFrames;
  uint tPrePlanned;

  PathFinder_RRT_Time(TimedConfigurationProblem &_TP) :TP(_TP){
    delta_buffer = arr(TP.C.getJointState().N);
    periodicDimensions = std::vector<bool>(TP.C.getJointState().N, false);

    for (auto *j: TP.C.activeJoints){
      if (j->type == rai::JT_hingeX || j->type == rai::JT_hingeY|| j->type == rai::JT_hingeZ){
        periodicDimensions[j->qIndex] = true;
      }
    }
  };

  arr projectToManifold(const arr &q, const double t);

  TimedPath plan(const arr &q0, const double t0, const arr &qT, double tGoalLowerBound=0, double tGoalUpperBound=-1);
  TimedPath plan(const arr &q0, const double t0, const TimedGoalSampler gs, double tGoalLowerBound=0, double tGoalUpperBound=-1);

  arr delta_buffer;
  double q_metric(const arr& d) const;
  double distance(const Node &n1, const Node &n2);

  Node* extend(Tree* tree, const Node &goal, const bool connect = false);

  arr getDelta(const arr &p1, const arr &p2);
  Node* steer(const Node &start, const Node &goal, bool reverse);

  TimedPath extractPath(Tree* t1, Tree* t2, Node* leafNode1, Node* leafNode2){
    //std::cout << "cost: " << leafNode1->cost + leafNode2->cost << std::endl;

    // traverse over both sides and export result to time, path
    std::vector<Node*> np;
    {
      Node* n = leafNode1;
      while(!!n->parent){
        np.push_back(n);
        n = n->parent;
      }
      np.push_back(n);

      std::reverse(np.begin(), np.end());
    }

    {
      Node* n = leafNode2;
      while(n->parent){
        np.push_back(n);
        n = n->parent;
      }
      np.push_back(n);
    }

    if (t1->reverse){
      std::reverse(np.begin(), np.end());
    }

    arr path;
    arr time;

    path.resize(0, leafNode1->q.N);
    for (Node* n: np){
      time.append(n->t);
      path.append(n->q);
    }

    return TimedPath(path, time);
  }
};
