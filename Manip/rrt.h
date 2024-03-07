#ifndef _RRT_H
#define _RRT_H

#include <stdlib.h>
#include <Algo/ann.h>
#include <Kin/kin.h>
#include <Kin/proxy.h>
#include <Gui/opengl.h>
#include <Kin/frame.h>

#include <PlanningSubroutines/ConfigurationProblem.h>

#include "methods.h"

/// just a data structure, no algorithms
struct RRT_SingleTree : GLDrawer {
  ANN ann;         //ann stores all points added to the tree in ann.X
  uintA parent;    //for each point we also store the index of the parent node
  rai::Array<ptr<QueryResult>> queries;
  arr disp3d;
  uint nearestID = UINT_MAX; //nearest node from the last 'getProposalToward' call!

  RRT_SingleTree(const arr& q0, const ptr<QueryResult>& q0_qr);

  //core method
  double getNearest(const arr& target);
  arr getProposalTowards(const arr& target, double stepsize);

  arr getNewSample(const arr& target, double stepsize, double p_sideStep, bool& isSideStep, const uint recursionDepth);

  //trivial
  uint add(const arr& q, uint parentID, const ptr<QueryResult>& _qr);

  //trivial access routines
  uint getParent(uint i){ return parent(i); }
  uint getNumberNodes(){ return ann.X.d0; }
  uint getDim(){ return ann.X.d1; }
  arr getNode(uint i){ return ann.X[i].copy(); }
  arr getLast(){ return ann.X[ann.X.d0-1].copy(); }
  arr getRandomNode(){ return ann.X[rnd(ann.X.d0)].copy(); }
  arr getPathFromNode(uint fromID);

  void glDraw(OpenGL &gl);

public:
  arr getSideStep(std::shared_ptr<QueryResult> qr);
};

//===========================================================================

class CollisionChecker{
  public: 
    CollisionChecker(const rai::Configuration &C) :cp(C){};

    virtual bool isValid(const arr &q){
      qr = cp.query(q);
      return qr->isFeasible;
    };

    bool isValidForRobot(const arr &q);
    bool isValidForAll(const arr &q);
    bool isValidForFrame(const arr &q);

    ConfigurationProblem cp;
    std::shared_ptr<QueryResult> qr;
};

class TorqueCollisionChecker: public CollisionChecker{
  public:
    TorqueCollisionChecker(const rai::Configuration &C, const std::string &ee, const std::string &obj)
      : CollisionChecker(C), ee_(ee), obj_(obj){};

    virtual bool isValid(const arr &q){
      bool collFeasible = CollisionChecker::isValid(q);

      // geometry collsion check
      if (!collFeasible){
        return false;
      }

      // torque in endeffector collision check
      if (cp.C[STRING(obj_)]->parent->name == STRING(ee_)){
        const arr objPos = cp.C[STRING(obj_)]->getPosition();
        const arr eePos = cp.C[STRING(ee_)]->getPosition();

        const arr diff = objPos - eePos;
        const double l = length(diff({0,1}));

        //std::cout << l << " " << diff({0,1}) << std::endl;

        if (l > limit_){
          return false;
        }
      }

      return true;
    }

    double limit_{.3};
    std::string ee_;
    std::string obj_;
};

///algorithms
struct PathFinder_RRT : PathFinder {
  double stepsize;
  uint verbose;
  bool intermediateCheck{true};

  uint maxIter{50000};

  uint n_backStep=0, n_backStepGood=0, n_sideStep=0, n_sideStepGood=0, n_forwardStep=0, n_forwardStepGood=0, n_rndStep=0, n_rndStepGood=0;

  std::shared_ptr<CollisionChecker> cc_;

  bool checkConnection(const arr &start, const arr &end, const double res=0.01);

  PathFinder_RRT(ConfigurationProblem& _P, const arr& goals, double _stepsize = .2, uint _verbose=0, bool _intermediateCheck=true) : PathFinder(_P, goals), stepsize(_stepsize), verbose(_verbose), intermediateCheck(_intermediateCheck) {
    cc_ = std::make_shared<CollisionChecker>(_P.C);
    //cc_ = std::make_shared<TorqueCollisionChecker>(_P.C, "a0_gripper", "table");
    //cc_ = std::make_shared<TorqueCollisionChecker>(_P.C, "a0_gripper", "stick");
  }

  void planForward(const arr& qT);
  arr planConnect(const arr& qT, double p_forwardStep=.5, double p_sideStep=.0, double p_backwardStep=.0); //default numbers: equivalent to standard bidirect

  bool growTreeTowardsRandom(RRT_SingleTree& rrt);
  bool growTreeToTree(RRT_SingleTree& rrt_A, RRT_SingleTree& rrt_B, double p_forwardStep, double p_sideStep, double p_backwardStep);

  virtual ptr<PathResult> run(double timeBudget=1.);

private:
  rai::Configuration DISP;
};

void revertPath(arr& path);

#endif
