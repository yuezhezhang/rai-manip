#include "methods.h"
#include <Control/control.h>

struct PathFixer : PathMethod {
  PathFixer(ConfigurationProblem& _P, const ptr<PathResult>& _initialPath) : PathMethod(_P, _initialPath) {}

  virtual ptr<PathResult> run();
};

struct CtrlTarget_PD;

struct PathFeedbackFollower : PathMethod {
  CtrlSet ctrlSet;
  CtrlSolver ctrl;
  std::shared_ptr<CtrlObjective> qAccCost, qVelCost, coll, follow;

  PathFeedbackFollower(ConfigurationProblem& _P, const ptr<PathResult>& _initialPath, double maxStep);

  virtual ptr<PathResult> run();

  bool step_old(int verbose);
};
