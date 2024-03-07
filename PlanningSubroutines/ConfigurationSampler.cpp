#include "ConfigurationSampler.h"

#include <Optim/constrained.h>

ptr<ConfigurationSampler::Result> ConfigurationSampler::run(){

  GoalStateProblem p_goal(p);

  p_goal.scaleCollisions=0.;
  arr x=p.q0;
  {
    arr dual;
    OptConstrained opt(x, dual, p_goal);
    opt.run();
  }
  p_goal.scaleCollisions=1e1;
  {
    arr dual;
    OptConstrained opt(x, dual, p_goal);
    opt.run();
  }

  //final query
  auto r = make_shared<Result>();
  r->qr = p.query(x);
  r->feasible = r->qr->isGoal && r->qr->isFeasible;
  r->goal = x;
  LOG(0) <<*r->qr;

  return r;
}

ptr<ConfigurationSampler::Result> TimedConfigurationSampler::run(const double t){
  tp.A.setToTime(tp.C, t);

  auto res = ConfigurationSampler::run();
  return res;
}
