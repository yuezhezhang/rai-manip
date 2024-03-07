#pragma once

#include "SkeletonProblem.h"

struct SkeletonSampler {
  SkeletonSampler(SkeletonProblem& P) : P(P) {}

  struct Result{
    arr frameStates;
    arr jointStates;
    rai::Graph report;
    double sos_sumOfSqr;
    double eq_sumOfAbs;
    double ineq_sumOfPos;
  };
  ptr<Result> sample(double noise=.01);

private:
  SkeletonProblem& P;
};

struct TimedSkeletonSampler: public SkeletonSampler{
  TimedSkeletonSampler(TimedSkeletonProblem &_TP): SkeletonSampler(_TP), TP(_TP){}
  ptr<Result> sample(const arr &ts, double noise=.01);
private:
  TimedSkeletonProblem &TP;
};
