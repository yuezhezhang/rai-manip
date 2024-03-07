#pragma once

#include "ConfigurationProblem.h"
#include "Animation.h"

#include <KOMO/komo.h>


struct SkeletonProblem : ConfigurationProblem {
  KOMO komo; //here (ab)used only to convert the skeleton into constraints, and store the constraints (in komo.objs)

  intA fixedConfigurations;
  arrA qForFixedConfigurations;
  rai::Array<FrameL> fixedFrames;

  SkeletonProblem(const rai::Configuration& _C, const bool _computeCollisions=true);
  SkeletonProblem(const rai::Configuration& _C, const Skeleton &S, const bool _computeCollisions=true);

  SkeletonProblem& setT(const uint T);
  SkeletonProblem& setSkeleton(const Skeleton& S){ komo.setSkeleton(S); return *this; }
  SkeletonProblem& setFixedConfigurations(const arrA& _qForFixedConfigurations);

  ///query the when setting configuration s (in animation time configurationTime(s)) to joint state x
  shared_ptr<QueryResult> query(const arr& x, const uint s);
};

struct TimedSkeletonProblem : public SkeletonProblem {
  rai::Animation A;

  TimedSkeletonProblem(const rai::Configuration& _C, const Skeleton &S, const rai::Animation &_A): SkeletonProblem(_C, S, true), A(_A){};

  TimedSkeletonProblem& setAnimation(const rai::Animation& _A){ A=_A; return *this; }

  shared_ptr<QueryResult> query(const arr& x, const arr &t, const uint s);
};
