#include "SkeletonProblem.h"

SkeletonProblem::SkeletonProblem(const rai::Configuration& _C, const bool _computeCollisions)
  : ConfigurationProblem(_C, _computeCollisions){
  komo.setModel(_C, _computeCollisions);
}
SkeletonProblem::SkeletonProblem(const rai::Configuration& _C, const Skeleton &S, const bool _computeCollisions)
  : SkeletonProblem(_C, _computeCollisions){
  uint T = 0;
  for(const SkeletonEntry& s:S) {
    if (s.phase0 > T) {T = s.phase0;}
    if (s.phase1 > T) {T = s.phase1;}
  }

  // TODO: why does this have to be in this exact order?
  setT(T);
  setSkeleton(S);

  komo.verbose = 0;
}

SkeletonProblem& SkeletonProblem::setT(const uint T){
  komo.setTiming(T, 1, 5., 1);
  //komo.add_qControlObjective({}, 1, 1e-1);
  //komo.add_collision(true, .0, 1e1);
  komo.add_collision(true, .01, 1e0);
  komo.add_jointLimits(true, 0., 1e1);
  //komo.addSquaredQuaternionNorms();

  return *this;
}

SkeletonProblem& SkeletonProblem::setFixedConfigurations(const arrA& _qForFixedConfigurations){
  qForFixedConfigurations = _qForFixedConfigurations;

  FrameL activeConfigurations;
  for (uint t=0; t<komo.timeSlices.d0-komo.k_order; ++t){
    if (qForFixedConfigurations(t).d0 == 0){
      for(auto* f:komo.timeSlices[t + komo.k_order]) {
        if(f->joint && f->joint->active &&
            f->joint->type != rai::JT_rigid && f->joint->type != rai::JT_free){
          activeConfigurations.append(f);
        }
      }
      continue;
    }

    FrameL F;
    for(auto* f:komo.timeSlices[t + komo.k_order]) {
      if(f->joint && f->joint->active &&
          f->joint->type != rai::JT_rigid && f->joint->type != rai::JT_free){
        F.append(f);
      }
    }
    komo.pathConfig.setJointState(qForFixedConfigurations(t), F);
    fixedFrames.append(F);
  }
  // set configuratiosn to active
  komo.pathConfig.selectJoints(activeConfigurations);
  
  return *this;
}

shared_ptr<QueryResult> SkeletonProblem::query(const arr& x, const uint s){
  objectives.clear();

  intA configs={(int)s};
  for(shared_ptr<GroundedObjective>& o: komo.objs){
    if(o->configs==configs){
      objectives.append(make_shared<Objective>(o->feat, o->type));
    }
  }

  return ConfigurationProblem::query(x);
}

shared_ptr<QueryResult> TimedSkeletonProblem::query(const arr& x, const arr &ts, const uint s){
  rai::Configuration C;
  C.copy(komo.world);

  for (uint i=0; i<ts.d0; ++i){
    uint t = ts(i);
    A.setToTime(C, t);

    const FrameL F = komo.timeSlices[i+komo.k_order];
    komo.pathConfig.setFrameState(C.getFrameState(), F);
  }

  return SkeletonProblem::query(x, s);
};
