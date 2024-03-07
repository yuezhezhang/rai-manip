#include "SkeletonSampler.h"

#include <Optim/constrained.h>

ptr<SkeletonSampler::Result> SkeletonSampler::sample(double noise){
  OptOptions opt;
  opt.stopTolerance = 1e-4;
  opt.stopFTolerance = 1e-4;
  opt.stopGTolerance = 1e-4;

  P.komo.verbose = 0;

  // set the fixed keyframes again
  if (P.qForFixedConfigurations.d0 != 0){
    for (uint t=0; t<P.komo.timeSlices.d0-P.komo.k_order; ++t){
      if (P.qForFixedConfigurations(t).d0 == 0){continue;}
      P.komo.pathConfig.setJointState(P.qForFixedConfigurations(t), P.fixedFrames(t));
    }
  }

  //P.komo.displayPath("pre");

  P.komo.reset();
  P.komo.optimize(noise, opt);

  //P.komo.displayPath("post");

  auto r = make_shared<Result>();

  r->frameStates = P.komo.getPath_frames();
  r->jointStates = P.komo.getPath(P.C.getJointIDs());

  r->report = P.komo.getReport(false);
  r->sos_sumOfSqr = r->report.get<double>("sos");
  r->eq_sumOfAbs = r->report.get<double>("eq");
  r->ineq_sumOfPos = r->report.get<double>("ineq");

  return r;
}

ptr<SkeletonSampler::Result> TimedSkeletonSampler::sample(const arr& ts, double noise){
  rai::Configuration C;
  C.copy(TP.komo.world);

  for (uint i=0; i<ts.d0; ++i){
    const uint t = ts(i);
    TP.A.setToTime(C, t);

    const FrameL F = TP.komo.timeSlices[i+TP.komo.k_order];
    TP.komo.pathConfig.setFrameState(C.getFrameState(), F);
  }

  return SkeletonSampler::sample(noise);
}
