#include "pathFixer.h"

#include <Algo/spline.h>
#include <Control/control.h>
#include <Kin/viewer.h>

ptr<PathResult> PathFixer::run(){
  bool feasible=true;
  arr path = initialPath->path;
  for(uint t=0;t<path.d0;t++){
//    cout <<"step " <<t;
    //bool f = makePoseFeasible(path[t](), P);
    //feasible &= f;
  }

  auto r = make_shared<PathResult>(*initialPath);
  r->path = path;
  r->feasible = feasible;
  return r;
}

double phaseVelocityControl(const arr& x_real, const arr& x_ref, const arr& v_ref, double upper = .05, double lower = .1){
  arr v = v_ref;
  double vel = length(v);

  double phaseVelocity = 1.;
  if(vel > 1e-8){
    v /= vel;
    double backlog = scalarProduct(v, x_ref - x_real);
    //    cout <<"spatial backlog: " <<backlog <<" \ttime backlog: " <<backlog/vel;
    double b = backlog/vel;
    if(b>upper){
      if(b>lower){
        phaseVelocity = 0.;
      }else{
        phaseVelocity = (b-lower)/(upper-lower);
      }
    }
    //    cout <<" \tphaseScale: " <<phaseScale <<endl;
  }

  return phaseVelocity;
}

/* TODO:
 * move getCollisions(y, J, FrameL&, margin) to rai::Configuration
 * negDistance feature can take larger FrameL!
 * have a basic NewtonIK solve problem of "closest to ref, but no collision"
 * have a KOMO IK solve that?
 */

PathFeedbackFollower::PathFeedbackFollower(ConfigurationProblem& _P, const ptr<PathResult>& _initialPath, double maxStep)
  : PathMethod(_P, _initialPath), ctrl(P.C, .1, 2) {

  qVelCost = ctrlSet.addObjective(make_feature(FS_qItself, {}, P.C, {1e-4}, {}, 1), OT_sos);
  qAccCost = ctrlSet.addObjective(make_feature(FS_qItself, {}, P.C, {1e-4}, {}, 2), OT_sos);
  coll = ctrlSet.addObjective(make_feature(FS_accumulatedCollisions, {}, P.C, {1e2}), OT_eq);

  follow = ctrlSet.addObjective(make_feature(FS_qItself, {}, P.C, {1e0}), OT_sos);
  follow ->setRef(make_shared<CtrlTarget_PathCarrot>(initialPath->path, maxStep, 1.));

  ctrl.set(ctrlSet);
}

#if 0
bool PathFeedbackFollower::step_old(int verbose){
  arr x_ref = S->eval(phase);
  arr v_ref = S->eval(phase, 1);

  //PD on pose
  pd->setTarget(x_ref, v_ref);
  if(phase<1.){
    pd->setTimeScale( .1 ); //5./double(S.times.N) ); //10.*tau);
  }else{
    pd->setTimeScale( .1 ); //5./double(S.times.N) ); //10.*tau);
  }
  NIY //  auto status = pd->step(x, tau, x);
  x = pd->y_ref;
  v = pd->v_ref;

  std::shared_ptr<QueryResult> qr = P.query(x);
  arr x0=x;
  uint count=0;
  while(false && qr->coll_y.N){
    //      qr->writeDetails(cout, P, 10.);
#if 0
    arr delta = qr->getBackwardStep(1.1, .1);
    double l = length(delta);
    double d = .1;
    if(l>d) delta *= d/l;
    x += delta;
    if(l<d) break;
    qr = P.query(x);
#else
    arr y,J;
    double margin=.05;
    qr->getViolatedContacts(y, J, margin);
    if(!y.N) break;
    y -= margin;

    if(sumOfAbs(y)<1e-3) break;

    arr Jinv = pseudoInverse(J, NoArr, 1e-4);
    arr delta = Jinv * (-.5 * y);
    delta += (eye(x.N) - Jinv * J) * (x0 - x); //null step (zero point regularization of Newton step)

    x += delta;

    qr = P.query(x);

    if(length(delta)<1e-3) break;
    if(count++ >10) break;
#endif

  }

  if(verbose>0){
    LOG(0) <<" phase:" <<phase <<" err:" <<length(x-x_ref) <<" coll:" <<min(cat({0.},qr->coll_y)) <<" qr:" <<*qr <<endl;
  }

  if(phase<1.){
    phase += tau; //phaseVelocityControl(x, x_ref, v_ref, .05, .1) * tau;
  }

  return (phase>=1./* && status==AS_converged*/);
}

bool PathFeedbackFollower::step(int verbose){
  arr x_ref = S->eval(phase);
  arr v_ref = S->eval(phase, 1);


  std::shared_ptr<QueryResult> qr = P.query(x);
  arr x0=x;
  uint count=0;
  while(false && qr->coll_y.N){
    //      qr->writeDetails(cout, P, 10.);
#if 0
    arr delta = qr->getBackwardStep(1.1, .1);
    double l = length(delta);
    double d = .1;
    if(l>d) delta *= d/l;
    x += delta;
    if(l<d) break;
    qr = P.query(x);
#else
    arr y,J;
    double margin=.05;
    qr->getViolatedContacts(y, J, margin);
    if(!y.N) break;
    y -= margin;

    if(sumOfAbs(y)<1e-3) break;

    arr Jinv = pseudoInverse(J, NoArr, 1e-4);
    arr delta = Jinv * (-.5 * y);
    delta += (eye(x.N) - Jinv * J) * (x0 - x); //null step (zero point regularization of Newton step)

    x += delta;

    qr = P.query(x);

    if(length(delta)<1e-3) break;
    if(count++ >10) break;
#endif

  }

  if(verbose>0){
    LOG(0) <<" phase:" <<phase <<" err:" <<length(x-x_ref) <<" coll:" <<min(cat({0.},qr->coll_y)) <<" qr:" <<*qr <<endl;
  }

  if(phase<1.){
    phase += tau; //phaseVelocityControl(x, x_ref, v_ref, .05, .1) * tau;
  }

  return (phase>=1./* && status==AS_converged*/);
}
#endif

ptr<PathResult> PathFeedbackFollower::run(){

  arr q = initialPath->path[0];
  P.C.setJointState(q);
  ctrl.komo.setConfiguration(-2, q);
  ctrl.komo.setConfiguration(-1, q);

  rai::ConfigurationViewer V;
  V.setConfiguration(P.C);

  int verbose=2;
  arr path(0u,q.N);
  path.append(q);

  for(uint t=0;t<1000;t++){
    ctrl.update(P.C);
    q = ctrl.solve();
    P.C.setJointState(q);
    V.setConfiguration(P.C, STRING("following t:" <<t), false);
    rai::wait(.05);
    path.append(q);
    if(follow->status==AS_converged) break;
  }

  auto r = make_shared<PathResult>(path);
  return r;
}
