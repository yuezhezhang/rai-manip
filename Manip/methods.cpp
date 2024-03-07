#include "methods.h"

#include <Optim/constrained.h>
#include <Algo/spline.h>
#include <KOMO/pathTools.h>
#include <Kin/F_qFeatures.h>

ptr<GoalSampler::Result> GoalSampler::run(){

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

arrA conv_arr2arrA(const arr& x){
  CHECK_EQ(x.nd, 2, "");
  arrA y(x.d0);
  for(uint i=0;i<y.N;i++) y(i) = x[i];
  return y;
}

PathOptimizer::PathOptimizer(ConfigurationProblem& p, const arr& initialPath, double duration)
  : P(p) {

  P.C.setJointState(initialPath[0]);
  setModel(P.C, true);
  setTiming(1., 64, duration);

  for(auto o: P.objectives)  addObjective({1.}, o->feat, {}, o->type);

  //-- other default objectives
  add_qControlObjective({}, 2, 1.);
  add_collision(true, .0, 1e1);
  addObjective({1.}, FS_qItself, {}, OT_sos, {1e2}, {}, 1); //slow at end

  reset();
  initWithWaypoints(conv_arr2arrA(initialPath), initialPath.d0, false);
}

ptr<PathOptimizer::Return> PathOptimizer::run(double timeBudget){
  //-- optimize
//  verbose=6;
//  animateOptimization=1;
  optimize(false);

//  while(displayTrajectory(.1, false)){}

  auto r = make_shared<Return>();
  r->path = pathConfig.getJointState();
  r->path.reshape(T,-1);

  r->report = getReport(false);
  r->sos_sumOfSqr = r->report.get<double>("sos");
  r->eq_sumOfAbs = r->report.get<double>("eq");
  r->ineq_sumOfPos = r->report.get<double>("ineq");

  return r;
}

void makeSpline(rai::Spline& S, const arr& path, double initialDuration){
  arr x = path;
  x.prepend(path[0]);
  x.append(path[-1]);
  arr times(x.d0);
  for(uint i=0;i<x.d0;i++) times.elem(i) = initialDuration*double(i)/double(times.N-1);

  S.set(2, x, times);
}

ptr<PathResult> TimeOptimizer::run(double initialDuration) {
  rai::Spline S;
  makeSpline(S, path, initialDuration);

  double maxVel=0., maxAcc=0.;
  for(double t=-.1;t<1.1;t+=.002){
//    arr q = S.eval(t, 0);
//    arr vel = S.eval(t, 1);
    double vel = length(S.eval(t,1));
    double acc = length(S.eval(t,2));
      //    cout <<t <<' ' <<vel <<' ' <<q <<endl;
    rai::maxEq(maxVel, vel);
    rai::maxEq(maxAcc, acc);
  }

  double v_speedUp = 2./maxVel;
  double a_speedUp = sqrt(3./maxAcc);
  double speedUp = rai::MIN(v_speedUp, a_speedUp);

  cout <<"maxVel:  " <<maxVel <<endl;
  cout <<"maxAcc:  " <<maxAcc <<endl;
  cout <<"speedUp: " <<speedUp <<endl;

  auto r = make_shared<PathResult>(path);
  r->duration = initialDuration/speedUp;
  std::cout << "duration: " << r->duration << std::endl;

  return r;
}

MultiAgentKOMO::MultiAgentKOMO(rai::Configuration& C, uint na, bool computeCollisions) : KOMO(){
  setModel(C, computeCollisions);
  setTiming(0., 64, 10., 2);
  addSquaredQuaternionNorms();
  setupConfigurations2();
  agentBusyTil.resize(na).setZero();
}

void MultiAgentKOMO::runForSingleAgent(uint agent, ConfigurationProblem& p, const arr& initialPath, double duration){
  uint s0 = agentBusyTil(agent);
  double t0 = conv_step2time(s0-1, stepsPerPhase);
  if(duration<1.) duration=1.;
  double phase = duration/(tau*stepsPerPhase);
  cout <<"agent: " <<agent <<" startStep: " <<s0 <<" startTime: " <<t0 <<" duration: " <<duration <<" phase: " <<phase <<endl;

  //-- resample the path to fit with KOMO resolution
  arr q = path_resample(initialPath, duration/(tau*initialPath.d0));
  uint sFinal = s0 + q.d0; //(plus 1)
  double tFinal = conv_step2time(sFinal-1, stepsPerPhase);
  cout <<"resampled length: " <<q.d0 <<endl;

  //extend KOMO to take longer path
  if(T<sFinal){
    T = sFinal;
    setupConfigurations2();
  }
  //extend the path to fit longer KOMO
  while(s0+q.d0<T) q.append(q[-1]);

  //-- select optimization variables
  selectJointsBySubtrees({STRING('a'<<agent<<"_world")}, {t0, 100.});

  FILE("z.g") <<world <<endl;

  //-- set path
  q.reshape(-1);
  set_x2(q);

#if 1
  //-- add objectives from the PathProblem into KOMO
  clearObjectives();
  for(auto o: p.objectives)  addObjective({tFinal}, o->feat, {}, o->type);
  allAgentObjectives.append(objectives);


  //-- other default objectives
  //constrain start state to q0 (because we can't use a prefix!)
  arr q0 = getConfiguration_q(s0);
  addObjective({t0}, FS_qItself, {}, OT_eq, {1e1}, q0, 0, +1, +2);
  //squared accelerations
//  setSquaredQAccVelHoming(t0, tFinal, 1., 0., 0., +3);
  add_qControlObjective({t0, tFinal}, 2, 1., {}, +3);
  //collisions
  addObjective({t0, tFinal}, FS_accumulatedCollisions, {}, OT_eq, {1e1});
  //slow at end
//  addObjective({tFinal}, FS_qItself, {}, OT_sos, {1e2}, {}, 1);
  ptr<Objective> slow = addObjective({tFinal}, make_shared<F_qItself>(world.getCtrlFramesAndScale()), {}, OT_sos, {1e2}, NoArr, 1);
  allAgentObjectives.append(slow);

  //-- optimize
  verbose=5;
//  animateOptimization=1;
  optimize(0.);

  //return value.. not used yet
  q = x;
  q.reshape(-1, q0.N);
#endif

  agentBusyTil(agent) = sFinal;

  //  displayTrajectory(.1, false);
}

/*void MultiAgentKOMO::pastePickAndPlaceAndRun(uint agent, uint object, const arr& goalFrameState, ptr<PathResult> initPickPath, ptr<PathResult> initPlacePath){
  uint s0 = agentBusyTil(agent);
  double t0 = conv_step2time(s0-1, stepsPerPhase);
  if(initPickPath->duration<1.) initPickPath->duration=1.;
  if(initPlacePath->duration<1.) initPlacePath->duration=1.;
  double duration = initPickPath->duration + initPlacePath->duration;
  double phase = duration/(tau*stepsPerPhase);
  cout <<"agent: " <<agent <<" startStep: " <<s0 <<" startTime: " <<t0 <<" duration: " <<duration <<" phase: " <<phase <<endl;

  //-- resample the path to fit with KOMO resolution
  arr q_pick = path_resample(initPickPath->path, initPickPath->duration/(tau*initPickPath->path.d0));
  arr q_place = path_resample(initPlacePath->path, initPlacePath->duration/(tau*initPlacePath->path.d0));
  uint s1 = s0 + q_pick.d0;
  uint s2 = s1 + q_place.d0; //(plus 1)
  uint s3 = s2 + stepsPerPhase/2;
  double t1 = conv_step2time(s1-1, stepsPerPhase);
  double t2 = conv_step2time(s2-1, stepsPerPhase);
  double t3 = conv_step2time(s3-1, stepsPerPhase);
  cout <<"resampled length: " <<q_pick.d0 <<endl;
  cout <<"resampled length: " <<q_place.d0 <<endl;

  //extend KOMO to take longer path
  if(T<s3){
    T = s3;
    setupConfigurations2();
  }
  //extend the path to fit longer KOMO
  while(s0 + q_pick.d0 + q_place.d0<T) q_place.append(q_place[-1]);

  //-- select optimization variables
  selectJointsBySubtrees({STRING('a'<<agent<<"_world")}, {t0, 100.});

  FILE("z.g") <<world <<endl;

  //-- set path
  arr _x = cat(q_pick, q_place);
  _x.reshape(-1);
  set_x2(_x);

  clearObjectives();

  //add switches; because configurations are created already, we need to apply switchtes in retrospect!
  addSwitch_mode(SY_initial, SY_stable, t1, t2, nullptr, STRING('a'<<agent<<"_gripper"), STRING('o'<<object));
  addSwitch(t2, true, rai::JT_rigid, rai::SWInit_copy, "bugaBase", STRING('o'<<object));
  addObjective({t2}, FS_pose, {STRING('o'<<object)}, OT_eq, {1e0}, NoArr, 1, +0, +1);

  retrospectChangeJointType(s1-1, T, world.getFrame(STRING('o'<<object))->ID, rai::JT_rigid);
  retrospectApplySwitches2();

  //initialize object pose after t2 with goal pose
  {
    uint t = conv_time2step(t2, stepsPerPhase);
    uint oid = world.getFrame(STRING('o'<<object))->ID;
    for(;t<T;t++){
      timeSlices(t+k_order, oid)->set_X()->set(goalFrameState[object+1]);
//      getConfiguration_t(t).getFrame(STRING('o'<<object))->set_X()->set(goalFrameState[object+1]);
    }
  }

  //touch
  addObjective({t1}, FS_distance, { STRING('a'<<agent<<"_finger1"), STRING('o'<<object)}, OT_eq, {1e1}, {});
  addObjective({t1}, FS_distance, { STRING('a'<<agent<<"_finger2"), STRING('o'<<object)}, OT_eq, {1e1}, {});
  //goal -- covered by the t2 switch constraint!
//  addObjective({t2}, FS_position, {STRING('o'<<object)}, OT_eq, {1e1}, goalFrameState(object+1,{0,2}), 0, -1, -1 );
//  addObjective({t2}, FS_quaternion, {STRING('o'<<object)}, OT_eq, {1e1}, goalFrameState(object+1,{3,6}), 0, -1, -1 );

  allAgentObjectives.append(objectives);

  //-- other default objectives
  //constrain start state to q0 (because we can't use a prefix!)
  arr q0 = getConfiguration_q(s0);
  addObjective({t0}, FS_qItself, {}, OT_eq, {1e1}, q0, 0, +1, +2);
  //squared accelerations
//  setSquaredQAccVelHoming(t0, t2, 1., 0., 0., +3);
  add_qControlObjective({t0,-1.}, 2, 1., {}, +3, +0);
  //collisions
  addObjective({t0, -1.}, FS_accumulatedCollisions, {}, OT_eq, {1e1});
  //slow at end
  ptr<Objective> slow1 = addObjective({t1}, make_shared<F_qItself>(world.getCtrlFramesAndScale()), {}, OT_sos, {1e2}, NoArr, 1);
  allAgentObjectives.append(slow1);
  ptr<Objective> slow2 = addObjective({t2}, make_shared<F_qItself>(world.getCtrlFramesAndScale()), {}, OT_sos, {1e2}, NoArr, 1);
  allAgentObjectives.append(slow2);
  //move away
  addObjective({t2,-1.}, FS_distance, { STRING('a'<<agent<<"_arm1"), STRING('o'<<object)}, OT_sos, {1e1}, {-.1}, 2);
//  add_qControlObjective({t3,-1.}, 1, 1e1);


#if 1
  //-- optimize
  verbose=5;
//  animateOptimization=1;
  optimize(0.);

  //return value.. not used yet
//  _x = x;
//  _x.reshape(-1, q0.N);
#endif

  agentBusyTil(agent) = s2;

  //  displayTrajectory(.1, false);
}*/

void MultiAgentKOMO::runForAllAgents(){

  //-- select ALL variables
  selectJointsBySubtrees({}, {}, true);

  FILE("z.g") <<world <<endl;
  world.report();

  //-- add objectives from the PathProblem into KOMO
  clearObjectives();
  //copy all agent goal objectives
  objectives = allAgentObjectives;

  //-- other default objectives
  //squared accelerations
  add_qControlObjective({}, 2, 1.);
  //collisions
  addObjective({}, FS_accumulatedCollisions, {}, OT_eq, {1e1});
  //slow at end
//  addObjective({tFinal}, FS_qItself, {}, OT_sos, {1e2}, {}, 1);

  //-- optimize
  verbose=5;
//  animateOptimization=1;
  optimize(0.);

//  displayTrajectory(.1, false);
}
