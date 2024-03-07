#include "rrt.h"

#include <GL/gl.h>
#include <Kin/viewer.h>

#include <chrono>

RRT_SingleTree::RRT_SingleTree(const arr& q0, const ptr<QueryResult>& q0_qr){
  add(q0, 0, q0_qr);
}

uint RRT_SingleTree::add(const arr& q, uint parentID, const ptr<QueryResult>& _qr){
  ann.append(q);
  parent.append(parentID);
  queries.append(_qr);
  disp3d.append(_qr->disp3d);
  disp3d.reshape(-1,3);

  CHECK_EQ(parent.N, ann.X.d0, "");
  CHECK_EQ(queries.N, ann.X.d0, "");
  CHECK_EQ(disp3d.d0, ann.X.d0, "");
  return parent.N-1;
}

double RRT_SingleTree::getNearest(const arr& target){
  //find NN
  nearestID = ann.getNN(target);
  return length(target - ann.X[nearestID]);
}

arr RRT_SingleTree::getProposalTowards(const arr& target, double stepsize){
  //find NN
  nearestID = ann.getNN(target);

  //compute default step
  arr delta = target - ann.X[nearestID]; //difference vector between q and nearest neighbor
  double dist = length(delta);
  if(dist>stepsize)  delta *= stepsize/dist;

  return getNode(nearestID) + delta;
}

arr RRT_SingleTree::getNewSample(const arr& target, double stepsize, double p_sideStep, bool& isSideStep, const uint recursionDepth){
  //find NN
  nearestID = ann.getNN(target);
  std::shared_ptr<QueryResult> qr = queries(nearestID);

  //compute default step
  arr delta = target - getNode(nearestID);
  double dist = length(delta);
  if(dist>stepsize) delta *= stepsize/dist;

  //without side stepping, we're done
  isSideStep = false;
  if(p_sideStep<=0. || recursionDepth >= 3) return getNode(nearestID) + delta;

  //check whether this is a predicted collision
  bool predictedCollision=false;
  if(qr->coll_y.N){
    arr y = qr->coll_y + qr->coll_J * delta;
    if(min(y)<0.) predictedCollision = true;
  }

  if(predictedCollision && p_sideStep>0. && rnd.uni()<p_sideStep){
    isSideStep=true;

    //compute new target
    arr d = qr->getSideStep();
    d *= rnd.uni(stepsize,2.) / length(d);
    arr targ = getNode(nearestID) + d;
    bool tmp;
    return getNewSample(targ, stepsize, p_sideStep, tmp, recursionDepth + 1);
  }else{
    return getNode(nearestID) + delta;
  }

  HALT("shouldn't be here");
  return NoArr;
}

arr RRT_SingleTree::getPathFromNode(uint fromID){
  arr path;
  uint node = fromID;
  for(;;){
    path.append(ann.X[node]);
    if(!node) break;
    node = getParent(node);
  }
  path.reshape(-1, ann.X.d1);
  return path;
}

void RRT_SingleTree::glDraw(OpenGL& gl){
  glColor(.0, .0, .0);
  glLineWidth(2.f);
  glBegin(GL_LINES);
  for(uint i=1;i<getNumberNodes();i++){
    glVertex3dv(&disp3d(parent(i),0));
    glVertex3dv(&disp3d(i,0));
  }
  glEnd();
  glLineWidth(1.f);
}

//===========================================================================

bool PathFinder_RRT::growTreeTowardsRandom(RRT_SingleTree& rrt){
  const arr start = rrt.ann.X[0];
  arr t(rrt.getNode(0).N);
  t = P.sample();

  arr q = rrt.getProposalTowards(t, stepsize);

  if(cc_->isValid(q)){
    if (!checkConnection(start, q)){
      return false;
    }

    rrt.add(q, rrt.nearestID, cc_->qr);
    return true;
  }
  return false;
}

bool PathFinder_RRT::growTreeToTree(RRT_SingleTree& rrt_A, RRT_SingleTree& rrt_B, double p_forwardStep, double p_sideStep, double p_backwardStep){
  bool isSideStep, isForwardStep;
  //decide on a target: forward or random
  arr t;
  if(rnd.uni()<p_forwardStep){
    t = rrt_B.getRandomNode();
    isForwardStep = true;
  }else{
    t.resize(rrt_A.getNode(0).N);
    t = P.sample();
    isForwardStep = false;
  }

  //sample configuration towards target, possibly sideStepping
  arr q = rrt_A.getNewSample(t, stepsize, p_sideStep, isSideStep, 0);

  //evaluate the sample
  bool valid = cc_->isValid(q);

  if(isForwardStep){  n_forwardStep++; if(valid) n_forwardStepGood++; }
  if(!isForwardStep){  n_rndStep++; if(valid) n_rndStepGood++; }
  if(isSideStep){  n_sideStep++; if(valid) n_sideStepGood++; }

  //if infeasible, make a backward step from the sample configuration
  if(!valid && p_backwardStep>0. && rnd.uni()<p_backwardStep){
    t = q + cc_->qr->getBackwardStep();
    q = rrt_A.getNewSample(t, stepsize, p_sideStep, isSideStep, 0);
    
    valid = cc_->isValid(q);
    n_backStep++; if(valid) n_backStepGood++;
    if(isSideStep){  n_sideStep++; if(valid) n_sideStepGood++; }
  };

  if(valid){
    const arr start = rrt_A.ann.X[rrt_A.nearestID];
    if (!checkConnection(start, q)){
      return false;
    }

    rrt_A.add(q, rrt_A.nearestID, cc_->qr);
    double dist = rrt_B.getNearest(q);
    const arr tmp = rrt_B.ann.X[rrt_B.nearestID];
    if(dist<stepsize && checkConnection(tmp, q)) return true;
  }
  return false;
}

//===========================================================================

void PathFinder_RRT::planForward(const arr& qT){
  DISP.clear();
  DISP.copy(P.C);

  arr q = P.q0;

  RRT_SingleTree rrt(P.q0, P.query(P.q0));
  DISP.gl()->add(rrt);

  bool success=false;

  uint i;
  for(i=0;i<100000;i++){
    //let rrt0 grow
    bool added = growTreeTowardsRandom(rrt);
    if(added){
      if(length(rrt.getLast() - qT)<stepsize) success = true;
    }
    if(success) break;

    //some output
    if (verbose > 1){
      if(!(i%100)){
        DISP.setJointState(q); //updates display (makes it slow)
        DISP.watch(false);
        cout <<"RRT samples=" <<i <<" tree size = " <<rrt.getNumberNodes() <<std::endl;
      }
    }
  }

  if(!success) return;

  if (verbose > 0){
    std::cout <<"\nSUCCESS!"
              <<"\n  tested samples=" <<2*i
              <<"\n  #tree-size=" <<rrt.getNumberNodes()
     << std::endl;
  }

  arr path = rrt.getPathFromNode(rrt.nearestID);
  revertPath(path);

  //display
  if(verbose > 1){
    std::cout << "path-length= " << path.d0 <<endl;
    DISP.proxies.clear();

    for(uint t=0;t<path.d0;t++){
      DISP.setJointState(path[t]);
      //DISP.watch();
      DISP.watch(false);
      rai::wait(.1);
    }
  }

  path >>FILE("z.path");
}

// TODO: This can probably be made more efficient for our use-case
// Produce the whole sequence in one go
double corput(uint n, uint base){
  double q=0;
  double bk=(double)1/base;

  while (n > 0) {
    q += (n % base)*bk;
    n /= base;
    bk /= base;
  }
  return q;
}

bool PathFinder_RRT::checkConnection(
    const arr &start,
    const arr &end,
    const double res){
  uint disc = uint(length(start-end) / res);
  if (disc < 2){
    disc = 2;
  }
  for (uint i=1; i<disc; ++i){
    double ind = corput(i, 2);
    const arr p = start + ind * (end-start);

    // TODO: change to check feasibility properly (with path constraints)
    if(!cc_->isValid(p)){
      return false;
    }

    //C.watch(false);
  }
  return true;
}

arr PathFinder_RRT::planConnect(const arr& qT, double p_forwardStep, double p_sideStep, double p_backwardStep){

  RRT_SingleTree rrt0(P.q0, P.query(P.q0));
  RRT_SingleTree rrtT(qT, P.query(qT));

  if(verbose>1){
    DISP.clear();
    DISP.copy(P.C);
    DISP.gl()->add(rrt0);
    DISP.gl()->add(rrtT);
  }

  bool success=false;

  uint i;
  for(i=0;i<maxIter;i++){
    success = growTreeToTree(rrt0, rrtT, p_forwardStep, p_sideStep, p_backwardStep);
    if(success) break;

    success = growTreeToTree(rrtT, rrt0, p_forwardStep, p_sideStep, p_backwardStep);
    if(success) break;

    //some output
    if (verbose > 1){
      if(!(i%100)){
        DISP.setJointState(rrt0.getLast());
        DISP.watch(false);
        std::cout <<"iteration " << i << " RRT queries=" <<P.evals <<" tree sizes = " <<rrt0.getNumberNodes()  <<' ' <<rrtT.getNumberNodes() <<std::endl;
      }
    }
  }

  if(!success) return NoArr;

  if (verbose > 0){
    cout <<"\nSUCCESS!" <<endl;
    cout <<"  RRT queries=" <<P.evals <<" tree sizes = " <<rrt0.getNumberNodes()  <<' ' <<rrtT.getNumberNodes() <<std::endl;
    cout <<"  forwardSteps: " <<(100.*n_forwardStepGood/n_forwardStep) <<"%/" <<n_forwardStep;
    cout <<"  backSteps: " <<(100.*n_backStepGood/n_backStep) <<"%/" <<n_backStep;
    cout <<"  rndSteps: " <<(100.*n_rndStepGood/n_rndStep) <<"%/" <<n_rndStep;
    cout <<"  sideSteps: " <<(100.*n_sideStepGood/n_sideStep) <<"%/" <<n_sideStep;
    cout <<endl;
  }

  arr path = rrt0.getPathFromNode(rrt0.nearestID);
  arr pathT = rrtT.getPathFromNode(rrtT.nearestID);

  revertPath(path);
  path.append(pathT);

  //display
  if (verbose > 1){
    cout <<"  path-length=" <<path.d0 <<endl;
    DISP.proxies.clear();
    DISP.watch(true);
    rai::wait(.2);
    for(uint t=0;t<path.d0;t++){
      DISP.setJointState(path[t]);
      DISP.watch();
      rai::wait(.2);
    }

    DISP.clear();
  }

  path >>FILE("z.path");
  return path;
}

void revertPath(arr& path){
  uint N = path.d0;
  arr x;
  for(uint i=0;i<N/2;i++){
    x = path[i];
    path[i] = path[N-1-i];
    path[N-1-i] = x;
  }
}

ptr<PathResult> PathFinder_RRT::run(double timeBudget){
  CHECK(goals.N, "no goal given")
  arr path = planConnect(goals, .5, .9, .9);

  return make_shared<PathResult>(path);
}

