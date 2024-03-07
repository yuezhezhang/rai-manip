#include <KOMO/komo.h>
#include <Kin/forceExchange.h>
#include <Kin/F_forces.h>

#include <Assembly/Dependency.h>

// this custom comparator is necessary to avoid
// rai::String comparisons, which are ill-defined
bool rai_str_cmp(const std::pair<double, rai::String> &a,
              const std::pair<double, rai::String> &b){
    return (a.first < b.first);
}

// Derive the base-classes for all the following cost-classes
// - enables usage of 'ground parts', neighbors, and already placed parts
struct AssemblyCost{
  rai::Configuration C;

  NeighborHandler nh;
  AssemblyCost(const NeighborHandler &_nh) : nh(_nh) {};

  // compute distance between the parts:
  // can be used as a secondary decision on which pieced to place
  // works by checking where the goal of the piece is
  double dists(const StringA &potentialParts){
    if(potentialParts.d0 == 0) {return 0;}

    double d = 0.;
    for (const rai::String p: potentialParts){
      rai::Frame *f = C[p];
      rai::Frame *r = C[STRING("m0_mee")];
      if (r){
        d += length(r->ensure_X().pos.getArr() - f->ensure_X().pos.getArr());
      }

      // sum up distance between all potential parts
      for (const rai::String q: potentialParts){
        rai::Frame *g = C[q];
        d += length(g->ensure_X().pos.getArr() - f->ensure_X().pos.getArr());
      }
    }
    return d;
  }

  // functor that is called to be evaluated on each action sequence to rank them
  double operator()(const StringA &placedParts, const StringA& potentialParts){
    // set potential parts to final position
    // TODO: add back in again
    /*
    for (const auto &o: potentialParts){
      rai::Frame *f = C[o];
      rai::Frame *g = C[STRING("o"<<o<<"g")];
      f->set_X() = g->ensure_X();
    }*/

    return this->operator()(placedParts, potentialParts) + dists(potentialParts)/100;
  }
};


// Stability evaluation:
// Sets up a system of equations describing the forces
// using komo (and solves it)
struct MaxTorqueCost: public AssemblyCost{
  // used to save the previous computational result as attempt to speed up the computation
  std::map<int, arr> prevForce;
  std::map<int, arr> prevPos;
  
  MaxTorqueCost(const NeighborHandler &_nh)
    : AssemblyCost(_nh){}

  uint _computeKey(uint a, uint b){
    if (a < b){
      return 1000 * a + b;
    }
    return 1000 * b + a;
  }

  double solve(uintA& freeObjects, uintA& pairs){
    while(C.forces.N) {delete C.forces.last();}

    //create force variables for pairs
    for(uint i=0;i<pairs.d0;i++){
      rai::Frame *a = C.frames(pairs(i,0));
      rai::Frame *b = C.frames(pairs(i,1));

      //only add pairs that interact with free objects
      if(!freeObjects.contains(a->ID) && !freeObjects.contains(b->ID)) {
        continue; 
      }

      rai::ForceExchange *c = new rai::ForceExchange(*a, *b, rai::FXT_poa);
      if (prevForce.count(_computeKey(a->ID, b->ID)) > 0){
        c->force = prevForce[_computeKey(a->ID, b->ID)];
        c->poa = prevPos[_computeKey(a->ID, b->ID)];
      }
      else{
        c->poa = (.5*(a->ensure_X().pos + b->ensure_X().pos)).getArr();
        c->force = 0.0;
      }
      C._state_q_isGood = false;
    }

    //create a komo instance
    KOMO komo;
    komo.setModel(C, false);
    komo.setIKOpt();

    for(rai::ForceExchange* con:C.forces){
      //    komo.addObjective({1., 1.}, make_shared<TM_Contact_POAisInIntersection_InEq>(from, to, .01), OT_ineq, {1e1});
      komo.addObjective({1., 1.}, make_shared<F_fex_Force>(), {con->a.name, con->b.name}, OT_sos, {1e-2});
      //    komo.addObjective({1., 1.}, make_shared<TM_Contact_ForceIsNormal>(from, to), OT_sos, {1e-1});
    }

    for(uint id:freeObjects){
      komo.addObjective({1., 1.}, make_shared<F_TotalForce>(false), {C.frames(id)->name}, OT_sos, {1e1});
    }

    komo.verbose = 0;
    // komo.animateOptimization= 2;
    komo.reset();
    OptOptions opt;
    opt.stopTolerance = 1e-1;
    opt.stopIters = 10;
    opt.stopEvals = 100;
    opt.damping   = 1;
    opt.maxStep = 10;
    komo.run(opt);

    // extract results
    double maxTorque = 0;
    double sumTorque = 0;
    arr avgTorqueArr = ARR(0, 0, 0);
    doubleA torques = consts<double>(0, C.frames.N);
    intA cnts = consts<int>(0, C.frames.N);

    for(auto* c: komo.pathConfig.forces) {
      const arr pt = (.5*(c->a.ensure_X().pos + c->b.ensure_X().pos)).getArr();
      const arr dist =  c->poa - pt;

      const arr force =  c->force;
      const arr moment = crossProduct(dist, force);
      const double torque = length(moment);

      torques(C[c->a.name]->ID) += torque;
      torques(C[c->b.name]->ID) += torque;

      cnts(C[c->a.name]->ID) += 1;
      cnts(C[c->b.name]->ID) += 1;

      sumTorque += torque;
      avgTorqueArr += moment;

      if (torque > maxTorque){
        maxTorque = torque;
      }

      prevForce[_computeKey(c->a.ID, c->b.ID)] = force;
      prevPos[_computeKey(c->a.ID, c->b.ID)] = c->poa;
    }

    const double avgTorque = sumTorque / pairs.d0;
    avgTorqueArr /= 1.0 * pairs.d0;
    std::cout << avgTorque << " " << maxTorque << " " << avgTorqueArr << std::endl;

    // return length(avgTorqueArr);
    // return sumTorque;
    return avgTorque;
    // return maxTorque;
  }

  // actual functor that is called
  double operator()(const StringA &placedParts, const StringA &potentialParts){
    StringA parts = placedParts;
    parts.append(potentialParts);

    // what are neighboring pieces of the to be placed parts?
    uintA pairs(uint(0), uint(2));
    for(const rai::String p: parts){
      for(const rai::String n: placedParts){
        if (nh.getNeighbours(p).contains(n)){
          auto f = C[p];
          auto r = C[n];

          if (!pairs.contains({f->ID, r->ID}) && !pairs.contains({r->ID, f->ID})){
            pairs.append({f->ID, r->ID});
          }
        }
      }
    }

    uintA freeObjects;
    for(const rai::String p: parts){
      auto f = C[p];
      if (!f) continue;
      if(f->shape){
        f->setJoint(rai::JT_rigid);
        f->setMass(.1);

        if(!nh.groundParts.contains(p)){
          freeObjects.append(f->ID);
        }
      }
    }

    double res = 0;
    if (freeObjects.N == 0){
      res = -100;
    }
    else{
      res = solve(freeObjects, pairs);
    }

    return res + potentialParts.d0*100;
  }
};

struct RandomCost: public AssemblyCost{
  RandomCost(const NeighborHandler &_nh) : AssemblyCost(_nh) {};

  double operator()(const StringA &placedParts, const StringA &potentialParts){
    return std::rand()/(RAND_MAX + 1.);
  };
  
};

struct ZCost: public AssemblyCost{
  ZCost(const NeighborHandler &_nh) : AssemblyCost(_nh) {};

  double operator()(const StringA &placedParts, const StringA &potentialParts){
    return nh.C[potentialParts(0)]->getPosition()(2);
  };
};

struct XYCost: public AssemblyCost{
  XYCost(const NeighborHandler &_nh) : AssemblyCost(_nh) {};

  double operator()(const StringA &placedParts, const StringA &potentialParts){
    //std::cout << nh.C[potentialParts(0)]->getPosition()({0,1}) << std::endl;
    //std::cout <<  std::fabs(nh.C[potentialParts(0)]->getPosition()(0)) << std::endl;
    return std::fabs(nh.C[potentialParts(0)]->getPosition()(0));
    //return length(nh.C[potentialParts(0)]->getPosition()({0,1}));
  };
};

struct MaxNeighborsCost: public AssemblyCost{
  MaxNeighborsCost(const NeighborHandler &_nh) : AssemblyCost(_nh) {};

  double operator()(const StringA &placedParts, const StringA &potentialParts){
    double cost = 0;
    doubleA pc(potentialParts.d0);

    uint cnt = 0;
    for(auto i: potentialParts){
      if(nh.groundParts.contains(i)) { // treat ground as 'neighbor'
        pc(cnt) += 1;
        cost += 3;
      }
      for (auto j: placedParts){
        if(nh.getNeighbours(i).contains(j)){ // count actual placed neighbors 
          pc(cnt) += 1;
          ++cost;
        }
      }

      ++cnt;
    }

    const double minCost = pc.min();
    const double avgCost = cost / potentialParts.d0;
    
    // return -minCost + potentialParts.d0*100;
    return -avgCost + potentialParts.d0*100;
  }
};

struct DisplacementCost: public AssemblyCost{};
