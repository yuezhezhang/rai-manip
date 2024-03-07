#include <Kin/kin.h>
#include <Core/array.h>

#include <Assembly/costs.h>

arr gnuplot(const double x){
  double r = std::sqrt(x);
  double g = x * x * x;
  double b = std::sin(x * 2 * RAI_PI);

  return ARR(r, g, b);
}

arr gnuplot2(const double x){
  double r = x / 0.32 - 0.78125;

  double g = 2 * x - 0.84;

  double b;
  if (x < 0.25){ b = 4 * x;}
  else if ((x >= 0.25) && (x < 0.92)){b = -2 * x + 1.84;}
  else{ b = x / 0.08 - 11.5;}
    
  arr c = ARR(r, g, b);
  return c;
}

arr matplotlib_cycle(const uint i, const bool render=true){
  // matplotlib
  arr c;
  
  if (!render){
    c = {0.12156862745098039, 0.4666666666666667, 0.7058823529411765,
      1.0, 0.4980392156862745, 0.054901960784313725,
      0.17254901960784313, 0.6274509803921569, 0.17254901960784313,
      0.8392156862745098, 0.15294117647058825, 0.1568627450980392,
      0.5803921568627451, 0.403921568627451, 0.7411764705882353,
      0.5490196078431373, 0.33725490196078434, 0.29411764705882354,
      0.8901960784313725, 0.4666666666666667, 0.7607843137254902,
      0.4980392156862745, 0.4980392156862745, 0.4980392156862745,
      0.7372549019607844, 0.7411764705882353, 0.13333333333333333,
      0.09019607843137255, 0.7450980392156863, 0.8117647058823529};
  }
  
  // for rendering
  else {
    c = {1., 0., 0.,
             0., 0., 1.,
             0.17254901960784313, 0.6274509803921569, 0.17254901960784313,
             0, 1, 0
    };
  }

  c.reshape(-1, 3);

  return c[i%c.d0]();
}

//https://en.wikipedia.org/wiki/Halton_sequence
double halton(int n) {
  int base = 10;
  double f = 1;
  double r = 0;
  while(n > 0)
  {
    f = f/(double)base;
    r = r + f*(n%base);
    n = floor(n/(double)base);
  }
  return r;
}

arr halton_color(const int n) {
  arr c;
  c.append(halton(3*n));
  c.append(halton(3*n+1));
  c.append(halton(3*n+2));
  return c;
}

arr changeSaturation(const arr &col, double change) {
  const double Pr = .299;
  const double Pg = .587;
  const double Pb = .114;

  double R = col(0);
  double G = col(1);
  double B = col(2);

  const double  P=sqrt(
  R*R*Pr+
  G*G*Pg+
  B*B*Pb ) ;

  R=P+(R-P)*change;
  G=P+(G-P)*change;
  B=P+(B-P)*change;

  return {R, G, B};
}

struct Assembly{
  struct AssemblyPart{
    enum State{unplaced=1, placed, transition};

    rai::String name;

    arr initialPose;
    arr finalPose;

    State state = unplaced;
    int tPlaced = -1;
  };

  Assembly(const rai::String &goalStateFile, const rai::String &initialStateFile="", const bool rotate=false){
    rai::Configuration Cgoal;
    Cgoal.addFile(goalStateFile);

    uint cnt = 0;
    for (auto f: Cgoal.frames){
      const double scale = rai::getParameter<double>("assembly/scaleModel", 1.);  // display stuff
      if (f->shape){
        if (!f->name.startsWith("_")) {f->name = STRING("b"<<cnt++);}

        f->setContact(1);

        // Scale the full model
        f->getShape().mesh().scale(scale);
        f->setPosition(f->getPosition()*scale);
      }
    }

    // check feasibility of final configuration
    Cgoal.ensure_q();
    //Cgoal.watch(true);

    Cgoal.stepSwift();
    const double d = Cgoal.getTotalPenetration();
    std::cout << "Total penetration: " << d << std::endl;

    if (d > .1){
      HALT("Model not buildable");
    }


    if (goalStateFile.contains("fit")){
      makeConvexHulls(Cgoal.frames);
    }

    rai::Configuration Cinit;
    if (initialStateFile.N == 0){
      // distribute the pieces from the final state in a circle
      uint cnt = 0;
      for (auto *f: Cgoal.frames){
        if (!f->shape || f->name.startsWith("_")) {continue;}

        // place piecse in a spiral
        const double r = cnt * 0.05 + 3;
        //const double r = cnt * 0.02 + 3;
        arr p;
        if (goalStateFile.contains("fit")){
          p = {r*std::sin(cnt/2.), r*std::cos(cnt/2.), 0.15};
        }
        else{
          p = {r*std::sin(cnt/2.), r*std::cos(cnt/2.), 0.1};
        }
        
        rai::Frame *n = new rai::Frame(Cinit, f);
        n->setPosition(p);

        // rotate pieces to be flat: approximate normal of the piece to be the one from
        // center of the piece to the absolute center
        if (goalStateFile.contains("fit")){
          auto mesh = f->getShape().mesh();
          mesh.fuseNearVertices(1e-5);

          arr W;
          arr v;
          arr Y;

          pca(Y, v, W, mesh.V);
          arr mat = W;
          transpose(mat);

          double w = std::sqrt(1.0 + mat(0,0) + mat(1,1) + mat(2, 2)) / 2.0;
          double w4 = (4.0 * w);
          double x = (mat(2, 1) - mat(1,2)) / w4 ;
          double y = (mat(0,2) - mat(2, 0)) / w4 ;
          double z = (mat(1,0) - mat(0,1)) / w4 ;

          n->setQuaternion({w, x, y, z});
        }
        ++cnt;
      }
    }
    else{
      Cinit.addFile(initialStateFile);
    }

    *this = Assembly(Cgoal, Cinit, false);
  };

  Assembly(rai::Configuration &Cgoal, rai::Configuration &Cinit, bool failIfNotSame=true){
    n = new NeighborHandler(Cgoal);
    for (auto *f: Cgoal.frames){
      if (!f->shape || f->name.startsWith("_")){continue;}
      if (!Cinit[f->name]){
        if (failIfNotSame) {
          HALT("Frame in goal but not in initial configuration");
        }
        else{
          std::cout << "Frame " << f->name << " in goal but not in initial configuration" << std::endl;;
          continue;
        }
      }
      
      AssemblyPart p;
      p.name = f->name;
      p.finalPose = f->getPose();
      p.initialPose = Cinit[f->name]->getPose();

      parts.insert({p.name, p});
    }
  }

  struct cmp {
    bool operator()(const rai::String& a, const rai::String& b) const {
      // TODO: implement less-than comparison in rai::String
      std::string c(a.p);
      std::string d(b.p);
      return c < d;
    }
  };

  std::map<rai::String, AssemblyPart, cmp> parts;

  NeighborHandler *n;

  StringA getPlacedParts(){
    rai::Array<rai::String> placed;

    for (auto e: parts){
      if (e.second.state == AssemblyPart::State::placed){
        placed.append(e.first);
      }
    }

    return placed;
  };

  // TODO: properly deal with 'agent'-how to rank
  StringA rank(const StringA &candidates, const arr& robotPos){
    const bool zRanking = rai::getParameter<bool>("assembly/rankZ", false);
    const bool xyRanking = rai::getParameter<bool>("assembly/rankXY", false);
    const bool rndRanking = rai::getParameter<bool>("assembly/rankRnd", false);

    const StringA placed = getPlacedParts();

    MaxNeighborsCost c(*n);
    RandomCost rnd(*n);
    ZCost z(*n);
    XYCost xy(*n);

    // permute tmp randomly to ensure that equal values are not ordered according to the ordering in candidates
    StringA candidates_shuffled = candidates;
    candidates_shuffled.permuteRandomly();

    std::vector<std::pair<double, rai::String>> tmp;
      //std::cout << "G" << std::endl;
    for (const rai::String obj: candidates_shuffled){
      //std::cout << "B" << std::endl;
      const double dist = length(robotPos({0, 2}) - parts[obj].finalPose({0,2}));
      
      double cost = 0;
      cost += c(placed, {obj});
      if (zRanking) {cost += z(placed, {obj}); }
      if (rndRanking) {cost += 1.5 * rnd(placed, {obj}); }
      if (xyRanking) {cost += 2. * xy(placed, {obj}); }
      
      tmp.emplace_back(cost, obj);
    }

    //std::cout << "G" << std::endl;
    // sort according to cost
    std::sort(tmp.begin(), tmp.end(), rai_str_cmp);

    StringA sorted;
    for(auto elem: tmp){
      sorted.append(elem.second);
    }

    return sorted;
  }

  StringA getFeasibleParts(){
    const StringA placed = getPlacedParts();

    std::map<rai::String, uint, cmp> numNeighborsPlaced;

    // iterate over all placed parts and add the neighbors
    for (const rai::String &str: placed){
      for (rai::String tmp: n->getNeighbours(str)){
        numNeighborsPlaced[tmp] += 1;
      }
    }

    // iterate over the ground parts and add them
    for (rai::String p: n->groundParts){
      numNeighborsPlaced[p] += 1;
    }

    StringA feasible;
    for (auto e: numNeighborsPlaced){
      if (e.second >= 1 && !placed.contains(e.first)){
        feasible.append(e.first);        
      }
    }

    return feasible;
  };

  bool removeParts(const rai::Array<rai::String> &partNames);
  bool removeParts(const uint num){
    const uint s = parts.size();
    if (num >= s){
      HALT("Attempting to remove too many elements");
    }

    auto begin = parts.begin();
    std::advance(begin, num+1);
    auto end = parts.end();

    parts.erase(begin, end);
  }

  /*bool setNeighbors(const rai::String &neighborsFile){
    n.neighbors << FILE(neighborsFile);
  };*/

  void setConfigurationToStart(rai::Configuration &C, bool removeIfNotInList=true){
    //std::cout << parts.size() << std::endl;

    for (auto *f: C.frames){
      if (!f->shape){continue;}
      if (parts.find(f->name) == parts.end()) {
        if (removeIfNotInList){
          std::cout << "Removing " << f->name << std::endl;
          delete f;
        }
        std::cout << "Skipping " << f->name << std::endl;
        continue;
      }

      const arr pose = parts[f->name].initialPose;
      f->setPose(pose);
    }

    //C.watch(true);
    
    C.swiftDelete();
    C.stepSwift();
  };

  void setConfigurationToPiece(rai::Configuration &C, const uint startPart=0, const bool removeIfNotInList=true){
    setConfigurationToStart(C, removeIfNotInList);

    for (uint i=0; i<startPart; ++i){
      // figure out which parts are urrently feasible
      const StringA feasible = getFeasibleParts();

      // rank them
      const StringA ranked = rank(feasible, {0., 0., 0.}); 
      
      // place the first one
      const rai::String obj = ranked(0);
      auto &part = parts[obj];
      part.state = AssemblyPart::State::placed;

      C[obj]->setColor(gnuplot2(1.*i/parts.size()));
      C[obj]->setContact(1);
      C[obj]->setPose(part.finalPose);
    }
  }

};
