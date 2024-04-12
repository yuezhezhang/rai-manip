#include <unordered_map>
#include "Animation.h"

#include <Kin/frame.h>

void rai::Animation::AnimationPart::write(std::ostream& os) const{
  os <<'\n';
  arr tmp(1);
  tmp = start;
  tmp.writeTagged(os, "start");
  os <<'\n';
  frameIDs.writeTagged(os, "frameIDs");
  os <<'\n';
  frameNames.writeTagged(os, "frameNames");
  os <<'\n';
  frameCols.writeTagged(os, "frameColors");
  os <<'\n';
  X.writeTagged(os, "poses");
}

void rai::Animation::AnimationPart::read(std::istream& is) {
  arr tmp;
  tmp.readTagged(is, "start");
  start = tmp(0);
  frameIDs.readTagged(is, "frameIDs");
  frameNames.readTagged(is, "frameNames");
  frameCols.readTagged(is, "frameColors");
  X.readTagged(is, "poses");
}

void rai::Animation::read(istream& is) {
  A.readTagged(is, "animation");
}

uint rai::Animation::getT() const{
  uint T=0;
  for(auto& a: A) {
    if(a.X.d0 + a.start>T) {
      T=a.X.d0 + a.start;
    }
  }

  return T;
}

void rai::Animation::setToTime(rai::Configuration& C, const double t, const double tIgnoreBefore) const{
  // TODO: make sure that the ordering is used to speed this up
  // - we know for a lot of things that they have not started yet
  if (A.d0 == 0){
    return;
  }

  // arr poses(0, 7);
  // uintA ids;

  // for (const auto *f: C.frames){
  //   const uint id = f->ID;
  //   arr tmp;
  //   uint tDiff;

  //   for(int i=A.d0-1; i>=0; --i){
  //     // this assumes that the newest animation part always has priority
  //     const auto &a = A(i);
  //     const double &start = a.start;

  //     if (t < start || start < 0 || a.X.d0 == 0) {
  //       //LOG(-1) <<"Warning: start < 0";
  //       continue;
  //     }
  //     /*if (std::ceil(t) > start + a.X.d0 + 30){
  //       //LOG(-1) <<  "Warning: time > anim_part";
  //       continue;
  //     }*/

  //     // enables us to reduce the things we have to look at if we are certain that
  //     // animation parts before a certain time can be disregarded
  //     if (tIgnoreBefore > start + a.X.d0){
  //       continue;
  //     }

  //     const uint dt = std::max({uint(t - (start + a.X.d0-1)), 0u});

  //     if (tmp.N == 0 || dt < tDiff){
  //       const int idx = a.frameIDs.findValue(id);
  //       if (idx == -1){
  //         continue;
  //       }
        
  //       if (t >= start + a.X.d0-1){
  //         // tmp = a.X[a.X.d0-1];
  //         tmp.referToDim(a.X, a.X.d0-1, idx);
  //         tDiff = t - (start + a.X.d0-1);
  //       }
  //       else{
  //         // we might want to interpolate here - quaternion interpolation could be a problem
  //         // tmp = a.X[uint(std::floor(t - start))];
  //         tmp.referToDim(a.X, uint(std::floor(t - start)), idx);
  //         tDiff = 0;
  //       }
  //     }
  //   }

  //   if (tmp.N > 0){
  //     poses.append(tmp);
  //     ids.append(id);
  //   }
  // }

  // auto start_time = std::chrono::high_resolution_clock::now();

  // datastructure to collect the positions
  std::unordered_map<uint, uint> timeMap;

  arr tmp;

  for(int i=A.d0-1; i>=0; --i){
    // this assumes that the newest animation part always has priority
    const auto &a = A(i);
    const double &start = a.start;

    if (t < start || start < 0 || a.X.d0 == 0) {
      //LOG(-1) <<"Warning: start < 0";
      continue;
    }
    /*if (std::ceil(t) > start + a.X.d0 + 30){
      //LOG(-1) <<  "Warning: time > anim_part";
      continue;
    }*/

    // enables us to reduce the things we have to look at if we are certain that
    // animation parts before a certain time can be disregarded
    if (tIgnoreBefore > start + a.X.d0){
      continue;
    }

    uint tDiff;
    if (t >= start + a.X.d0-1){
      // tmp = a.X[a.X.d0-1];
      tmp.referToDim(a.X, a.X.d0-1);
      tDiff = t - (start + a.X.d0-1);
    }
    else{
      // we might want to interpolate here - quaternion interpolation could be a problem
      // tmp = a.X[uint(std::floor(t - start))];
      tmp.referToDim(a.X, uint(std::floor(t - start)));
      tDiff = 0;
    }

    uint cnt = 0;
    for (const auto &id: a.frameIDs){
      if (timeMap.count(id) == 0 || timeMap[id] > tDiff){
        poseMap[id].referToDim(tmp, cnt);
        timeMap[id] = tDiff;
      }
      cnt++;
    }
  }

  
  arr poses(timeMap.size(), 7);
  uintA set_ids(timeMap.size());

  uint cnt = 0;
  for (const auto &e: timeMap){
    //if (length(C.getFrameState(uintA({e.first})) - e.second) < 1e-3) {continue;}

    set_ids(cnt) = e.first;
    poses[cnt]() = poseMap[e.first];

    // ids.append(e.first);
    // poses.append(e.second);

    ++cnt;
  }

  C.setFrameState(poses, set_ids);

  // // go over all the objects that have not been set now
  // arr unset_ids;
  // for (const auto f: C.frames){
  //   if (set_ids.contains(f->ID)){continue;}
  //   unset_ids.append(f->ID);

  //   // std::cout << f->name << std::endl;
  // }

  // timeMap.clear();
  
  // for (uint i=0; i<A.d0; ++i){
  //   const auto &a = A(i);
  //   const double &start = a.start;

  //   if (start < t){
  //     continue;
  //   }

  //   for (const uint id: unset_ids){
  //     if (a.frameIDs.contains(id)){
  //       uint ind = a.frameIDs.findValue(id);
  //       if (timeMap.count(id) == 0 || start < timeMap[id]){
  //         timeMap[id] = start;

  //         tmp.referToDim(a.X, 0);
  //         poseMap[id].referToDim(tmp, ind);
  //       }
  //     }
  //   }
  // }

  // {
  //   arr poses(timeMap.size(), 7);
  //   uintA ids(timeMap.size());

  //   uint cnt = 0;
  //   for (const auto &e : timeMap) {
  //     // if (length(C.getFrameState(uintA({e.first})) - e.second) < 1e-3)
  //     // {continue;}

  //     ids(cnt) = e.first;
  //     poses[cnt]() = poseMap[e.first];

  //     // ids.append(e.first);
  //     // poses.append(e.second);

  //     ++cnt;
  //   }

  //   C.setFrameState(poses, ids);
  // }

  // auto end_time = std::chrono::high_resolution_clock::now();
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end_time - start_time ).count();

  /*{
    std::ofstream f;
    f.open("./out/Anim_time.txt", std::ios_base::app);
    f << duration << std::endl;
  }*/
}

void rai::Animation::play(rai::Configuration& C, bool pause, bool exp){
  rai::ConfigurationViewer V;
  V.setConfiguration(C);
  V.watch(STRING("\"Multi-Robot Assembly.\" Here we go."));
  uint T=getT();

  for(uint t=0;t<T;t++){
    setToTime(C, t);
    V.setConfiguration(C, STRING("Animation t:" <<t), pause);

    if (exp){
      // rotate camera
      if (true){
        const double r = std::sqrt(V.displayCamera().X.pos.x*V.displayCamera().X.pos.x + V.displayCamera().X.pos.y*V.displayCamera().X.pos.y);
        V.displayCamera().X.pos = rai::Vector(r*std::sin(t*.0025), r*std::cos(t*.0025), V.displayCamera().X.pos.z);
        V.displayCamera().focusOrigin();
        V.displayCamera().upright();
      }

      // std::string fname = "z.vid/world"+std::to_string(t)+".dae";
      // C.writeCollada(fname.c_str());
      //rai::wait(dispWait);
      V.savePng();
    }

    rai::wait(.1);
  }

}

void rai::Animation::drawTrajectories(rai::Configuration& C){
  for(int i=0; i<A.d0; ++i){
    // this assumes that the newest animation part always has priority
    auto &a = A(i);
    DrawPaths *dp = new DrawPaths(a.X);
    C.gl()->add(*dp);
  }
  C.watch(true);
}

void rai::Animation::write(ostream& os) const {
  A.writeTagged(os, "animation");
}

