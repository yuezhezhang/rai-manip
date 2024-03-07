#pragma once

#include <Core/array.h>

struct TimedPath{
  arr path;
  arr time;

  TimedPath(const arr &_path, const arr &_time): path(_path), time(_time){};

  arr resample(const arr &times, rai::Configuration &C){
    auto periodicDimensions = std::vector<bool>(C.getJointState().N, false);

    for (auto *j: C.activeJoints){
      if (j->type == rai::JT_hingeX || j->type == rai::JT_hingeY|| j->type == rai::JT_hingeZ){
        periodicDimensions[j->qIndex] = true;
      }
    }
    
    arr p(0u, path.d1);
    for(double t: times){
      arr pos = getPos(t, C, periodicDimensions);
      p.append(pos);
    }

    return p;
  };

  arr getPos(const double t, rai::Configuration &C, const std::vector<bool> &periodicDimensions){
    if (t <= time(0)){
      return path[0];
    }
    if (t >= time(time.d0-1)){
      return path[time.d0-1];
    }

    for(uint i=0; i<time.N-1; ++i){
      if (time(i) <= t && t < time(i+1)){
        const arr p1 = path[i]();
        const arr p2 = path[i+1]();
        arr delta = p2 - p1;
        for (uint l=0; l<delta.N; ++l){
          if (periodicDimensions[l]){
            // this is an angular joint -> we need to check the other direction
            const double start = p2(l);
            const double end = p1(l);
            delta(l) = std::fmod(start - end + 3.*RAI_PI, 2*RAI_PI) - RAI_PI;
          }
        }

        return path[i] + (t - time(i))/(time(i+1) - time(i))*(delta);
        //return path[i];
      }
    }

    // To ensure there is a point returned
    std::cout << "Something went wrong!" << std::endl;
    return path[time.d0-1];

  };
};

