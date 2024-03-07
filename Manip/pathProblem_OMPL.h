#pragma once

#include <ompl/base/StateValidityChecker.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/spaces/TimeStateSpace.h>

#include <ompl/base/State.h>
#include <ompl/base/ScopedState.h>

// for the interface between rai and ompl
#include <PlanningSubroutines/ConfigurationProblem.h>

namespace ob = ompl::base;

struct ConfigurationProblem_OMPL : ob::StateValidityChecker {
  ConfigurationProblem& P;
  //arr x;
  mutable uint numQueries = 0;

  ConfigurationProblem_OMPL(ConfigurationProblem& _P, const ob::SpaceInformationPtr& si)
    : ob::StateValidityChecker(si),
      P(_P) {}

  virtual ~ConfigurationProblem_OMPL() = default;

  virtual bool isValid(const ompl::base::State *state) const{
    const double* _state = state->as<ob::RealVectorStateSpace::StateType>()->values;

    this->numQueries++;
    if(!(numQueries%100)) cout <<"#queries: " <<numQueries <<endl;

    arr x(_state, P.q0.N, true);
    auto qr = P.query(x);

    return qr->isFeasible;
  }
};

struct TimedConfigurationProblem_OMPL : ob::StateValidityChecker {
  //ConfigurationProblem& P;
  TimedConfigurationProblem &TP;
  //arr q;
  mutable uint numQueries = 0;

  TimedConfigurationProblem_OMPL(TimedConfigurationProblem& _TP, const ob::SpaceInformationPtr& si)
    : ob::StateValidityChecker(si),
      TP(_TP) {}

  virtual ~TimedConfigurationProblem_OMPL() = default;

  virtual bool isValid(const ompl::base::State *state) const{
    //const double _state = state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(0)->position;
    const double* _state = state->as<ob::CompoundState>()->as<ob::RealVectorStateSpace::StateType>(0)->values;
    const double time = state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;
    //const double* _state = state->as<ob::RealVectorStateSpace::StateType>()->values;
  
    this->numQueries++;
    //if(!(numQueries%100)) cout <<"#queries: " <<numQueries <<endl;

    arr q(_state, TP.q0.N, true);
    //auto qr = P.query(x);
    // Check if start q is feasible
    try{
      auto qr = TP.query(q, time);
      //std::cout << qr->isFeasible << " " << q << " " << time << std::endl;
      //TP.A.setToTime(TP.C, time);
      //std::cout << time << std::endl;
      if (true || time != 0) {
        //std::cout << time << std::endl; 
        //if (qr->isFeasible) TP.C.watch(false);
        /*else if (!qr->isFeasible) {
          qr->writeDetails(cout, TP.C);
          TP.C.watch(true);
        }*/
      }

      //if (qr->isFeasible) TP.C.watch(true);

      //std::cout << qr->isFeasible << std::endl;

      return qr->isFeasible;
    }
    catch(...){
      std::cout << "lul" << std::endl;
      return false;
    }
  }
};
