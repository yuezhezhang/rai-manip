#pragma once

#include "ConfigurationProblem.h"

struct ConfigurationSampler {
  ConfigurationSampler(ConfigurationProblem& p) : p(p) {}

  struct Result{
    bool feasible;
    arr goal;
    ptr<QueryResult> qr;
  };
  ptr<Result> run();

private:
  ConfigurationProblem& p;
};

struct TimedConfigurationSampler: ConfigurationSampler {
  TimedConfigurationSampler(TimedConfigurationProblem& _tp) : ConfigurationSampler(_tp), tp(_tp) {}
  ptr<Result> run(const double t);

private:
  TimedConfigurationProblem& tp;
};
