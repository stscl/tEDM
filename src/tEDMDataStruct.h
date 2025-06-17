#ifndef tEDMDataStruct_H
#define tEDMDataStruct_H

#include <vector>

struct PartialCorRes {
  int first;
  double second;
  double third;

  // Default constructor
  PartialCorRes() : first(0), second(0.0), third(0.0) {}

  // Constructor to initialize all members
  PartialCorRes(int f, double s, double t) : first(f), second(s), third(t) {}
};

struct DeLongPlacementsRes {
  double theta;
  std::vector<double> X;
  std::vector<double> Y;

  // Default constructor
  DeLongPlacementsRes() : theta(0.0), X(), Y() {}

  // Parameterized constructor
  DeLongPlacementsRes(double t, const std::vector<double>& x, const std::vector<double>& y)
    : theta(t), X(x), Y(y) {}
};

#endif // tEDMDataStruct_H
