#include "KLFitter/ResDoubleCB_3.h"

KLFitter::ResDoubleCB_3::ResDoubleCB_3(const char * filename) :
  KLFitter::ResDoubleCBBase(filename, ((4 * 4) + (1 * 2))) {
}

KLFitter::ResDoubleCB_3::~ResDoubleCB_3() = default;

double KLFitter::ResDoubleCB_3::GetMean(double x) {
  return evaluate(x, 0);
}

double KLFitter::ResDoubleCB_3::GetSigma(double x) {
  return evaluate(x, 4);
}

double KLFitter::ResDoubleCB_3::GetAlphaL(double x) {
  return evaluate(x, 8);
}

double KLFitter::ResDoubleCB_3::GetAlphaR(double x) {
  return evaluate(x, 12);
}

double KLFitter::ResDoubleCB_3::GetNL(double x) {
  return fParameters.at(16);
}

double KLFitter::ResDoubleCB_3::GetNR(double x) {
  return fParameters.at(17);
}

double KLFitter::ResDoubleCB_3::evaluate(const double x, const size_t start) {
  double output = 0;
  for (size_t offset = 0; offset <= 3; ++offset) {
    const size_t idx = start + offset;
    output += fParameters.at(idx) / std::pow(x, idx);
  }
  return output;
}
