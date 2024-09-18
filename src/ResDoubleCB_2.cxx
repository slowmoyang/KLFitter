#include "KLFitter/ResDoubleCB_2.h"

KLFitter::ResDoubleCB_2::ResDoubleCB_2(const char * filename) :
  KLFitter::ResDoubleCBBase(filename, (3 * 4) + (1 * 2)) {
}

KLFitter::ResDoubleCB_2::~ResDoubleCB_2() = default;

double KLFitter::ResDoubleCB_2::GetMean(double x) {
  return evaluate(x, 0);
}

double KLFitter::ResDoubleCB_2::GetSigma(double x) {
  return evaluate(x, 3);
}

double KLFitter::ResDoubleCB_2::GetAlphaL(double x) {
  return evaluate(x, 6);
}

double KLFitter::ResDoubleCB_2::GetAlphaR(double x) {
  return evaluate(x, 9);
}

double KLFitter::ResDoubleCB_2::GetNL(double x) {
  return fParameters.at(12);
}

double KLFitter::ResDoubleCB_2::GetNR(double x) {
  return fParameters.at(13);
}

double KLFitter::ResDoubleCB_2::evaluate(const double x, const size_t start) {
  double output = 0;
  for (size_t offset = 0; offset <= 2; ++offset) {
    const size_t idx = start + offset;
    output += fParameters.at(idx) / std::pow(x, idx);
  }
  return output;
}
