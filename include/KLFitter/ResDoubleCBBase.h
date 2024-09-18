#ifndef KLFITTER_RESDOUBLECBBASE_H_
#define KLFITTER_RESDOUBLECBBASE_H_

#include <iostream>
#include <vector>

#include "KLFitter/ResolutionBase.h"

#include "RooRealVar.h"
#include "RooCrystalBall.h"


namespace KLFitter {
class ResDoubleCBBase : public ResolutionBase {
 public:
  explicit ResDoubleCBBase(const char * filename, const int npar);
  explicit ResDoubleCBBase(std::vector<double> const& parameters);

  virtual ~ResDoubleCBBase();

  virtual double GetMean(double x) = 0;
  virtual double GetAlphaL(double x) = 0;
  virtual double GetAlphaR(double x) = 0;
  virtual double GetNL(double x) = 0;
  virtual double GetNR(double x) = 0;

  double logp(double x, double xmeas, bool *good, double /*par*/ = 0) override;

  static bool CheckCrystalBallSanity(double *sigma,
                                     double *n_left,
                                     double *n_right) {
    static const double eps = 0.000000001;

    if (*sigma <= 0.) {
      *sigma = 0. + eps;
      return false;
    }

    if (*n_left <= 1.) {
      *n_left = 1. + eps;
      return false;
    }

    if (*n_right <= 1.) {
      *n_right = 1. + eps;
      return false;
    }


    return true;
  }

  // random varaible
  RooRealVar x_{"x", "", -2, 2};
  // parameter
  RooRealVar mean_{"mean", "", -10, 10};
  RooRealVar sigma_{"sigma", "", 1e-3, 1e2};
  RooRealVar alpha_left_{"alpha_left", "", 1e-3, 1e2};
  RooRealVar alpha_right_{"alpha_right", "", 1e-3, 1e2};
  RooRealVar n_left_{"n_left", "", 1 + 1e-3, 1e2};
  RooRealVar n_right_{"n_right", "", 1 + 1e-3, 1e2};

  RooCrystalBall pdf_{"model", "", x_, mean_, sigma_, alpha_left_, n_left_, alpha_right_, n_right_};

};
}  // namespace KLFitter

#endif  // KLFITTER_RESDOUBLECBBASE_H_
