/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

#include "KLFitter/ResDoubleCBBase.h"

// Needed for CrystalBall
#include "Math/Math.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResDoubleCBBase::ResDoubleCBBase(const char * filename, const int npar) :
  KLFitter::ResolutionBase(npar) {
  // read parameters from file
  ReadParameters(filename, npar);
}

// ---------------------------------------------------------
KLFitter::ResDoubleCBBase::ResDoubleCBBase(std::vector<double> const& parameters) :
  KLFitter::ResolutionBase(parameters) {
}

// ---------------------------------------------------------
KLFitter::ResDoubleCBBase::~ResDoubleCBBase() = default;


// ---------------------------------------------------------
double KLFitter::ResDoubleCBBase::logp(double x, double xmeas, bool *good, double /*par*/) {
  double mean = GetMean(x);
  double sigma = GetSigma(x);
  double alpha_left = GetAlphaL(x);
  double alpha_right = GetAlphaR(x);
  double n_left = GetNL(x);
  double n_right = GetNR(x);

  // sanity checks for n and sigma
  *good = CheckCrystalBallSanity(&sigma, &n_left, &n_right);

  const double dx = (x - xmeas) / x;

  x_.setVal(dx);
  mean_.setVal(mean);
  sigma_.setVal(sigma);
  alpha_left_.setVal(alpha_left);
  n_left_.setVal(n_left);
  alpha_right_.setVal(alpha_right);
  n_right_.setVal(n_right);

  const double prob = pdf_.getVal(RooArgSet(x_));
  return std::log(prob);
}
