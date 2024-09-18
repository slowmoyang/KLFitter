#include "TreeReaderTTAllJets.cxx"

#include "KLFitter/Fitter.h"
#include "KLFitter/DetectorCMS_PhaseII_noPU_v0.h"
#include "KLFitter/LikelihoodTopAllHadronic.h"
#include "KLFitter/Permutations.h"

#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"

#include <iostream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <set>


void run(const TString input_path,
         const TString detector_dir,
         const TString output_path,
         const bool is_matched_tt,
         const int max_jets) {

  auto input_file = TFile::Open(input_path, "READ");
  if (input_file == nullptr) {
    throw std::runtime_error("cannot open the input ROOT file");
  }
  auto input_tree = dynamic_cast<TTree*>(input_file->Get("tree"));
  TreeReaderTTAllJets reader{input_tree};

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////
  auto output_file = TFile::Open(output_path, "RECREATE");
  auto output_tree = input_tree->CloneTree(0);
  output_tree->SetDirectory(output_file);

  const std::vector<std::string> fit_status_vec = {
      "total",
      "FitAbortedDueToNaNMask",
      "MinuitDidNotConvergeMask",
      "AtLeastOneFitParameterAtItsLimitMask",
      "InvalidTransferFunctionAtConvergenceMask",
  };

  TH1D h_fit_status{"fit_status", "", static_cast<int>(fit_status_vec.size()), -0.5, fit_status_vec.size() - 0.5};
  h_fit_status.SetDirectory(output_file);
  for (size_t idx = 0; idx < fit_status_vec.size(); ++idx) {
      const int bin = static_cast<int>(idx) + 1;
      const std::string label = fit_status_vec.at(idx);
      h_fit_status.GetXaxis()->SetBinLabel(bin, label.c_str());
  }

  std::vector<std::string> float_branch_name_vec = {
      "nll",
      "event_probability",
  };
  const std::vector<std::string> gen_part_name_vec = {
      "t0", "b0", "w0", "q00", "q01",
      "t1", "b1", "w1", "q10", "q11",
  };
  const std::vector<std::string> component_vec = {
      "pt", "eta", "phi", "mass",
  };
  const std::vector<std::string> prefix_vec = {
      "rec", "fit"
  };
  for (const auto& prefix : prefix_vec) {
    for (const auto& gen_part_name : gen_part_name_vec) {
      for (const auto& component : component_vec) {
        const auto name = prefix + "_" + gen_part_name + "_" + component;
        float_branch_name_vec.push_back(name);
      }
    }
  }

  std::vector<std::string> int_branch_name_vec = {
      "num_permutations",
  };

  std::vector<std::string> bool_branch_name_vec = {
      "correct_b0",
      "correct_b1",
      "correct_w0",
      "correct_w1",
      "correct_t0",
      "correct_t1",
      "correct",
      "found",
  };

  std::map<std::string, float> float_branch_map;
  std::cout << "  - float:" << std::endl;
  for (const auto& each : float_branch_name_vec) {
    std::cout << "    - " << each << std::endl;
    float_branch_map.emplace(each, 0.0f);
    output_tree->Branch(each.c_str(), &float_branch_map.at(each));
  }


  std::map<std::string, int> int_branch_map;
  std::cout << "  - int:" << std::endl;
  for (const auto& each : int_branch_name_vec) {
    std::cout << "    - " << each << std::endl;
    int_branch_map.emplace(each, false);
    output_tree->Branch(each.c_str(), &int_branch_map.at(each));
  }


  std::map<std::string, bool> bool_branch_map;
  std::cout << "  - bool:" << std::endl;
  for (const auto& each : bool_branch_name_vec) {
    std::cout << "    - " << each << std::endl;
    bool_branch_map.emplace(each, false);
    output_tree->Branch(each.c_str(), &bool_branch_map.at(each));
  }

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////
  KLFitter::Fitter fitter{};

  KLFitter::DetectorCMS_PhaseII_noPU_v0 detector{detector_dir.Data()};
  if (not fitter.SetDetector(&detector)) {
    throw std::runtime_error("failed to set detector");
  }

  KLFitter::LikelihoodTopAllHadronic likelihood{};
  likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kVetoHybridNoFit);
  likelihood.SetFlagTopMassFixed(false);
  if (not fitter.SetLikelihood(&likelihood)) {
    throw std::runtime_error("Failed to set likelihood");
  }

  /////////////////////////////////////////////////////////////////////////////
  // analyse
  /////////////////////////////////////////////////////////////////////////////
  const std::vector<std::string> parton_name_vec = {"b0", "b1", "q00", "q01", "q10", "q11"};

    size_t count_correct_t0 = 0;
    size_t count_correct_t1 = 0;
    size_t count_correct_t = 0;
    size_t count_correct_tt = 0;


  // const auto total = input_tree->GetEntries();
  const int total = 1000;
  for (int idx = 0; idx < total; ++idx) {
    std::cout << "Event: " << idx + 1 << " out of " << total << std::endl;


    //#########################################################################
    // data loading
    //#########################################################################
    const Event event = reader.getEntry(idx);

    //#########################################################################
    // preprocessing
    //#########################################################################
    KLFitter::Particles particles{};
    for (int jet_idx = 0; jet_idx < max_jets; jet_idx++) {
      const TLorentzVector& jet = event.jet_vec.at(jet_idx);
      particles.AddParticle(
          /*particle=*/jet,
          /*DetEta=*/jet.Eta(),
          /*ptype=*/KLFitter::Particles::kParton,
          /*name=*/"",
          /*measuredindex=*/jet_idx,
          /*isBtagged=*/event.jet_btag_vec.at(jet_idx)
      );
    }

    // NOTE Add particles to the likelihood.
    if (not fitter.SetParticles(&particles)) {
      throw std::runtime_error("failed to add particles to KLFitter");
    }

    //#########################################################################
    // fit
    //#########################################################################
    double min_nll = 1e9;
    double best_event_probability = -1;
    bool found = false;

    // predicted partons with fitted momenta
    std::map<std::string, TLorentzVector> fit_parton_map{};
    // predicted partons' jet indices
    std::map<std::string, int> rec_parton_jet_idx_map{};


    const int num_permutations = fitter.Permutations()->NPermutations();
    for (int permutation_idx  = 0; permutation_idx < num_permutations; permutation_idx++) {
      fitter.Fit(permutation_idx);

      //-----------------------------------------------------------------------
      // fit diagnosis
      //-----------------------------------------------------------------------
      const unsigned int status = fitter.ConvergenceStatus();

      h_fit_status.Fill(0);
      if ((status & fitter.FitAbortedDueToNaNMask) != 0) {
        h_fit_status.Fill(1);
        continue;
      }
      if ((status & fitter.MinuitDidNotConvergeMask) != 0) {
        h_fit_status.Fill(2);
        continue;
      }
      if ((status & fitter.AtLeastOneFitParameterAtItsLimitMask) != 0) {
        h_fit_status.Fill(3);
        continue;
      }
      if ((status & fitter.InvalidTransferFunctionAtConvergenceMask) != 0) {
        h_fit_status.Fill(4);
        continue;
      }

      //-----------------------------------------------------------------------
      //
      //-----------------------------------------------------------------------
      const double nll = -fitter.Likelihood()->LogLikelihood(fitter.Likelihood()->GetBestFitParameters());

      if ((not found) or (nll < min_nll)) {
        found = true;
        min_nll = nll;
        auto best_model_particles = fitter.Likelihood()->ParticlesModel();
        auto best_permuted_particles = fitter.Likelihood()->PParticlesPermuted();
        best_event_probability = std::exp(fitter.Likelihood()->LogEventProbability());

        for (size_t parton_idx = 0; parton_idx < parton_name_vec.size(); ++parton_idx) {
          const auto parton_name = parton_name_vec.at(parton_idx);
          fit_parton_map[parton_name] = *((*best_model_particles).Parton(parton_idx));
        }

        for (size_t parton_idx = 0; parton_idx < parton_name_vec.size(); ++parton_idx) {
          const auto parton_name = parton_name_vec.at(parton_idx);
          const auto jet_idx = (**best_permuted_particles).JetIndex(parton_idx);
          rec_parton_jet_idx_map[parton_name] = jet_idx;
        }
      }
    } // permutation loop

    //#########################################################################
    // postprocessing
    //#########################################################################

    // reset
    for (auto it = float_branch_map.begin(); it != float_branch_map.end(); ++it) {
      it->second = 0;
    }
    for (auto it = int_branch_map.begin(); it != int_branch_map.end(); ++it) {
      it->second = 0;
    }
    for (auto it = bool_branch_map.begin(); it != bool_branch_map.end(); ++it) {
      it->second = false;
    }

    std::cout << "best NLL: " << min_nll << std::endl;
    std::cout << "found: " << found << std::endl;
    if (found) {
      //-----------------------------------------------------------------------
      // postprocessing
      //-----------------------------------------------------------------------
      fit_parton_map.emplace("w0", fit_parton_map.at("q00") + fit_parton_map.at("q01"));
      fit_parton_map.emplace("w1", fit_parton_map.at("q10") + fit_parton_map.at("q11"));
      fit_parton_map.emplace("t0", fit_parton_map.at("b0") + fit_parton_map.at("w0"));
      fit_parton_map.emplace("t1", fit_parton_map.at("b1") + fit_parton_map.at("w1"));

      // predicted partons at the detector (or reconstruction) level
      std::map<std::string, TLorentzVector> rec_parton_map{};
      // predicted partons' truth labels
      std::map<std::string, int> label_map{};
      for (size_t parton_idx = 0; parton_idx < parton_name_vec.size(); ++parton_idx) {
        const auto parton_name = parton_name_vec.at(parton_idx);
        const int jet_idx = rec_parton_jet_idx_map.at(parton_name);

        rec_parton_map.emplace(parton_name, event.jet_vec.at(jet_idx));
        label_map.emplace(parton_name, event.jet_label_vec.at(jet_idx));
      }

      rec_parton_map.emplace("w0", rec_parton_map.at("q00") + rec_parton_map.at("q01"));
      rec_parton_map.emplace("w1", rec_parton_map.at("q10") + rec_parton_map.at("q11"));
      rec_parton_map.emplace("t0", rec_parton_map.at("b0") + rec_parton_map.at("w0"));
      rec_parton_map.emplace("t1", rec_parton_map.at("b1") + rec_parton_map.at("w1"));

      // q00 and q01 are assinged to same w
      const bool correct_w0 = (label_map.at("q00") == label_map.at("q01")) and ((label_map.at("q00") == 1) or (label_map.at("q00") == 4));
      // q10 and q11 are assinged to same w
      const bool correct_w1 = (label_map.at("q10") == label_map.at("q11")) and ((label_map.at("q10") == 1) or (label_map.at("q10") == 4));

      // b0 is assinged to b or b~
      const bool correct_b0 = (label_map.at("b0") == 0) or (label_map.at("b0") == 3);
      // b1 is assinged to b or b~
      const bool correct_b1 = (label_map.at("b1") == 0) or (label_map.at("b1") == 3);

      const bool correct_t0 = correct_b0 and correct_w0 and ((label_map.at("b0") + 1) == label_map.at("q00"));
      const bool correct_t1 = correct_b1 and correct_w1 and ((label_map.at("b1") + 1) == label_map.at("q10"));

      const bool correct = correct_t0 and correct_t1;

      if (correct_t0) count_correct_t0++;
      if (correct_t1) count_correct_t1++;
      if (correct_t0 or correct_t1) count_correct_t++;
      if (correct_t0 and correct_t1) count_correct_tt++;

      for (size_t parton_idx = 0; parton_idx < parton_name_vec.size(); ++parton_idx) {
        const auto parton_name = parton_name_vec.at(parton_idx);

        std::cout << "  - " << parton_name << ": " << label_map.at(parton_name) << std::endl;
      }

      std::cout << "  - correct_t0: " << correct_t0 << std::endl;
      std::cout << "  - correct_t0: " << correct_t1 << std::endl;
      std::cout << "  - correct: " << correct << std::endl;
      std::cout << "  - processed: " << (idx + 1)
                << ", correct_t0: " << count_correct_t0
                << ", correct_t1: " << count_correct_t1
                << ", correct_t: " << count_correct_t
                << ", correct_tt: " << count_correct_tt
                << std::endl;


      //-----------------------------------------------------------------------
      // fill
      //-----------------------------------------------------------------------
      for (auto it = fit_parton_map.begin(); it != fit_parton_map.end(); ++it) {
        const std::string name = it->first;
        const TLorentzVector gen_part = it->second;
        const auto prefix = "fit_" + name + "_";
        float_branch_map.at(prefix + "pt") = gen_part.Pt();
        float_branch_map.at(prefix + "eta") = gen_part.Eta();
        float_branch_map.at(prefix + "phi") = gen_part.Phi();
        float_branch_map.at(prefix + "mass") = gen_part.M();
      }

      for (auto it = rec_parton_map.begin(); it != rec_parton_map.end(); ++it) {
        const std::string name = it->first;
        const TLorentzVector gen_part = it->second;

        const auto prefix = "rec_" + name + "_";

        float_branch_map.at(prefix + "pt") = gen_part.Pt();
        float_branch_map.at(prefix + "eta") = gen_part.Eta();
        float_branch_map.at(prefix + "phi") = gen_part.Phi();
        float_branch_map.at(prefix + "mass") = gen_part.M();
      }

      float_branch_map.at("nll") = min_nll;
      float_branch_map.at("event_probability") = best_event_probability;

      int_branch_map.at("num_permutations") = num_permutations;

      bool_branch_map.at("correct_b0") = correct_b0;
      bool_branch_map.at("correct_b1") = correct_b1;
      bool_branch_map.at("correct_w0") = correct_w0;
      bool_branch_map.at("correct_w1") = correct_w1;
      bool_branch_map.at("correct_t0") = correct_t0;
      bool_branch_map.at("correct_t1") = correct_t1;
      bool_branch_map.at("correct") = correct;
      bool_branch_map.at("found") = found;
    }

    output_tree->Fill();
  } // event loop

  output_tree->Print();
  output_file->Write();
  output_file->Close();
  input_file->Close();
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 6) {
      throw std::runtime_error("usage: tt-alljets input_path detector_dir output_path is_matched_tt max_jets");
    }

    const TString input_path{argv[1]};
    const TString detector_dir{argv[2]};
    const TString output_path{argv[3]};
    const bool is_matched_tt = std::atoi(argv[4]) == 1;
    const int max_jets = std::atoi(argv[5]);

    run(
        /*input_path=*/input_path,
        /*detector_dir=*/detector_dir,
        /*output_path=*/output_path,
        /*is_matched_tt=*/is_matched_tt,
        /*max_jets=*/max_jets
    );

  } catch (const std::exception& err) {
    std::cerr << "😱😱😱: " << err.what() << std::endl;
    return 1;

  }

  return 0;
}
