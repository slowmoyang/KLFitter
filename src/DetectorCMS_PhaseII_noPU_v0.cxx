#include "KLFitter/DetectorCMS_PhaseII_noPU_v0.h"
#include "KLFitter/ResolutionBase.h"
#include "KLFitter/ResDoubleCB_2.h"
#include "KLFitter/ResDoubleCB_3.h"

#include <cmath>
#include <iostream>


namespace KLFitter {
// ---------------------------------------------------------
DetectorCMS_PhaseII_noPU_v0::DetectorCMS_PhaseII_noPU_v0(std::string folder) : DetectorBase() {
  std::cout << "Using TF with crystal ball for jets ..." << std::endl;

  const std::string b_jet_dir = folder + "/" + "b-jet";
  const std::string light_jet_dir = folder + "/" + "light-jet";

  // energy resolution
  m_res_energy_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_3{(light_jet_dir + "/eta-0p00-1p50.txt").c_str()});
  m_res_energy_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_2{(light_jet_dir + "/eta-1p50-1p75.txt").c_str()});
  m_res_energy_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_2{(light_jet_dir + "/eta-1p75-2p15.txt").c_str()});
  m_res_energy_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_2{(light_jet_dir + "/eta-2p15-3p00.txt").c_str()});

  m_res_energy_b_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_3{(b_jet_dir + "/eta-0p00-1p50.txt").c_str()});
  m_res_energy_b_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_2{(b_jet_dir + "/eta-1p50-1p75.txt").c_str()});
  m_res_energy_b_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_2{(b_jet_dir + "/eta-1p75-2p15.txt").c_str()});
  m_res_energy_b_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleCB_2{(b_jet_dir + "/eta-2p15-3p00.txt").c_str()});

}

// ---------------------------------------------------------
DetectorCMS_PhaseII_noPU_v0::~DetectorCMS_PhaseII_noPU_v0() = default;

// ---------------------------------------------------------
ResolutionBase* DetectorCMS_PhaseII_noPU_v0::ResEnergyLightJet(double eta) {
  const double abs_eta = fabs(eta);

  if (abs_eta < m_jet_eta_bin_1) {
    return m_res_energy_light_jet_eta1.get();

  } else if (abs_eta < m_jet_eta_bin_2) {
    return m_res_energy_light_jet_eta2.get();

  } else if (abs_eta < m_jet_eta_bin_3) {
    return m_res_energy_light_jet_eta3.get();

  } else if (abs_eta <= m_jet_eta_bin_4) {
    return m_res_energy_light_jet_eta4.get();

  } else {
    std::cout << "DetectorCMS_PhaseII_noPU_v0::ResEnergyLightJet(). Eta range exceeded." << std::endl;
    return nullptr;

  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorCMS_PhaseII_noPU_v0::ResEnergyBJet(double eta) {
  const double abs_eta = fabs(eta);

  if (abs_eta < m_jet_eta_bin_1) {
    return m_res_energy_b_jet_eta1.get();

  } else if (abs_eta < m_jet_eta_bin_2) {
    return m_res_energy_b_jet_eta2.get();

  } else if (abs_eta < m_jet_eta_bin_3) {
    return m_res_energy_b_jet_eta3.get();

  } else if (abs_eta <= m_jet_eta_bin_4) {
    return m_res_energy_b_jet_eta4.get();

  } else {
    std::cout << "DetectorCMS_PhaseII_noPU_v0::ResEnergyBJet(). Eta range exceeded." << std::endl;
    return nullptr;

  }
}

}  // namespace KLFitter
