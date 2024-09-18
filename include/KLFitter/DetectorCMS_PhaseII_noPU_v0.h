#ifndef KLFITTER_DetectorCMS_PhaseII_noPU_v0_H_
#define KLFITTER_DetectorCMS_PhaseII_noPU_v0_H_

#include <memory>
#include <string>

#include "KLFitter/DetectorBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

class DetectorCMS_PhaseII_noPU_v0 : public DetectorBase {
 public:
  /**
   * The default constructor.
   * @param folder The folder with transfer function parameters.
   */
  explicit DetectorCMS_PhaseII_noPU_v0(std::string folder = "");

  /// The (defaulted) destructor.
  ~DetectorCMS_PhaseII_noPU_v0();

  ResolutionBase* ResEnergyLightJet(double eta = 0.) override;
  ResolutionBase* ResEnergyBJet(double eta = 0.) override;

  ResolutionBase* ResEnergyGluonJet(double /*eta*/) {return nullptr;}
  ResolutionBase* ResEnergyElectron(double /*eta*/) {return nullptr;}
  ResolutionBase* ResEnergyMuon(double /*eta*/) {return nullptr;}


  ResolutionBase* ResEnergyPhoton(double /*eta*/) {return nullptr;}
  ResolutionBase* ResMissingET() {return nullptr;}
  ResolutionBase* ResEtaLightJet(double /*eta*/) {return nullptr;}
  ResolutionBase* ResEtaBJet(double /*eta*/) {return nullptr;}
  ResolutionBase* ResPhiLightJet(double /*eta*/) {return nullptr;}
  ResolutionBase* ResPhiBJet(double /*eta*/) {return nullptr;}

 private:
  /// The energy resolution of light jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta4;

  /// The energy resolution of b jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_b_jet_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_b_jet_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_b_jet_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_b_jet_eta4;

  /// The eta binning for jets
  const double m_jet_eta_bin_1{1.50};
  const double m_jet_eta_bin_2{1.75};
  const double m_jet_eta_bin_3{2.15};
  const double m_jet_eta_bin_4{3.00};

};
}  // namespace KLFitter

#endif  // KLFITTER_DetectorCMS_PhaseII_noPU_v0_H_
