#ifndef TREEREADERTTALLJETS_H_
#define TREEREADERTTALLJETS_H_

#include "TBranch.h"
#include "TTree.h"
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

struct Event {
  std::vector<TLorentzVector> jet_vec;
  std::vector<int> jet_btag_vec;
  std::vector<int> jet_label_vec;


  Event(std::vector<TLorentzVector> _jet_vec,
        std::vector<int> _jet_btag_vec,
        std::vector<int> _jet_label_vec)

    : jet_vec(_jet_vec),
      jet_btag_vec(_jet_btag_vec),
      jet_label_vec(_jet_label_vec) {

  }
};





class TreeReaderTTAllJets {
 public:
  TreeReaderTTAllJets(TTree*);
  ~TreeReaderTTAllJets();

  Event getEntry(long);


 private:
  TTree* tree_;

  // branches

  std::vector<float>* jet_pt_;
  std::vector<float>* jet_eta_;
  std::vector<float>* jet_phi_;
  std::vector<float>* jet_mass_;
  std::vector<int>* jet_btag_;
  std::vector<int>* jet_parton_idx_;
  bool is_jet_parton_matched_;

  TBranch* b_jet_pt_;
  TBranch* b_jet_eta_;
  TBranch* b_jet_phi_;
  TBranch* b_jet_mass_;
  TBranch* b_jet_btag_;
  TBranch* b_jet_parton_idx_;
  TBranch* b_is_jet_parton_matched_;


};


TreeReaderTTAllJets::TreeReaderTTAllJets(TTree* tree)
  : tree_(tree) {
  jet_pt_ = nullptr;
  jet_eta_ = nullptr;
  jet_phi_ = nullptr;
  jet_mass_ = nullptr;
  jet_btag_ = nullptr;
  jet_parton_idx_ = nullptr;
  is_jet_parton_matched_ = false;

  tree_->SetMakeClass(1);
  tree_->SetBranchAddress("jet_pt", &jet_pt_, &b_jet_pt_);
  tree_->SetBranchAddress("jet_eta", &jet_eta_, &b_jet_eta_);
  tree_->SetBranchAddress("jet_phi", &jet_phi_, &b_jet_phi_);
  tree_->SetBranchAddress("jet_mass", &jet_mass_, &b_jet_mass_);
  tree_->SetBranchAddress("jet_btag", &jet_btag_, &b_jet_btag_);
  tree_->SetBranchAddress("jet_parton_idx", &jet_parton_idx_, &b_jet_parton_idx_);
  tree_->SetBranchAddress("is_jet_parton_matched", &is_jet_parton_matched_, &b_is_jet_parton_matched_);

}


TreeReaderTTAllJets::~TreeReaderTTAllJets() {
  /* if (tree_ != nullptr) { */
  /*   delete tree_->GetCurrentFile(); */
  /*   delete tree_; */
  /* } */
}







Event TreeReaderTTAllJets::getEntry(const long entry) {
  jet_pt_->clear();
  jet_eta_->clear();
  jet_phi_->clear();
  jet_mass_->clear();
  jet_btag_->clear();
  jet_parton_idx_->clear();
  is_jet_parton_matched_ = false;

  tree_->GetEntry(entry);

  const size_t jet_vec_size = jet_pt_->size();

  std::vector<TLorentzVector> jet_vec;
  jet_vec.reserve(jet_vec_size);
  for (size_t idx = 0; idx < jet_vec_size; ++idx) {
    TLorentzVector jet;
    jet.SetPtEtaPhiM(
        jet_pt_->at(idx),
        jet_eta_->at(idx),
        jet_phi_->at(idx),
        jet_mass_->at(idx)
    );
    jet_vec.push_back(jet);
  }

  std::vector<int> jet_label{};
  jet_label.reserve(jet_parton_idx_->size());
  // b0:  0 --> 0
  // q00: 1 --> 1
  // q01: 2 --> 1
  // b1:  3 --> 3
  // q10: 4 --> 4
  // q11: 5 --> 4
  // other: -1 -> -1

  for (const int parton_idx : *jet_parton_idx_) {
    int label = parton_idx;
    if (parton_idx == 2) {
      label = 1;

    } else if (parton_idx == 5) {
      label = 4;

    }

    jet_label.push_back(label);
  }

  return Event{
    jet_vec,
    *jet_btag_,
    jet_label
  };
}




#endif  // TREEREADERTTALLJETS_H_
