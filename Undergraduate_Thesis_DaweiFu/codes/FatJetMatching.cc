#include "DeepNTuples/FatJetHelpers/interface/FatJetMatching.h"

#include <unordered_set>
#include "TString.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace deepntuples;


std::pair<FatJetMatching::FatJetLabel,const reco::GenParticle*> FatJetMatching::higgs_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{

  auto higgs = getFinal(parton);

  if (debug_){
    using namespace std;
    cout << "jet: " << jet->polarP4() << endl;
    cout << "H:   "; printGenParticleInfo(higgs, -1);
  }

  bool is_hVV = false;
  if (higgs->numberOfDaughters() >= 3) {
    // e.g., h->Vqq or h->qqqq
    is_hVV = true;
  }else {
    // e.g., h->VV*
    for (const auto &p : higgs->daughterRefVector()){
      auto pdgid = std::abs(p->pdgId());
      if (pdgid == ParticleID::p_Wplus || pdgid == ParticleID::p_Z0){
        is_hVV = true;
        break;
      }
    }
  }

  if (is_hVV){
    // h->WW or h->ZZ
    std::vector<const reco::GenParticle*> hVV_daus;
    int n_el = 0, n_mu = 0, n_tau = 0, n_quarks = 0;
    bool is_lhe_hvv = false;
    for (unsigned idau=0; idau<higgs->numberOfDaughters(); ++idau){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(idau));
      auto pdgid = std::abs(higgs->daughter(idau)->pdgId());
      if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) {
        ++n_quarks;
        hVV_daus.push_back(dau);
      }else if (pdgid >= ParticleID::p_eminus && pdgid <= ParticleID::p_nu_tau){
        if (pdgid == ParticleID::p_eminus) ++n_el;
        else if (pdgid == ParticleID::p_muminus) ++n_mu;
        else if (pdgid == ParticleID::p_tauminus) ++n_tau;
        hVV_daus.push_back(dau);
      }else{
        // const auto d = getDaughterQuarks(getFinal(dau));
        // hVV_daus.insert(hVV_daus.end(), d.begin(), d.end());
        is_lhe_hvv = true;
        auto daufinal = getFinal(dau);
        for (unsigned j=0; j<daufinal->numberOfDaughters(); ++j){
          const auto *ddau = dynamic_cast<const reco::GenParticle*>(daufinal->daughter(j));
          auto dpdgid = std::abs(ddau->pdgId());
          if (dpdgid >= ParticleID::p_d && dpdgid <= ParticleID::p_b){
            ++n_quarks;
            hVV_daus.push_back(ddau);
          }else if (dpdgid >= ParticleID::p_eminus && dpdgid <= ParticleID::p_nu_tau){
            if (dpdgid == ParticleID::p_eminus) ++n_el;
            else if (dpdgid == ParticleID::p_muminus) ++n_mu;
            else if (dpdgid == ParticleID::p_tauminus) ++n_tau;
            hVV_daus.push_back(ddau);
          }
        }
      }
    }
    auto is_neutrino = [&](int pdgid){
      pdgid = std::abs(pdgid);
      return pdgid == ParticleID::p_nu_e || pdgid == ParticleID::p_nu_mu || pdgid == ParticleID::p_nu_tau;
    };
    // auto printDaus = [&](std::vector<const reco::GenParticle*>& parts, const pat::Jet* jet){
    //   using namespace std;
    //   cout << "Found " << parts.size() << " quarks from Higgs decay" << endl;
    //   for (const auto * gp : parts){
    //     using namespace std;
    //     printGenParticleInfo(gp, -1);
    //     cout << " ... dR(q, jet) = " << reco::deltaR(*gp, *jet) << endl;
    //   }
    // };
    // if four daughters are all quarks
    if (n_quarks == 4){
      int n_daus_in_jet = 0, idx = 0;
      std::vector<int> idx_daus_in_jet;
      for (const auto *gp : hVV_daus){
        auto dr = reco::deltaR(*gp, *jet);
        if (dr < distR){
          ++n_daus_in_jet;
          idx_daus_in_jet.push_back(idx);
        }
        ++idx;
      }
      if (n_daus_in_jet >= 4){
        // std::cout << "4q!" << std::endl;
        return std::make_pair(FatJetLabel::H_ww4q, higgs);
      }
      else if (n_daus_in_jet >= 3){
        // std::cout << "3q!" << std::endl;
        return std::make_pair(FatJetLabel::H_ww3q, higgs);
      }
      else if (n_daus_in_jet >= 2 && is_lhe_hvv){
        if (idx_daus_in_jet == std::vector<int>({0,1}) || idx_daus_in_jet == std::vector<int>({2,3})){
          // std::cout << "2q from same W!" << std::endl;
          return std::make_pair(FatJetLabel::H_ww2qsame, higgs);
        }else {
          // std::cout << "2q from separate Ws!" << std::endl;
          return std::make_pair(FatJetLabel::H_ww2qsep, higgs);
        }
      }
    }else if (n_quarks == 2 && n_tau == 0) { // lvqq
      int n_daus_in_jet = 0;
      for (const auto *gp : hVV_daus){
        if (is_neutrino(gp->pdgId())){
          continue;
        }
        auto dr = reco::deltaR(*gp, *jet);
        if (dr < distR){
          ++n_daus_in_jet;
        }
      }
      if (n_daus_in_jet >= 3){
        // std::cout << "lvqq!" << std::endl;
        if (n_el > 0){
          return std::make_pair(FatJetLabel::H_wwevqq, higgs);
        }else if (n_mu > 0){
          return std::make_pair(FatJetLabel::H_wwmvqq, higgs);
        }
      }
    }else if (n_quarks == 2 && n_tau == 1) { // tauvqq
      int tau_pos = 0;
      int n_daus_in_jet = 0;
      for(std::size_t igp = 0; igp < hVV_daus.size(); ++igp) {
        auto gp = hVV_daus[igp];
        if (is_neutrino(gp->pdgId())){
          continue;
        }
        if (std::abs(gp->pdgId()) == ParticleID::p_tauminus){
          tau_pos = igp;
          continue;
        }
        auto dr = reco::deltaR(*gp, *jet);
        if (dr < distR){
          ++n_daus_in_jet;
        }
      }
      auto tau = getFinal(hVV_daus[tau_pos]);
      if (n_daus_in_jet >= 2){ // both two quarks are in
        bool is_leptau = false;
        for (unsigned j=0; j<tau->numberOfDaughters(); ++j){
          const auto *tdau = dynamic_cast<const reco::GenParticle*>(tau->daughter(j));
          auto tdau_pdgid = std::abs(tdau->pdgId());
          if (tdau_pdgid == ParticleID::p_eminus || tdau_pdgid == ParticleID::p_muminus){
            is_leptau = true;
            auto dr = reco::deltaR(*tdau, *jet);
            if (dr < distR){
              // std::cout << "leptau!" << std::endl;
              if (tdau_pdgid == ParticleID::p_eminus){
                return std::make_pair(FatJetLabel::H_wwleptauevqq, higgs);
              }else if (tdau_pdgid == ParticleID::p_muminus){
                return std::make_pair(FatJetLabel::H_wwleptaumvqq, higgs);
              }
            }
          }
        }
        if (!is_leptau){ // hadronic taus
          auto dr = reco::deltaR(*tau, *jet);
          if (dr < distR){
            // std::cout << "hadtau!" << std::endl;
            return std::make_pair(FatJetLabel::H_wwhadtauvqq, higgs);
          }
        }
      }
    }else {
      // std::cout << "**undefined!" << std::endl;
      return std::make_pair(FatJetLabel::Invalid, higgs);
      // throw std::logic_error("[FatJetMatching::higgs_label] Illegal H->WW mode");
    }
    // std::cout << "**unmatch!" << std::endl;
    return std::make_pair(FatJetLabel::Invalid, higgs);

    // if (debug_){
    //   using namespace std;
    //   cout << "Found " << hVV_daus.size() << " quarks from Higgs decay" << endl;
    //   for (const auto * gp : hVV_daus){
    //     using namespace std;
    //     printGenParticleInfo(gp, -1);
    //     cout << " ... dR(q, jet) = " << reco::deltaR(*gp, *jet) << endl;
    //   }
    // }

    // unsigned n_quarks_in_jet = 0;
    // for (const auto *gp : hVV_daus){
    //   auto dr = reco::deltaR(*gp, *jet);
    //   if (dr < distR){
    //     ++n_quarks_in_jet;
    //   }
    // }
    // if (n_quarks_in_jet >= 4){
    //   return std::make_pair(FatJetLabel::H_qqqq, higgs);
    // }

  }else if (isHadronic(higgs)) {
    // direct h->qq

    auto hdaus = getDaughterQuarks(higgs);
    if (hdaus.size() < 2) throw std::logic_error("[FatJetMatching::higgs_label] Higgs decay has less than 2 quarks!");
//    if (zdaus.size() >= 2)
    {
      double dr_q1    = reco::deltaR(jet->p4(), hdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), hdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(hdaus.at(0), hdaus.at(1));
      }
      auto pdgid_q1 = std::abs(hdaus.at(0)->pdgId());
      auto pdgid_q2 = std::abs(hdaus.at(1)->pdgId());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        cout << "pdgid(q1)        : " << pdgid_q1 << endl;
        cout << "pdgid(q2)        : " << pdgid_q2 << endl;
      }

      if (dr_q1<distR && dr_q2<distR){
        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          return std::make_pair(FatJetLabel::H_bb, higgs);
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
          return std::make_pair(FatJetLabel::H_cc, higgs);
        }else {
          return std::make_pair(FatJetLabel::H_qq, higgs);
        }
      }
    }
  }else {
    // test h->tautau
    std::vector<const reco::GenParticle*> taus;
    for (unsigned i=0; i<higgs->numberOfDaughters(); ++i){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(i));
      if (std::abs(dau->pdgId()) == ParticleID::p_tauminus){
        taus.push_back(dau);
      }
    }
    if (taus.size() == 2){
      // higgs -> tautau
      // use first version or last version of the tau in dr?
      double dr_tau1    = reco::deltaR(jet->p4(), taus.at(0)->p4());
      double dr_tau2    = reco::deltaR(jet->p4(), taus.at(1)->p4());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, tau1)    : " << dr_tau1 << endl;
        cout << "deltaR(jet, tau2)    : " << dr_tau2 << endl;
      }

      auto isHadronicTau = [](const reco::GenParticle* tau){
        for (const auto &dau : tau->daughterRefVector()){
          auto pdgid = std::abs(dau->pdgId());
          if (pdgid==ParticleID::p_eminus || pdgid==ParticleID::p_muminus){
            return false;
          }
        }
        return true;
      };

      auto tau1 = getFinal(taus.at(0));
      auto tau2 = getFinal(taus.at(1));
      if (dr_tau1<distR && dr_tau2<distR){
        if (isHadronicTau(tau1) && isHadronicTau(tau2)) {
          return std::make_pair(FatJetLabel::H_tautau, higgs);
        }
      }
    }
  }

  return std::make_pair(FatJetLabel::Invalid, nullptr);

}


std::pair<FatJetMatching::FatJetLabel,const reco::GenParticle*> FatJetMatching::qcd_label(const pat::Jet* jet, const reco::GenParticleCollection& genParticles, double distR)
{

  const reco::GenParticle *parton = nullptr;
  double minDR = 999;
  for (const auto &gp : genParticles){
    if (gp.status() != 23) continue;
    auto pdgid = std::abs(gp.pdgId());
    if (!(pdgid<ParticleID::p_t || pdgid==ParticleID::p_g)) continue;
    auto dr = reco::deltaR(gp, *jet);
    if (dr<distR && dr<minDR){
      minDR = dr;
      parton = &gp;
    }
  }
  if (debug_){
    using namespace std;
    if (parton){
      cout << "parton"; printGenParticleInfo(parton, -1);
      cout << "dr(jet, parton): " << minDR << endl;
    }
  }

  auto n_bHadrons = jet->jetFlavourInfo().getbHadrons().size();
  auto n_cHadrons = jet->jetFlavourInfo().getcHadrons().size();

  if (n_bHadrons>=2) {
    return std::make_pair(FatJetLabel::QCD_bb, parton);
  }else if (n_bHadrons==1){
    return std::make_pair(FatJetLabel::QCD_b, parton);
  }else if (n_cHadrons>=2){
    return std::make_pair(FatJetLabel::QCD_cc, parton);
  }else if (n_cHadrons==1){
    return std::make_pair(FatJetLabel::QCD_c, parton);
  }

  return std::make_pair(FatJetLabel::QCD_others, parton);
}
