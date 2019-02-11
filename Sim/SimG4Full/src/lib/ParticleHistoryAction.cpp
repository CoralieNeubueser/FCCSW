#include "SimG4Full/ParticleHistoryAction.h"

#include "SimG4Common/EventInformation.h"

#include "G4EventManager.hh"
#include "G4LorentzVector.hh"

namespace sim {

ParticleHistoryAction::ParticleHistoryAction(double aEnergyCut, bool selectTaggedOnly)
    : m_energyCut(aEnergyCut), m_selectTaggedOnly(selectTaggedOnly) {}

void ParticleHistoryAction::PreUserTrackingAction(const G4Track* aTrack) {
  // copy track
  auto startTrack = new G4Track(*aTrack);
  auto trackID = aTrack->GetTrackID();
  // set ids, which are not copied
  startTrack->SetTrackID(trackID);
  startTrack->SetParentID(aTrack->GetParentID());
  // add it to map
  m_trackMap[trackID] = startTrack;
}

void ParticleHistoryAction::PostUserTrackingAction(const G4Track* aTrack) {
  auto g4EvtMgr = G4EventManager::GetEventManager();
  auto evtinfo = dynamic_cast<sim::EventInformation*>(g4EvtMgr->GetUserInformation());
  if (selectParticle(*aTrack, m_energyCut, m_selectTaggedOnly)) {
    auto search = m_trackMap.find(aTrack->GetTrackID());
    if (search != m_trackMap.end()) {
      evtinfo->addParticle(search->second);
    } else {
      throw std::runtime_error("ParticleHistoryaction::G4Track entry does not exist in map");
    }
  }
}

bool ParticleHistoryAction::selectParticle(const G4Track& aTrack, double aEnergyCut, bool selectTaggedOnly) {
  G4LorentzVector p4(aTrack.GetMomentum(), aTrack.GetTotalEnergy());
  if (p4.e() < aEnergyCut) {
    return false;
  }
  // in case only tagged tracks, with a user information set, should be selected
  if (selectTaggedOnly && (!aTrack.GetUserInformation())) {
    return false;
  }
  // all other cases
  return true;
}
}
