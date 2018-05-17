#include "DetSensitive/AggregateBirksLawCalorimeterSD.h"

// FCCSW
#include "DetCommon/DetUtils.h"

// DD4hep
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4VolumeManager.h"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// Geant4
#include "G4SDManager.hh"

namespace det {
AggregateBirksLawCalorimeterSD::AggregateBirksLawCalorimeterSD(const std::string& aDetectorName,
                                               const std::string& aReadoutName,
                                               const dd4hep::Segmentation& aSeg)
    : G4VSensitiveDetector(aDetectorName), 
      m_calorimeterCollection(nullptr), 
      m_seg(aSeg),
      // variables for birks law
      m_material("Polystyrene"),
      m_birk1(0.0130 * CLHEP::g / (CLHEP::MeV * CLHEP::cm2)),
      m_birk2(9.6e-6 * CLHEP::g / (CLHEP::MeV * CLHEP::cm2) * CLHEP::g / (CLHEP::MeV * CLHEP::cm2)) {
  // name of the collection of hits is determined byt the readout name (from XML)
    collectionName.insert(aReadoutName);
}

AggregateBirksLawCalorimeterSD::~AggregateBirksLawCalorimeterSD() {}

void AggregateBirksLawCalorimeterSD::Initialize(G4HCofThisEvent* aHitsCollections) {
  // create a collection of hits and add it to G4HCofThisEvent
  // deleted in ~G4Event
  m_calorimeterCollection =
      new G4THitsCollection<dd4hep::sim::Geant4CalorimeterHit>(SensitiveDetectorName, collectionName[0]);
  aHitsCollections->AddHitsCollection(G4SDManager::GetSDMpointer()->GetCollectionID(m_calorimeterCollection),
                                      m_calorimeterCollection);
}

bool AggregateBirksLawCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  // check if energy was deposited
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;

  // Corection of energy deposit following Birks law
  G4double response = 0.;
  G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
  G4double charge = aStep->GetPreStepPoint()->GetCharge();
  if ((charge != 0.) && (m_material.compare(material->GetName()) == 0)) {
    G4double rkb = m_birk1;
    // --- correction for particles with more than 1 charge unit ---
    // --- based on alpha particle data ---
    if (std::fabs(charge) > 1.0) rkb *= 7.2 / 12.6;

    if (aStep->GetStepLength() != 0) {
      G4double dedx = edep / (aStep->GetStepLength()) / (material->GetDensity());
      response = edep / (1. + rkb * dedx + m_birk2 * dedx * dedx);
    } else {
      response = edep;
    }
  } else {
    response = edep;
  }
  edep = response;

  // as in dd4hep::sim::Geant4GenericSD<Calorimeter>
  CLHEP::Hep3Vector prePos = aStep->GetPreStepPoint()->GetPosition();
  CLHEP::Hep3Vector postPos = aStep->GetPostStepPoint()->GetPosition();
  CLHEP::Hep3Vector midPos = 0.5 * (postPos + prePos);
  dd4hep::Position pos(midPos.x(), midPos.y(), midPos.z());
  // check the cell ID
  uint64_t id = utils::cellID(m_seg, *aStep);
  dd4hep::sim::Geant4CalorimeterHit* hit = nullptr;
  dd4hep::sim::Geant4CalorimeterHit* hitMatch = nullptr;
  // Check if there is already some energy deposit in that cell
  for (int i = 0; i < m_calorimeterCollection->entries(); i++) {
    hit = dynamic_cast<dd4hep::sim::Geant4CalorimeterHit*>(m_calorimeterCollection->GetHit(i));
    if (hit->cellID == id) {
      hitMatch = hit;
      hitMatch->energyDeposit += edep;
      return true;
    }
  }
  // if not, create a new hit
  // deleted in ~G4Event
  hitMatch = new dd4hep::sim::Geant4CalorimeterHit(pos);
  hitMatch->cellID = id;
  hitMatch->energyDeposit = edep;
  m_calorimeterCollection->insert(hitMatch);
  return true;
}
}
