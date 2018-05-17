#ifndef DETSENSITIVE_AGGREGATECALORIMETERSD_H
#define DETSENSITIVE_AGGREGATECALORIMETERSD_H

// DD4hep
#include "DDG4/Geant4Hits.h"
#include "DDSegmentation/Segmentation.h"

// Geant
#include "G4THitsCollection.hh"
#include "G4VSensitiveDetector.hh"

/** AggregateBirksLawCalorimeterSD DetectorDescription/DetSensitive/src/AggregateBirksLawCalorimeterSD.h AggregateBirksLawCalorimeterSD.h
 *
 *  Sensitive detector for calorimeter (aggregates energy deposits within each cell).
 *  It is based on dd4hep::sim::Geant4GenericSD<Calorimeter> (but it is not identical).
 *  In particular, the position of the hit is set to G4Step::GetPreStepPoint() position.
 *  No timing information is saved (energy deposits are aggregated in the cells)
 *
 *  Birks law reduces the energy deposited in the scintillator, following:
 *
 *  Note : the material is assumed ideal, which means that impurities and aging effects are not taken into account
 *  edep = destep / (1. + RKB*dedx + C*(dedx)**2)
 *  the basic units of the coefficient are g/(MeV*cm**2) and de/dx is obtained in MeV/(g/cm**2)
 *  values from NIM 80 (1970) 239-244
 *  RKB = 0.013  g/(MeV*cm**2)  and  C = 9.6e-6  g**2/((MeV**2)(cm**4))
 *
 *  @author    Anna Zaborowska
 *  @author    Coralie Neubueser
 */

namespace det {
class AggregateBirksLawCalorimeterSD : public G4VSensitiveDetector {
public:
  /** Constructor.
   *  @param aDetectorName Name of the detector
   *  @param aReadoutName Name of the readout (used to name the collection)
   *  @param aSeg Segmentation of the detector (used to retrieve the cell ID)
   */
  AggregateBirksLawCalorimeterSD(const std::string& aDetectorName,
                         const std::string& aReadoutName,
                         const dd4hep::Segmentation& aSeg);
  /// Destructor
  virtual ~AggregateBirksLawCalorimeterSD();
  /** Initialization.
   *  Creates the hit collection with the name passed in the constructor.
   *  The hit collection is registered in Geant.
   *  @param aHitsCollections Geant hits collection.
   */
  virtual void Initialize(G4HCofThisEvent* aHitsCollections) final;
  /** Process hit once the particle hit the sensitive volume.
   *  Checks if the energy deposit is larger than 0, calculates the position and cellID,
   *  saves that into the hit collection.
   *  If there is already entry in the same cell, the energy is accumulated.
   *  Otherwise new hit is created.
   *  @param aStep Step in which particle deposited the energy.
   */
  virtual bool ProcessHits(G4Step* aStep, G4TouchableHistory*) final;

private:
  /// Collection of calorimeter hits
  G4THitsCollection<dd4hep::sim::Geant4CalorimeterHit>* m_calorimeterCollection;
  /// Segmentation of the detector used to retrieve the cell Ids
  dd4hep::Segmentation m_seg;
  // Variables needed for the calculation of birks law
  const std::string m_material;
  const double m_birk1;
  const double m_birk2;
};
}

#endif /* DETSENSITIVE_AGGREGATECALORIMETERSD_H */
