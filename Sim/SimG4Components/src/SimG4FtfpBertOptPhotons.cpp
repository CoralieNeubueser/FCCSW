#include "SimG4FtfpBertOptPhotons.h"

// Geant4
#include "FTFP_BERT.hh"
#include "G4VModularPhysicsList.hh"
#include "G4OpticalPhysics.hh"

DECLARE_TOOL_FACTORY(SimG4FtfpBertOptPhotons)

SimG4FtfpBertOptPhotons::SimG4FtfpBertOptPhotons(const std::string& aType, const std::string& aName, const IInterface* aParent)
    : AlgTool(aType, aName, aParent) {
  declareInterface<ISimG4PhysicsList>(this);
}

SimG4FtfpBertOptPhotons::~SimG4FtfpBertOptPhotons() {}

StatusCode SimG4FtfpBertOptPhotons::initialize() { return AlgTool::initialize(); }

StatusCode SimG4FtfpBertOptPhotons::finalize() { return AlgTool::finalize(); }

G4VModularPhysicsList* SimG4FtfpBertOptPhotons::physicsList() {
  auto pL = new FTFP_BERT;
  pL->RegisterPhysics( new G4OpticalPhysics );
  // ownership passed to SimG4Svc which will register it in G4RunManager. To be deleted in ~G4RunManager()
  return pL;
}
