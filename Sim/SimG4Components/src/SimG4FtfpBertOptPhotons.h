#ifndef SIMG4COMPONENTS_G4FTFPBERTOPTPHOTONS_H
#define SIMG4COMPONENTS_G4FTFPBERTOPTPHOTONS_H
// Gaudi
#include "GaudiKernel/AlgTool.h"

// FCCSW
#include "SimG4Interface/ISimG4PhysicsList.h"

/** @class SimG4FtfpBert SimG4Components/src/SimG4FtfpBert.h SimG4FtfpBert.h
 *
 *  FTFP_BERT physics list tool.
 * + scecialised detailed optical photons.
 *  @author Anna Zaborowska, Coralie Neubüser
 */

class SimG4FtfpBertOptPhotons : public AlgTool, virtual public ISimG4PhysicsList {
public:
  explicit SimG4FtfpBertOptPhotons(const std::string& aType, const std::string& aName, const IInterface* aParent);
  virtual ~SimG4FtfpBertOptPhotons();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize();
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize();
  /** Get the physics list.
   *  @return pointer to G4VModularPhysicsList (ownership is transferred to the caller)
   */
  virtual G4VModularPhysicsList* physicsList();
};

#endif /* SIMG4COMPONENTS_G4FTFPBERTOPTPHOTONS_H */
