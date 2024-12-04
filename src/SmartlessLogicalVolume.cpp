//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================

// Framework include files
#include <DDG4/Factories.h>
// Geant4 include files
#include <G4Material.hh>
#include <G4LogicalVolume.hh>

// Forward declarations
namespace  {  class SmartlessLogicalVolume;  }

/// Namespace example name of the user
namespace dd4hep   {

  /// Class member specialization to create a G4LogicalCrystalVolume
  /** Class member specialization to create a G4LogicalCrystalVolume
   *  
   *  \author  M.Frank
   *  \version 1.0
   *  \ingroup DD4HEP_SIMULATION
   */
  template <> G4LogicalVolume* 
  Geant4LogicalVolumeFactory<SmartlessLogicalVolume>::create(dd4hep::Detector& /* description */,
							      Volume      volume,
							      G4VSolid*   solid,
							      G4Material* material)
  {

    auto* ptr = new G4LogicalVolume(solid, material, volume.name());
    printout(ALWAYS,"SmartlessLogicalVolume",
	     "====> Created specialize logical volume [SmartlessLogicalVolume]: %s",
	     volume.name());
    return ptr;
  }
}      // End namespace dd4hep

DECLARE_GEANT4LOGICALVOLUME(SmartlessLogicalVolume)
