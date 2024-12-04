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

#include <string>
#include <charconv> // For std::from_chars

// Forward declarations
namespace  {  class SmartlessLogicalVolume;  }

/// Namespace example name of the user
namespace dd4hep   {


  /// C++17 way of translating string into integer
  std::errc safeStringToInt(const std::string& str, int & value) {
      auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), value);
      return ec;
  }

  /** Class member specialization to create a a G4LogicalVolume specifying the Smartless parameter
   *  
   *  \author  A. Tolosa-Delgado
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

    const std::string default_prop = "";
    std::string volume_prop  = volume.getProperty("SetSmartless",default_prop);
    if( default_prop != volume_prop)
    {
      // try conversion from string to integer
      auto s = volume_prop;
      int Smartless;
      bool IsConversionDone = std::from_chars(s.data(), s.data() + s.size(), Smartless).ec == std::errc{} ;
      if ( IsConversionDone )
      {
        ptr->SetSmartless(Smartless);
        printout(ALWAYS,"SmartlessLogicalVolume",
	     "====> Created specialize logical volume [SmartlessLogicalVolume] %s with Smartless parameter %d",
	     volume.name(), Smartless);
      }

    }
    return ptr;
  }
}      // End namespace dd4hep

DECLARE_GEANT4LOGICALVOLUME(SmartlessLogicalVolume)
