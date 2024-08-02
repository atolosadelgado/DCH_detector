/** ======= DCHdigi ==========
 * Gaudi Algorithm for DCH digitization
 *
 *
 * @author Alvaro Tolosa-Delgado, Brieuc Francois
 * @date   2024-08
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires a collection of SimTrackerHits <br>
 * This code uses DD4hep length natural unit (cm), but EDM4hep data is (usually) in mm. Please be careful with units.  <br>
 * <h4>Output</h4>
 * Processor produces collection of ParticleID<br>
 * @param DCH_simhits The name of input collection, type edm4hep::SimTrackerHitCollection <br>
 * (default name empty) <br>
 * @param GeoSvcName Geometry service name <br>
 * (default value GeoSvc)
 * @param DCH_name ARC subdetector name <br>
 * (default value DCH_v2) <br>
 * <br>
 */

#ifndef DCHDIGI_H
#define DCHDIGI_H

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/EventHeaderCollection.h"

// EDM4HEP extension
#include "extension/DriftChamberDigiV2Collection.h"


#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// DD4hep
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"

#include <string>

#include "DDRec/DCH_info.h"

#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"

constexpr double MM_TO_CM = 0.1;

using colltype_in  = edm4hep::SimTrackerHitCollection;
using colltype_out = extension::DriftChamberDigiV2Collection;

struct DCHdigi final
    : k4FWCore::MultiTransformer<
          std::tuple<colltype_out>(
              const colltype_in&, const edm4hep::EventHeaderCollection&)> {
  DCHdigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<colltype_out> operator()(
      const colltype_in& ,
      const edm4hep::EventHeaderCollection&   ) const override;

private:

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc>                        m_geoSvc;

  /// Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io);

  /// Send error message to logger and throw exception
  void ThrowException(std::string s) const;

  int CalculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return m_decoder->get ( id,"layer" ) + dch_data->nlayersPerSuperlayer*m_decoder->get ( id,"superlayer" ) + 1;
  }

  int CalculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return m_decoder->get ( id,"nphi" );
  }

  // the following functions should be upstreamed to the data extension at DD4hep
  // to avoid code duplication and keep it centralized
  TVector3 Calculate_hitpos_to_wire_vector(int ilayer, int nphi, const TVector3 & hit_position /*in cm*/) const;
  TVector3 Calculate_wire_vector_ez       (int ilayer, int nphi) const;
  TVector3 Calculate_wire_z0_point        (int ilayer, int nphi) const;
  double   Calculate_wire_phi_z0          (int ilayer, int nphi) const;

  /// Declare here variables to be initialized at when creating the algorithm and then used within the event loop

  /// Pointer to drift chamber data extension
  dd4hep::rec::DCH_info * dch_data = {nullptr};

  TH1D * hDpw;


};

DECLARE_COMPONENT(DCHdigi);

#endif
