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
#include <random>

#include "DDRec/DCH_info.h"

#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"

/// constant to convert from mm (EDM4hep) to DD4hep (cm)
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


  //------------------------------------------------------------------
  //          machinery for geometry

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc>                        m_geoSvc;

  /// Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Pointer to drift chamber data extension
  dd4hep::rec::DCH_info * dch_data = {nullptr};

  //------------------------------------------------------------------
  //          machinery for smearing the position

  /// along the sense wire position resolution in mm
  Gaudi::Property<float> m_z_resolution{this, "zResolution", 10.0,
                               "Spatial resolution in the z direction (from reading out the wires at both sides) [mm]"};
  /// xy resolution in mm
  Gaudi::Property<float> m_xy_resolution{this, "xyResolution", 10., "Spatial resolution in the xy direction [mm]"};

  /// create seed using the uid
  SmartIF<IUniqueIDGenSvc>                   m_uidSvc;
  /// use thread local engine from C++ standard
  inline static thread_local std::mt19937_64 m_engine;


  // Operator std::normal_distribution<T>::operator()(Generator& g) is a non-const member function and thus cannot be called for a constant object. So we defined the distribution as mutable.
  // Gaussian random number generator used for the smearing of the z position, in cm!
  mutable std::normal_distribution<double> m_gauss_z_cm;
  // Gaussian random number generator used for the smearing of the xy position, in cm!
  mutable std::normal_distribution<double> m_gauss_xy_cm;


  //------------------------------------------------------------------
  //        ancillary functions

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io);

  /// Send error message to logger and then throw exception
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

  //------------------------------------------------------------------
  //        debug information
  /// Flag to create output file with debug histgrams
  Gaudi::Property<bool> m_create_debug_histos{this, "create_debug_histograms", true, "Create output file with histograms for debugging"};

  /// histogram to store distance from hit position to the wire
  TH1D * hDpw;

  /// histogram to store distance from hit projection to the wire (should be zero)
  TH1D * hDww;

  /// histogram to store smearing along the wire
  TH1D * hSz;

  /// histogram to store smearing perpendicular the wire
  TH1D * hSxy;



};

DECLARE_COMPONENT(DCHdigi);

#endif
