#include "DCHdigi.h"

// ROOT
#include "TVector3.h"

// STL
#include <iostream>
#include <sstream>

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       DCHdigi constructor       ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file", {"default name for the collection"}),
DCHdigi::DCHdigi(const std::string& name, ISvcLocator* svcLoc)
: MultiTransformer(name, svcLoc,
    {
        KeyValues("DCH_simhits", {""}),
        KeyValues("HeaderName", {"EventHeader"}),
    },
    {
        KeyValues("DCH_DigiCollection", {"DCH_DigiCollection"})
    }
  )
{
    m_geoSvc = serviceLocator()->service(m_geoSvcName);
    m_uidSvc = serviceLocator()->service(m_uidSvcName);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi::initialize() {

    if (!m_uidSvc)
    ThrowException( "Unable to get UniqueIDGenSvc" );

    m_gauss_z_cm  = std::normal_distribution<double>(0., m_z_resolution.value()*MM_TO_CM );
    m_gauss_xy_cm = std::normal_distribution<double>(0., m_xy_resolution.value()*MM_TO_CM);

    //-----------------
    // Retrieve the subdetector
    std::string DCH_name(m_DCH_name.value());
    if ( 0 == m_geoSvc->getDetector()->detectors().count(DCH_name) )
    {
        ThrowException( "Detector <<" + DCH_name + ">> does not exist." );
    }

    dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  retrieve data extension     //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    this->dch_data = DCH_DE.extension<dd4hep::rec::DCH_info>();

    ///////////////////////////////////////////////////////////////////////////////////

    //-----------------
    // Retrieve the readout associated with the detector element (subdetector)
    dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(DCH_name);
    if(not dch_sd.isValid() )
    ThrowException("No valid Sensitive Detector was found for detector <<" + DCH_name + ">>.");

    dd4hep::Readout dch_readout = dch_sd.readout();
    // set the cellID decoder
    m_decoder = dch_readout.idSpec().decoder();
    //-----------------


    std::stringstream ss;
    PrintConfiguration(ss);
    info() << ss.str().c_str() <<endmsg;
    if( m_create_debug_histos.value() )
    {
        hDpw = new TH1D("hDpw", "Distance hit to the wire, in cm", 100,0,1);
        hDpw->SetDirectory(0);
        hDww = new TH1D("hDww", "Distance hit projection to the wire, in cm. Should be zero", 100,0,1);
        hDww->SetDirectory(0);
        hSz  = new TH1D("hSz", "Smearing along the wire, in cm", 100,0,5*m_z_resolution.value());
        hSz->SetDirectory(0);
        hSxy = new TH1D("hSxy", "Smearing perpendicular the wire, in cm", 100,0,5*m_xy_resolution.value());
        hSxy->SetDirectory(0);
    }
    return StatusCode::SUCCESS;
  }



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::tuple<colltype_out>
DCHdigi::operator()(const colltype_in& input_sim_hits,
    const edm4hep::EventHeaderCollection&  headers) const {

    // initialize seed for random engine
    this->PrepareRandomEngine( headers );

    debug() << "Input Sim Hit collection size: " << input_sim_hits.size() << endmsg;

    // Create the collections we are going to return
    colltype_out output_digi_hits;

    //loop over hit collection
    for (const auto& input_sim_hit : input_sim_hits) {
        dd4hep::DDSegmentation::CellID cellid = input_sim_hit.getCellID();
        int ilayer = this->CalculateLayerFromCellID(cellid );
        int nphi   = this->CalculateNphiFromCellID(cellid );
        TVector3 hit_position { input_sim_hit.getPosition()[0]*MM_TO_CM ,
                                input_sim_hit.getPosition()[1]*MM_TO_CM ,
                                input_sim_hit.getPosition()[2]*MM_TO_CM };


        // -------------------------------------------------------------------------
        //      calculate hit position projection into the wire
        TVector3 hit_to_wire_vector = this->Calculate_hitpos_to_wire_vector(ilayer, nphi,hit_position);
        TVector3 hit_projection_on_the_wire = hit_position + hit_to_wire_vector;
        if( m_create_debug_histos.value() )
        {
            double distance_hit_wire = hit_to_wire_vector.Mag();
            hDpw->Fill(distance_hit_wire);
            // the distance from the hit projection and the wire should be zero
            // TVector3 dummy_vector = this->Calculate_hitpos_to_wire_vector(ilayer, nphi,hit_projection_on_the_wire);
            // hDww->Fill( dummy_vector.Mag() );
        }
        TVector3 wire_direction_ez = this->Calculate_wire_vector_ez(ilayer, nphi);

        // -------------------------------------------------------------------------
        //       smear the position
        double smearing_z = m_gauss_z_cm( m_engine );
        if( m_create_debug_histos.value() ) hSz->Fill( smearing_z );

        hit_projection_on_the_wire +=  smearing_z*(wire_direction_ez.Unit());
        if( m_create_debug_histos.value() )
        {
            // the distance from the hit projection and the wire should be zero
            TVector3 dummy_vector = this->Calculate_hitpos_to_wire_vector(ilayer, nphi,hit_projection_on_the_wire);
            hDww->Fill( dummy_vector.Mag() );
        }

        double smearing_xy = m_gauss_xy_cm( m_engine );
        if( m_create_debug_histos.value() ) hSxy->Fill( smearing_xy );


        std::int32_t type = 0;
        std::int32_t quality = 0;
        float eDepError =0;
        // length units back to mm
        edm4hep::Vector3d positionSW = {hit_projection_on_the_wire.x()/MM_TO_CM,
                                        hit_projection_on_the_wire.y()/MM_TO_CM,
                                        hit_projection_on_the_wire.z()/MM_TO_CM };
        edm4hep::Vector3d directionSW =  {  wire_direction_ez.x()/MM_TO_CM,
                                            wire_direction_ez.y()/MM_TO_CM,
                                            wire_direction_ez.z()/MM_TO_CM };
        float distanceToWire = hit_to_wire_vector.Mag();
        std::uint32_t clusterCount = 0;
        std::uint32_t clusterSize = 0;

        auto & i = input_sim_hit;
        output_digi_hits.create( i.getCellID(),
            type,
            quality,
            i.getTime(),
            i.getEDep(),
            eDepError,
            positionSW,
            directionSW,
            distanceToWire,
            clusterCount,
            clusterSize
            );
        }                                                       // end loop over hit collection


    /////////////////////////////////////////////////////////////////
    return std::make_tuple<colltype_out>(std::move(output_digi_hits));
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi::finalize()
{
    if( m_create_debug_histos.value() )
    {
        std::unique_ptr<TFile> ofile{TFile::Open ( "dch_digi_alg_debug.root", "recreate" ) };
        ofile->cd();
        hDpw->Write();
        hDww->Write();
        hSxy->Write();
        hSz->Write();
    }

    return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       ThrowException       ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi::ThrowException(std::string s) const {
    error() << s.c_str()  << endmsg;
    throw std::runtime_error(s);
  }

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi::PrintConfiguration(std::ostream& io)
{
    io << "DCHdigi will use the following components:\n";
    io << "\tGeometry Service: "                  << m_geoSvcName.value().c_str()           << "\n";
    io << "\tDetector name: "                     << m_DCH_name.value().c_str()             << "\n";
    io << "\t\t|--Volume bitfield: "              << m_decoder->fieldDescription().c_str()  << "\n";
    io << "\t\t|--Number of layers: "             << dch_data->database.size()              << "\n";
    return;
}

void DCHdigi::PrepareRandomEngine(const edm4hep::EventHeaderCollection&  headers) const
{
    uint32_t evt_n = headers[0].getEventNumber();
    uint32_t run_n = headers[0].getRunNumber();
    size_t seed = m_uidSvc->getUniqueID(evt_n, run_n, this->name() );
    m_engine.seed(seed);
    // test random engine...
    m_engine.discard(10);
}

///////////////////////////////////////////////////////////////////////////////////////
/////       Ancillary functions for calculating the distance to the wire       ////////
///////////////////////////////////////////////////////////////////////////////////////
TVector3 DCHdigi::Calculate_wire_vector_ez(int ilayer, int nphi) const
{
    auto & l = this->dch_data->database.at(ilayer);

    // See original paper Hoshina et al, Computer Physics Communications 153 (2003) 3
    // eq. 2.9, for the definition of ez, vector along the wire

    // initialize some variables
    int stereosign = l.StereoSign();
    double rz0 = l.radius_sw_z0;
    double dphi = dch_data->twist_angle;
    // kappa is the same as in eq. 2.9
    double kappa = (1./dch_data->Lhalf)*tan(dphi/2);

    //--- calculating wire position
    // the points p1 and p2 correspond to the ends of the wire

    // point 1
    // double x1 = rz0; // m
    // double y1 = 0.; // m
    // double z1 = 0.; // m
    double x1 = rz0;                                        // m
    double y1 = -stereosign*rz0*kappa*dch_data->Lhalf;      // m
    double z1 = -dch_data->Lhalf;                           // m

    TVector3 p1 (x1,y1,z1);


    // point 2
    double x2 = rz0;                                        // m
    double y2 = stereosign*rz0*kappa*dch_data->Lhalf;       // m
    double z2 = dch_data->Lhalf;                            // m

    TVector3 p2 (x2,y2,z2);

    // calculate phi rotation of whole twisted tube, ie, rotation at z=0
    double phi_z0 = Calculate_wire_phi_z0(ilayer,nphi);
    p1.RotateZ(phi_z0);
    p2.RotateZ(phi_z0);

    //--- end calculating wire position

    return (p2-p1).Unit();

}

TVector3 DCHdigi::Calculate_wire_z0_point(int ilayer, int nphi) const
{
    auto & l = this->dch_data->database.at(ilayer);
    double rz0 = l.radius_sw_z0;
    TVector3 p1 (rz0,0,0);
    double phi_z0 = Calculate_wire_phi_z0(ilayer,nphi);
    p1.RotateZ(phi_z0);
    return p1;
}

// calculate phi rotation of whole twisted tube, ie, rotation at z=0
double DCHdigi::Calculate_wire_phi_z0(int ilayer, int nphi) const
{
    auto & l = this->dch_data->database.at(ilayer);
    int ncells = l.nwires/2;
    double phistep = TMath::TwoPi()/ncells;
    double phi_z0 = (nphi + 0.25*(l.layer%2))*phistep;
    return phi_z0;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Calculate vector from hit position to wire   /////////////////
///////////////////////////////////////////////////////////////////////////////////////
TVector3 DCHdigi::Calculate_hitpos_to_wire_vector(int ilayer, int nphi, const TVector3 & hit_position /*in cm*/) const
{
    // Solution distance from a point to a line given here:
    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    TVector3 n = this->Calculate_wire_vector_ez(ilayer, nphi);
    TVector3 a = this->Calculate_wire_z0_point (ilayer, nphi);
    // Remember using cm as natural units of DD4hep consistently!
    // TVector3 p {hit_position.x()*MM_TO_CM,hit_position.y()*MM_TO_CM,hit_position.z()*MM_TO_CM};

    TVector3 a_minus_p = a - hit_position;
    double a_minus_p_dot_n = a_minus_p.Dot( n );
    TVector3 scaled_n = a_minus_p_dot_n * n;
    //hit_to_wire_vector = a_minus_p - scaled_n;
    return (a_minus_p - scaled_n);
}
