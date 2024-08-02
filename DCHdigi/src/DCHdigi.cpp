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
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi::initialize() {


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

	return StatusCode::SUCCESS;
}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::tuple<colltype_out>
DCHdigi::operator()(const colltype_in& input_sim_hits,
                         const edm4hep::EventHeaderCollection&  /*headers*/) const {


	debug() << "Input Sim Hit collection size: " << input_sim_hits.size() << endmsg;

	// Create the collections we are going to return
	colltype_out output_digi_hits;

 // // //    //loop over hit collection
	// // // for (const auto& input_sim_hit : input_sim_hits) {
 // // //      dd4hep::DDSegmentation::CellID id = input_sim_hit.getCellID();
 // // //      std::cout << "New ARC hit with cell ID: " << id << std::endl;
 // // //
 // // //
 // // //
 // // //      // emplace output objects
 // // //      // output_digi_hits.create(....);
 // // //
	// // // }// end loop over hit collection


	/////////////////////////////////////////////////////////////////
	return std::make_tuple<colltype_out>(std::move(output_digi_hits));
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi::finalize() { return StatusCode::SUCCESS; }

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

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Calculate vector from hit position to wire   /////////////////
///////////////////////////////////////////////////////////////////////////////////////
TVector3 DCHdigi::Calculate_hitpos_to_wire_vector(int ilayer, int nphi, const TVector3 & p) const
{
    auto & l = this->dch_data->database.at(ilayer);

    int stereosign = l.StereoSign();
    double rz0 = l.radius_sw_z0;
    double dphi = dch_data->twist_angle;
    double kappa = (1./dch_data->Lhalf)*tan(dphi/2);
    int ncells = l.nwires/2;
    double phistep = TMath::TwoPi()/ncells;
    double phi_z0 = (nphi + 0.25*(l.layer%2))*phistep;

    // point 1
    // double x1 = rz0; // m
    // double y1 = 0.; // m
    // double z1 = 0.; // m
    double x1 = rz0; // m
    double y1 = -stereosign*rz0*kappa*dch_data->Lhalf; // m
    double z1 = -dch_data->Lhalf; // m

    TVector3 p1 (x1,y1,z1);


    // point 2
    double x2 = rz0; // m
    double y2 = stereosign*rz0*kappa*dch_data->Lhalf; // m
    double z2 = dch_data->Lhalf; // m

    TVector3 p2 (x2,y2,z2);

    p1.RotateZ(phi_z0);
    p2.RotateZ(phi_z0);
    // Solution distance from a point to a line given here:
    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    TVector3 n = (p2-p1).Unit();
    TVector3 a = p1;
    // TVector3 p {hit.position.x()*MM_TO_CM,hit.position.y()*MM_TO_CM,hit.position.z()*MM_TO_CM};

    TVector3 a_minus_p = a - p;
    double a_minus_p_dot_n = a_minus_p.Dot( n );
    TVector3 scaled_n = a_minus_p_dot_n * n;
    //p_to_wire_vector = a_minus_p - scaled_n;
    return (a_minus_p - scaled_n);
}
