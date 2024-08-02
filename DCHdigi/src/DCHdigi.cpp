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
                            KeyValues("DCH_DigiCollection", {""})
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
