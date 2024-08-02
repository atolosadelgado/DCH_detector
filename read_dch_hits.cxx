
#include <string>
std::string ddsim_output_file_str = "dch_proton_10GeV.root"; /*"dch_proton_10GeV.root"; "dch_proton_10GeV_limitsON.root"*/

std::string dch_compact_file_str = "./compact/DCH_standalone_o1_v02.xml";
std::string dch_detector_name_str = "DCH_v2";
std::string dch_encoding_str =  "system:5,superlayer:5,layer:4,nphi:11,stereosign:-1";

bool draw_wires = true;

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "vector"
// #include "edm4hep/SimTrackerHitData.h"
// #include "podio/ObjectID.h"
// #include "edm4hep/EventHeaderData.h"
// #include "edm4hep/MCParticleData.h"
// // use alias to simplify the code
// using ArcData_t = edm4hep::SimTrackerHitData;
// using MCParticle_t = edm4hep::MCParticleData;

#include "DDSegmentation/BitFieldCoder.h"
#include "DDRec/DCH_info.h"

// DD4hep
#include "DD4hep/Detector.h"
//-- for cellID-position conversion
#include "DDRec/CellIDPositionConverter.h"
//-- for cellID decoding
#include "DDSegmentation/BitFieldCoder.h"

// link ROOT interpreter with DD4hep libraries
R__LOAD_LIBRARY(DDG4IO);
R__LOAD_LIBRARY(DDCore);
R__LOAD_LIBRARY(DDRec);
#include "DDG4/Geant4Data.h"
#include <DDG4/Geant4Particle.h>

#include <cmath>
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TPolyLine3D.h"


TH3D * hXYZ = new TH3D("hXYZ","", 100,-2000, 2000, 100,-2000, 2000, 100,-2000, 2000 );
TH1D * hDpw = new TH1D("hDpw", "Distance hit to the wire, in cm", 100,0,1);
TH1D * hDpw2 = new TH1D("hDpw2", "Distance hit to the wire, in cm", 100,0,1);

using dd4hep::rec::DCH_info;
int DistanceFromPointPToCentralSWireOfDCHcell(int ilayer, int nphi, DCH_info * dch_data, const TVector3 & p, TVector3 & p_to_wire_vector)
{
    auto & l = dch_data->database.at(ilayer);

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
    p_to_wire_vector = a_minus_p - scaled_n;
    return 0;
}

struct mywire_t
{
   double L = {0}; //m
   double twisted_angle = {0}; //radians
   TVector3 p1;
   TVector3 p2;
   TVector3 wireDirection();
   TVector3 wirePoint();
   int      CalculateDistanceFromPointPToCentralSWire(const TVector3 & p, TVector3 & p_to_wire_vector);
   void set(double rz0, int stereosign, double angle);
};

void mywire_t::set(double rz0, int stereosign, double angle)
{
   double dphi = twisted_angle;
   double kappa = 2/L*tan(dphi/2);

   // point 1
   // double x1 = rz0; // m
   // double y1 = 0.; // m
   // double z1 = 0.; // m
   double x1 = rz0; // m
   double y1 = -stereosign*rz0*kappa*L/2; // m
   double z1 = -L/2; // m

   p1 = TVector3(x1,y1,z1);


   // point 2
   double x2 = rz0; // m
   double y2 = stereosign*rz0*kappa*L/2; // m
   double z2 = L/2; // m

   p2 = TVector3(x2,y2,z2);

   p1.RotateZ(angle);
   p2.RotateZ(angle);


}

TVector3 mywire_t::wireDirection()
{
    return (p2-p1).Unit();
}

TVector3 mywire_t::wirePoint()
{
    return p1;
}

int mywire_t::CalculateDistanceFromPointPToCentralSWire(const TVector3 & p, TVector3& p_to_wire_vector)
{
    // Solution distance from a point to a line given here:
    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    TVector3 n = this->wireDirection();
    TVector3 a = this->wirePoint();
    // TVector3 p {hit.position.x()*MM_TO_CM,hit.position.y()*MM_TO_CM,hit.position.z()*MM_TO_CM};

    TVector3 a_minus_p = a - p;
    double a_minus_p_dot_n = a_minus_p.Dot( n );
    TVector3 scaled_n = a_minus_p_dot_n * n;
    p_to_wire_vector = a_minus_p - scaled_n;
    return 0;
}





TPolyLine3D *sw_polyline;
TCanvas *c1 = new TCanvas();
//watch out units, dd4hep natural units is cm, but edm4hep is mm, stay in cm
constexpr double MM_TO_CM = 0.1;
void read_dch_hits_edm4hep();
void read_dch_hits_dd4hep();


void read_dch_hits()
{
    // read_dch_hits_edm4hep();
    read_dch_hits_dd4hep();
}

void read_dch_hits_dd4hep()
{

    // configure the histogram
    c1->Draw();
    hXYZ->SetDirectory(0);
    hXYZ->SetMarkerColor(2);
    hXYZ->SetMarkerStyle(20);
    hXYZ->SetMarkerSize(0.7);
    hXYZ->Draw();

    const auto & theDetector = & ( dd4hep::Detector::getInstance() );
    theDetector->fromXML ( dch_compact_file_str );
    dd4hep::DetElement DCH_DE = theDetector->detector ( dch_detector_name_str );
    dd4hep::rec::DCH_info * dch_data = DCH_DE.extension<dd4hep::rec::DCH_info>();
    // int ilayer_max = (--dch_data->database.end())->first;
    // int ilayer_min = dch_data->database.begin()->first;

    dd4hep::DDSegmentation::BitFieldCoder decoder ( dch_encoding_str );

    std::unique_ptr<TFile> ifile{TFile::Open ( ddsim_output_file_str.c_str() , "read" ) };
    if ( !ifile || ifile->IsZombie() ) {
        std::cerr << "Problem opening the file testFile.root" << std::endl;
        return;
    }

    TTreeReader aReader ( "EVENT", ifile.get() );
    TTreeReaderValue<std::vector<dd4hep::sim::Geant4Tracker::Hit*> >    DCHCollection ( aReader, "DCHCollection" );
    TTreeReaderValue<std::vector<dd4hep::sim::Geant4Particle*> > MCParticleCollection ( aReader, "MCParticles" );

    //watch out units, dd4hep natural units is cm, but edm4hep is mm, stay in cm
    auto ff_fdw_r = [&](double z, dd4hep::rec::DCH_info_struct::DCH_info_layer & l){
                    double r_z0 = l.radius_fdw_z0;
                    double alpha = dch_data->stereoangle_z0(r_z0);
                    double r2 = pow(r_z0, 2) + pow( z*tan(alpha),2);
                    return sqrt(r2);
                };
    auto ff_fuw_r = [&](double z, dd4hep::rec::DCH_info_struct::DCH_info_layer & l){
                    double r_z0 = l.radius_fuw_z0;
                    double alpha = dch_data->stereoangle_z0(r_z0);
                    double r2 = pow(r_z0, 2) + pow( z*tan(alpha),2);
                    return sqrt(r2);
                };
    mywire_t awire;
    awire.L = 2*(dch_data->Lhalf);
    awire.twisted_angle =  dch_data->twist_angle;
            
    // Loop over the entries of the tree
    while ( aReader.Next() ) {

        std::cout << Form("Event number %lld", aReader.GetCurrentEntry() ) << std::endl;

        for ( auto & hit_ptr : ( *DCHCollection ) ) {

            auto & hit = *hit_ptr;

            hXYZ->Fill(hit.position.x(), hit.position.y(), hit.position.z());
            std::cout << Form("%g\t%g\t%g",hit.position.x()*MM_TO_CM, hit.position.y()*MM_TO_CM, hit.position.z()*MM_TO_CM ) << std::endl;
            double r = sqrt( pow(hit.position.x(),2) + pow(hit.position.y(), 2) )*MM_TO_CM;

            // lets calcualte the radius of the corresponding ilayer
            dd4hep::DDSegmentation::CellID id = hit.cellID;
            int ilayer = decoder.get ( id,"layer" ) + 8*decoder.get ( id,"superlayer" ) + 1;
            auto & l = dch_data->database.at(ilayer);
            // result in cm, natural units of dd4hep
            double r_z_up   = ff_fuw_r(hit.position.z()*MM_TO_CM, l);
            double r_z_down = ff_fdw_r(hit.position.z()*MM_TO_CM, l);
            bool is_radius_correct = (r_z_down < r && r< r_z_up);

            // if the radius is not between the up/down radius of the layer, print error
            if( not is_radius_correct)
            {
                std::cout << "Error, point not contained within the layer radii: " << r_z_down << "\t" << r << "\t" << r_z_up << std::endl;
            }

            auto ff_phi_z = [&](double z_over_Lhalf ){ return atan (z_over_Lhalf * tan( dch_data->twist_angle/2./dd4hep::rad /* 15*TMath::DegToRad()*/ )); };
            int ncells = l.nwires/2;
            double phistep = TMath::TwoPi()/ncells;
            int nphi = decoder.get ( id,"nphi" );
            double phi_z0 = (nphi + 0.25*(l.layer%2))*phistep;

            double sim_phi = atan2( hit.position.y(), hit.position.x() );

            int stereosign = l.StereoSign(); //decoder.get(id, "stereosign");
            double calculated_phi = phi_z0  + stereosign*ff_phi_z( hit.position.z()*MM_TO_CM / dch_data->Lhalf );


            calculated_phi = std::remainder(calculated_phi,TMath::TwoPi());
            sim_phi = std::remainder(sim_phi,TMath::TwoPi());

            double phi_pull = (calculated_phi - sim_phi)/phistep;

            // if phi distance is bigger than phi step (or pull>1), means error
            if( 1 < fabs(phi_pull) )
            {
                std::cout << "Error calculating phi in layer " << ilayer << ", pull = " << phi_pull << std::endl;
            }
            // std::cout << ilayer << "\t" << nphi << std::endl;
            // std::cout << ilayer << "\t" << calculated_phi/phistep << "\t" << sim_phi/phistep << "\t" << l.StereoSign() << std::endl;




            awire.set( l.radius_sw_z0 , stereosign , phi_z0);
            TVector3 p {hit.position.x()*MM_TO_CM,hit.position.y()*MM_TO_CM,hit.position.z()*MM_TO_CM};
            TVector3 p_to_wire;
            awire.CalculateDistanceFromPointPToCentralSWire(p, p_to_wire);
            double distance_hit_wire = p_to_wire.Mag();
            std::cout << Form("Distance to the wire: %g", distance_hit_wire ) << std::endl;
            hDpw->Fill(distance_hit_wire);

            TVector3 rvector;
            DistanceFromPointPToCentralSWireOfDCHcell(ilayer, nphi, dch_data,p,rvector);
            hDpw2->Fill(rvector.Mag());


            if( draw_wires )
            {
                // Draw wire
                sw_polyline = new TPolyLine3D(2);
                sw_polyline->SetPoint(0,awire.p1.x()/dd4hep::mm, awire.p1.y()/dd4hep::mm, awire.p1.z()/dd4hep::mm);
                sw_polyline->SetPoint(1,awire.p2.x()/dd4hep::mm, awire.p2.y()/dd4hep::mm, awire.p2.z()/dd4hep::mm);
                c1->cd();
                sw_polyline->Draw();
            }




        }
        // if drawing wires, process only 1 event
        if( draw_wires ) break;

    }
}



// void read_dch_hits_edm4hep()
// {
//
//     c1->Draw();
//     hXYZ->SetDirectory(0);
//     hXYZ->Draw();
//     const auto & theDetector = & ( dd4hep::Detector::getInstance() );
//     theDetector->fromXML ( "k4geo/lcgeoTests/compact/DCH_standalone_o1_v02.xml" );
//     dd4hep::DetElement DCH_DE = theDetector->detector ( "DCH_v2" );
//     dd4hep::rec::DCH_info * dch_data = DCH_DE.extension<dd4hep::rec::DCH_info>();
//     int ilayer_max = (--dch_data->database.end())->first;
//     int ilayer_min = dch_data->database.begin()->first;
//
//     dd4hep::DDSegmentation::BitFieldCoder decoder (  "system:5,superlayer:5,layer:4,nphi:11,stereosign:-1" );
//
//     //std::unique_ptr<TFile> ifile{TFile::Open ( "dch_geantino.root", "read" ) };
//     std::unique_ptr<TFile> ifile{TFile::Open ( "dch_proton_20GeV.root", "read" ) };
//     if ( !ifile || ifile->IsZombie() ) {
//         std::cerr << "Problem opening the file testFile.root" << std::endl;
//         return;
//     }
//
//     // Create a TTreeReader to read the tree named "myTree"
//     TTreeReader aReader ( "events", ifile.get() );
//
//     // Create TTreeReaderValues for the branches "branch1" and "branch2"
//     TTreeReaderValue<std::vector<ArcData_t> >    ArcCollection ( aReader, "DCHCollection" );
//     // TTreeReaderValue<std::vector<MCParticle_t> > MCParticleCollection ( aReader, "ArcCollection" );
//
//     //watch out units, dd4hep natural units is cm, but edm4hep is mm, stay in mm
//     auto ff_fdw_r = [&](double z, dd4hep::rec::DCH_info_struct::DCH_info_layer & l){
//                     double r_z0 = l.radius_fdw_z0;
//                     double alpha = dch_data->stereoangle_z0(r_z0);
//                     double r2 = pow(r_z0, 2) + pow( z*tan(alpha),2);
//                     return sqrt(r2);
//                 };
//     auto ff_fuw_r = [&](double z, dd4hep::rec::DCH_info_struct::DCH_info_layer & l){
//                     double r_z0 = l.radius_fuw_z0;
//                     double alpha = dch_data->stereoangle_z0(r_z0);
//                     double r2 = pow(r_z0, 2) + pow( z*tan(alpha),2);
//                     return sqrt(r2);
//                 };
//     mywire_t awire;
//     awire.L = 2*(dch_data->Lhalf);
//     awire.twisted_angle =  dch_data->twist_angle;
//
//     // Loop over the entries of the tree
//     while ( aReader.Next() ) {
//         //std::cout << DCHCollection->size() << std::endl;
//         for ( auto & hit : ( *ArcCollection ) ) {
//
//             hXYZ->Fill(hit.position.x, hit.position.y, hit.position.z);
//             // std::cout << Form("%g\t%g\t%g",hit.position.x, hit.position.y, hit.position.z) << std::endl;
//             double r = sqrt( pow(hit.position.x,2) + pow(hit.position.y, 2) )*MM_TO_CM;
//             // aproach 1
//             // lets compare layer number
//             // // find layer:
//             // int hitlayer = -1;
//             // for( auto & [ilayer, l] : (dch_data->database) )
//             // {
//             //
//             //     double r_fuw = ff_fuw_r(hit.position.z, l);
//             //     double r_fdw = ff_fdw_r(hit.position.z, l);
//             //     // std::cout << hit.position.z << "\t" << l.radius_fuw_z0/dd4hep::mm << std::endl;
//             //     if( r_fdw < r && r_fuw > r)
//             //     {
//             //         hitlayer = ilayer;
//             //         break;
//             //     }
//             // }
//             // dd4hep::DDSegmentation::CellID id = hit.cellID;
//             // int ilayer = decoder.get ( id,"layer" ) + 8*decoder.get ( id,"superlayer" ) + 1;
//             // std::cout << ilayer - hitlayer << std::endl;
//
//             // aproach 2, lets calcualte the radius of the corresponding ilayer
//             // lets assume ilayer is ok...
//             dd4hep::DDSegmentation::CellID id = hit.cellID;
//             int ilayer = decoder.get ( id,"layer" ) + 8*decoder.get ( id,"superlayer" ) + 1;
//             auto & l = dch_data->database.at(ilayer);
//             // result in cm, natural units of dd4hep
//             double r_z_up   = ff_fuw_r(hit.position.z*MM_TO_CM, l);
//             double r_z_down = ff_fdw_r(hit.position.z*MM_TO_CM, l);
//             // std::cout << r_z_down << "\t" << r << "\t" << r_z_up << std::endl;
//             // std::cout << (r_z_down < r && r< r_z_up) << std::endl;
//
//             auto ff_phi_z = [&](double z_over_Lhalf ){ return atan (z_over_Lhalf * tan( 15*TMath::DegToRad() )); };
//             int ncells = l.nwires/2;
//             double phistep = TMath::TwoPi()/ncells;
//             int nphi = decoder.get ( id,"nphi" );
//             double phi_z0 = (nphi + 0.25*(l.layer%2))*phistep;
//             // double phi_z0_minus = (nphi -1 + 0.25*(l.layer%2))*phistep;
//             // double phi_z0_plus  = (nphi +1 + 0.25*(l.layer%2))*phistep;
//             // if(nphi == ncells-1)
//             // {
//             //     phi_z0_plus = (0 + 0.25*(l.layer%2))*phistep;
//             // }
//             // if(nphi == 0)
//             // {
//             //     phi_z0_minus = (ncells -1 + 0.25*(l.layer%2))*phistep;
//             // }
//             double sim_phi = atan2( hit.position.y, hit.position.x ); //atan( hit.position.y/hit.position.x) + 3.1415*( (hit.position.x < 0) && (hit.position.y < 0));
//
//             int stereosign = l.StereoSign(); //decoder.get(id, "stereosign");
//             double calculated_phi = phi_z0  + stereosign*ff_phi_z( hit.position.z*MM_TO_CM / dch_data->Lhalf );
//
//
//             calculated_phi = std::remainder(calculated_phi,TMath::TwoPi());
//             sim_phi = std::remainder(sim_phi,TMath::TwoPi());
//
//             double phi_pull = (calculated_phi - sim_phi)/phistep;
//
//             // std::cout << ilayer << "\t" << phi_pull << std::endl;
//             // std::cout << ilayer << "\t" << nphi << std::endl;
//             std::cout << ilayer << "\t" << calculated_phi/phistep << "\t" << sim_phi/phistep << "\t" << l.StereoSign() << std::endl;
//
//
//
//
//             awire.set( l.radius_sw_z0 , stereosign , phi_z0);
//
//             sw_polyline = new TPolyLine3D(2);
//             sw_polyline->SetPoint(0,awire.p1.x()/dd4hep::mm, awire.p1.y()/dd4hep::mm, awire.p1.z()/dd4hep::mm);
//       	    sw_polyline->SetPoint(1,awire.p2.x()/dd4hep::mm, awire.p2.y()/dd4hep::mm, awire.p2.z()/dd4hep::mm);
//       	    c1->cd();
//       	    sw_polyline->Draw();
//
//
//         }
//         break;
//
//     }
// }
