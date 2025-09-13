#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TBranch.h"
#include "THStack.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <filesystem>
#include "TKey.h"
#include "THnSparse.h"

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameters
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

//#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;
namespace fs = std::filesystem;

int nBinsPlane = 4; //Number of planes
int nBinsRPC = 72; //Number of RPCs
int nBinsBoard = 936; //Number of LBs

//o2::ccdb::CcdbApi api; //CCDB API as global object
o2::mid::Mapping mapping; //MID mapping object to construct ccdb object

//Function to upload to ccdb for PbPb data -> only get events with centrality > 50%
//void createCCDBPbPb(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB, THnSparse *hTotLB, long int startValidity, long int endValidity, int runNumber, string period) {
void createCCDBPbPb(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB, THnSparse *hTotLB, long int startValidity, long int endValidity, int firstRun, int lastRun, string period, bool merge, int trackGoal) {
    bool debug = false;

    //Variables used in the loop
    int LB936 = 0; //LB 1-> 936
    int plane = 0; //plane
    int counter = 0; //temporary variable to hold number of LB

    /*TH1F *hTestTotCounts = new TH1F(("hTestTotCounts_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBP = new TH1F(("hTestBP_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestNBP = new TH1F(("hTestNBP_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBoth = new TH1F(("hTestBoth_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);

    TH1F *hTestTotCountsAfter = new TH1F(("hTestTotCountsAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBPAfter = new TH1F(("hTestBPAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestNBPAfter = new TH1F(("hTestNBPAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBothAfter = new TH1F(("hTestBothAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);*/

    TH1F *hTestTotCounts = new TH1F(("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBP = new TH1F(("hTestBP_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestNBP = new TH1F(("hTestNBP_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBoth = new TH1F(("hTestBoth_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);

    TH1F *hTestTotCountsAfter = new TH1F(("hTestTotCountsAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBPAfter = new TH1F(("hTestBPAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestNBPAfter = new TH1F(("hTestNBPAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBothAfter = new TH1F(("hTestBothAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);


    std::vector<o2::mid::ChEffCounter> counterVector; //Vector of efficiency counter structs
    struct o2::mid::ChEffCounter entry; //Struct for efficiency counters
    o2::mid::ChamberEfficiency effMap; //Chamber efficiency map

    //Get projection of LB counts after pt cut
    TH1D* totLBCountsProj = hTotLB->Projection(0); 
    TH1D* BothPlanesLBCountsProj = hFiredBothPlanesLB->Projection(0);
    TH1D* BPLBCountsProj = hFiredBPLB->Projection(0);
    TH1D* NBPLBCountsProj = hFiredNBPLB->Projection(0);

    //Before cuts
    cout << "Runs between: " << firstRun << " and: " << lastRun << " before cuts" << endl;
    cout << "Total counts for LB "  << totLBCountsProj->GetEntries() << endl;
    cout << "Total counts BP LB "  << BPLBCountsProj->GetEntries() << endl;
    cout << "Total counts NBP LB "  << NBPLBCountsProj->GetEntries() << endl;
    cout << "Total counts Both LB "  << BothPlanesLBCountsProj->GetEntries() << endl; 
    
    //vector of struct, push back the struct populated with counts
    for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
        for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
            for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                
                plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                hTestTotCounts->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LB936));
                hTestBP->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LB936));
                hTestNBP->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LB936));
                hTestBoth->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LB936));

                counter++;
            }
        }
    }

    counter = 0;

    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    //Only accept muon track if pt > 2 GeV/c
    hTotLB->GetAxis(2)->SetRange(15,150);
    hFiredBothPlanesLB->GetAxis(2)->SetRange(15,150);
    hFiredBPLB->GetAxis(2)->SetRange(15,150);
    hFiredNBPLB->GetAxis(2)->SetRange(15,150);

    //Get projection of LB counts after pt cut
    totLBCountsProj = hTotLB->Projection(0); 
    BothPlanesLBCountsProj = hFiredBothPlanesLB->Projection(0);
    BPLBCountsProj = hFiredBPLB->Projection(0);
    NBPLBCountsProj = hFiredNBPLB->Projection(0); 
 
    //After cuts
    cout << "Runs between: " << firstRun << " and: " << lastRun << " after cuts" << endl;
    cout << "Total counts for LB "  << totLBCountsProj->GetEntries() << endl;
    cout << "Total counts BP LB "  << BPLBCountsProj->GetEntries() << endl;
    cout << "Total counts NBP LB "  << NBPLBCountsProj->GetEntries() << endl;
    cout << "Total counts Both LB "  << BothPlanesLBCountsProj->GetEntries() << endl; 

    //vector of struct, push back the struct populated with counts
    for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
        for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
            for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                
                //Debug printout
                if (debug) {
                    cout << "det ID " << ide << " col " << icol << " line " << iline << " LB " << mapping.getBoardId(iline,icol,ide) << endl;
                    cout << "LB: " << mapping.getBoardId(iline,icol,ide) << "\t unique FEEID:" << o2::mid::detparams::makeUniqueFEEId(ide, icol, iline) << "\t";
                }
                
                plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                entry.deId =  uint8_t(ide);
                entry.columnId = uint8_t(icol);
                entry.lineId = uint8_t(iline);

                //entry.counts[0] = hFiredBPLB->GetBinContent(LB936); //BP
                //entry.counts[1] = hFiredNBPLB->GetBinContent(LB936); //NBP
                //entry.counts[2] = hFiredBothPlanesLB->GetBinContent(LB936); //Both
                //entry.counts[3] = hTotLB->GetBinContent(LB936); //Total

                entry.counts[0] = BPLBCountsProj->GetBinContent(LB936); //BP
                entry.counts[1] = NBPLBCountsProj->GetBinContent(LB936); //NBP
                entry.counts[2] = BothPlanesLBCountsProj->GetBinContent(LB936); //Both
                entry.counts[3] = totLBCountsProj->GetBinContent(LB936); //Total

                counterVector.push_back(entry);

                hTestTotCountsAfter->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LB936));
                hTestBPAfter->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LB936));
                hTestNBPAfter->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LB936));
                hTestBothAfter->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LB936));

                counter++;

                //cout << LB936 << "\t" << (BothPlanesLBCountsProj->GetBinContent(LB936)/totLBCountsProj->GetBinContent(LB936))*100 << "\t" << 
                //(BPLBCountsProj->GetBinContent(LB936)/totLBCountsProj->GetBinContent(LB936))*100 << "\t" << 
                //(NBPLBCountsProj->GetBinContent(LB936)/totLBCountsProj->GetBinContent(LB936))*100 << endl;
                
                //Debug printout
                if (debug) {
                    cout << LB936 << "\t both: " << BothPlanesLBCountsProj->GetBinContent(LB936) << "\t tot: " << totLBCountsProj->GetBinContent(LB936) << endl;
                } 
            }
        }
    }

    //Delete objects here to avoid saturating RAM
    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    TFile *fOut;
    if (merge) {
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str())) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str());
        }
        fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/" +to_string(trackGoal) + "/o2-mid-ChEffCounter_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    }
    else {
        fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/o2-mid-ChEffCounter_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    }
    fOut->WriteObjectAny(&counterVector, "std::vector<o2::mid::ChEffCounter>","ccdb-object");
    fOut->Close();

    //TFile *fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/details_run_" + to_string(runNumber) + ".root").c_str(),"RECREATE");
    TFile *fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/details_run_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    fTest->cd();
    hTestTotCounts->Write("tot");
    hTestBP->Write("BP");
    hTestNBP->Write("NBP");
    hTestBoth->Write("Both");
    hTestTotCountsAfter->Write("totAfter");
    hTestBPAfter->Write("BPAfter");
    hTestNBPAfter->Write("NBPAfter");
    hTestBothAfter->Write("BothAfter");
    fTest->Close();
} //end of createCCDB for PbPb

//////////

//Function to upload to ccdb
//void createCCDB(TH1F *hFiredBPLB, TH1F *hFiredNBPLB, TH1F *hFiredBothPlanesLB, TH1F *hTotLB, long int startValidity, long int endValidity, int runNumber) {
//void createCCDB(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB, THnSparse *hTotLB, long int startValidity, long int endValidity, int runNumber, string period) {
void createCCDB(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB, THnSparse *hTotLB, long int startValidity, long int endValidity, int firstRun, int lastRun, string period, bool merge, int trackGoal) {

    bool debug = false;

    //Variables used in the loop
    int LB936 = 0; //LB 1-> 936
    int plane = 0; //plane

    int counter = 0; //temporary variable to hold number of LB

    /*TH1F *hTestTotCounts = new TH1F(("hTestTotCounts_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBP = new TH1F(("hTestBP_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestNBP = new TH1F(("hTestNBP_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBoth = new TH1F(("hTestBoth_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);

    TH1F *hTestTotCountsAfter = new TH1F(("hTestTotCountsAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBPAfter = new TH1F(("hTestBPAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestNBPAfter = new TH1F(("hTestNBPAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);
    TH1F *hTestBothAfter = new TH1F(("hTestBothAfter_"+to_string(runNumber)).c_str(),("hTestTotCounts_"+to_string(runNumber)).c_str(),936,0.5,936.5);*/

    TH1F *hTestTotCounts = new TH1F(("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBP = new TH1F(("hTestBP_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestNBP = new TH1F(("hTestNBP_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBoth = new TH1F(("hTestBoth_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);

    TH1F *hTestTotCountsAfter = new TH1F(("hTestTotCountsAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBPAfter = new TH1F(("hTestBPAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestNBPAfter = new TH1F(("hTestNBPAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);
    TH1F *hTestBothAfter = new TH1F(("hTestBothAfter_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),("hTestTotCounts_"+to_string(firstRun)+"_"+to_string(lastRun)).c_str(),936,0.5,936.5);

    std::vector<o2::mid::ChEffCounter> counterVector; //Vector of efficiency counter structs
    struct o2::mid::ChEffCounter entry; //Struct for efficiency counters
    o2::mid::ChamberEfficiency effMap; //Chamber efficiency map

    //Get projection of LB counts before pt cut
    TH1D* totLBCountsProj = hTotLB->Projection(0); 
    TH1D* BothPlanesLBCountsProj = hFiredBothPlanesLB->Projection(0);
    TH1D* BPLBCountsProj = hFiredBPLB->Projection(0);
    TH1D* NBPLBCountsProj = hFiredNBPLB->Projection(0); 

    //Before cuts
    cout << "Runs between: " << firstRun << " and: " << lastRun << " before cuts" << endl;
    cout << "Total counts for LB "  << totLBCountsProj->GetEntries() << endl;
    cout << "Total counts BP LB "  << BPLBCountsProj->GetEntries() << endl;
    cout << "Total counts NBP LB "  << NBPLBCountsProj->GetEntries() << endl;
    cout << "Total counts Both LB "  << BothPlanesLBCountsProj->GetEntries() << endl; 

    //vector of struct, push back the struct populated with counts
    for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
       for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
           for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                hTestTotCounts->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LB936));
                hTestBP->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LB936));
                hTestNBP->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LB936));
                hTestBoth->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LB936));

                counter++;
            }
        }
    }

    counter = 0;

    /*delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;*/

    //Only accept muon track if pt > 2 GeV/c
    hTotLB->GetAxis(1)->SetRange(15,150);
    hFiredBothPlanesLB->GetAxis(1)->SetRange(15,150);
    hFiredBPLB->GetAxis(1)->SetRange(15,150);
    hFiredNBPLB->GetAxis(1)->SetRange(15,150);

    //Get projection of LB counts after pt cut
    totLBCountsProj = hTotLB->Projection(0); 
    BothPlanesLBCountsProj = hFiredBothPlanesLB->Projection(0);
    BPLBCountsProj = hFiredBPLB->Projection(0);
    NBPLBCountsProj = hFiredNBPLB->Projection(0); 

    //After cuts
    cout << "Runs between: " << firstRun << " and: " << lastRun << " after cuts" << endl;
    cout << "Total counts for LB "  << totLBCountsProj->GetEntries() << endl;
    cout << "Total counts BP LB "  << BPLBCountsProj->GetEntries() << endl;
    cout << "Total counts NBP LB "  << NBPLBCountsProj->GetEntries() << endl;
    cout << "Total counts Both LB "  << BothPlanesLBCountsProj->GetEntries() << endl;

    //vector of struct, push back the struct populated with counts
    for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
        for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
            for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                
                //Debug printout
                if (debug) {
                    cout << "det ID " << ide << " col " << icol << " line " << iline << " LB " << mapping.getBoardId(iline,icol,ide);
                    cout << "LB: " << mapping.getBoardId(iline,icol,ide) << "\t unique FEEID:" << o2::mid::detparams::makeUniqueFEEId(ide, icol, iline) << "\t";
                }
                
                plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                entry.deId =  uint8_t(ide);
                entry.columnId = uint8_t(icol);
                entry.lineId = uint8_t(iline);

                //entry.counts[0] = hFiredBPLB->GetBinContent(LB936); //BP
                //entry.counts[1] = hFiredNBPLB->GetBinContent(LB936); //NBP
                //entry.counts[2] = hFiredBothPlanesLB->GetBinContent(LB936); //Both
                //entry.counts[3] = hTotLB->GetBinContent(LB936); //Total

                entry.counts[0] = BPLBCountsProj->GetBinContent(LB936); //BP
                entry.counts[1] = NBPLBCountsProj->GetBinContent(LB936); //NBP
                entry.counts[2] = BothPlanesLBCountsProj->GetBinContent(LB936); //Both
                entry.counts[3] = totLBCountsProj->GetBinContent(LB936); //Total

                counterVector.push_back(entry);

                hTestTotCountsAfter->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LB936));
                hTestBPAfter->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LB936));
                hTestNBPAfter->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LB936));
                hTestBothAfter->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LB936));

                counter++;
                
                //Debug printout
                if (debug) {
                    cout << LB936 << "\t both: " << hFiredBothPlanesLB->GetBinContent(LB936) << "\t tot: " << hTotLB->GetBinContent(LB936) << endl;
                } 
            }
        }
    }

    //Delete objects here to avoid saturating RAM
    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    /*delete hTestTotCounts;
    delete hTestBP;
    delete hTestNBP;
    delete hTestBoth;

    delete hTestTotCountsAfter;
    delete hTestBPAfter;
    delete hTestNBPAfter;
    delete hTestBothAfter;*/

    delete hTotLB;
    delete hFiredBothPlanesLB;
    delete hFiredBPLB;
    delete hFiredNBPLB;

    //TFile *fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/LHC22_pass7_skimmed/ccdb/o2-mid-ChEffCounter_"+to_string(runNumber)+".root").c_str(),"RECREATE");
    cout << "Opening file for CCDB" << endl;
    //TFile *fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/o2-mid-ChEffCounter_"+to_string(runNumber)+".root").c_str(),"RECREATE");
    //TFile *fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/o2-mid-ChEffCounter_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    TFile *fOut;
    
    if (merge) {
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str());
        }

        fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/" +to_string(trackGoal) + "/o2-mid-ChEffCounter_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    }
    else {
        fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/o2-mid-ChEffCounter_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    }
    fOut->WriteObjectAny(&counterVector, "std::vector<o2::mid::ChEffCounter>","ccdb-object");
    fOut->Close();

    cout << "Opening file for test purposes" << endl;
    //TFile *fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/details_run_" + to_string(runNumber) + ".root").c_str(),"RECREATE");
    TFile *fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/details_run_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    fTest->cd();
    hTestTotCounts->Write("tot");
    hTestBP->Write("BP");
    hTestNBP->Write("NBP");
    hTestBoth->Write("Both");
    hTestTotCountsAfter->Write("totAfter");
    hTestBPAfter->Write("BPAfter");
    hTestNBPAfter->Write("NBPAfter");
    hTestBothAfter->Write("BothAfter");
    fTest->Close();

    cout << "Closing test file" << endl;
    
} //end of createCCDB


//---------------//
//               //
// Main function //
//               //
//---------------//

void produceObjects() { //Main function

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    int trackGoal = 7e+6;
    int tracks = 0, cumulativeTracks = 0;

    //To keep track of the number of merges and the number of times the track number was below the goal (to calculate average run number)
    int mergeCounter = 0, avgCalculation = 0;
    //To sum runs to calculate "average" run number
    float mergeRun = 0;

    //Each element of this vector is the average run number for which the merge is done
    //e.g. the merge is done between runs 340, 342 and 345 -> (340+342+345)/3=342
    vector<float> averageRun;

    //First and last run numbers in each merge
    long int first, last;
    //Start time of first and last runs in the merge
    long int startTime, endTime;

    bool open = false; //used to keep track if the .dat file (updated to contain the list of .root files to be merged) has been opened once
    bool assigned = false; //used to keep track if the first run of the merge has been assigned or not (to find proper timestamp range)

    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    bool createByRun = false;
    bool createByRunPbPb = false;

    //Merge more files until desired number of tracks is reached (if merge == true)
    bool merge = true;

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    //string period = "LHC22_pass7_skimmed";
    //string period = "LHC23_pass4_skimmed";
    string period = "LHC23_PbPb_pass4";
    //string period = "LHC24_pass1_skimmed";
    /string period = "LHC25ad_pass2"; //pO
    //string period = "LHC25ae_pass2"; //O-O
    //string period = "LHC25af_pass2"; //Ne-Ne
    
    string globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/";

    string fileName = "AnalysisResults.root";

    //Check period name and decide if it's Pb-Pb or pp
    if (period == "LHC23_pass4_skimmed" || period == "LHC24_pass1_skimmed" || period == "LHC22_pass7_skimmed") { //pp data
        createByRun = true;
        createByRunPbPb = false;
    }

    else if (period == "LHC23_PbPb_pass3_I-A11" || period == "LHC23_PbPb_pass4") { //Pb-Pb data
        createByRun = false;
        createByRunPbPb = true;
    }

    else {
        cout << "Wrong period name inserted, please check it" << endl;
        return;
    }

    //Load hadd.C macro to merge the root files from different runs
    //gROOT->ProcessLine(".L /home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/hadd.C");
    gROOT->ProcessLine(".L /home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/haddTHNsparse.C");

    //Output file for the merge of root files if the number of tracks reaches the desired goal
    ofstream hMergeRuns;
    
    string runNumbers = globalPath+"run_IR_Bfield.txt";
    //string runDates = globalPath+"run_dates.txt";

    //Open txt file of runs + start/end of runs
    ifstream hRun;
    hRun.open(runNumbers.c_str());

    //Get start and end of each run
    bool isIn;
    float IR; //IR is also provided as information in this file (for pp it is not used)
    string bField; //plus/minus B-field info is given in this file.
    int run;
    long int start, end;
    vector<long int> vStart, vEnd;
    vector<int> vRun;

    while (hRun >> isIn >> run >> IR >> bField >> start >> end){
        if (isIn) {
            vRun.push_back(run);
            vStart.push_back(start);
            vEnd.push_back(end);
        }
    }

    cout << vRun.size() << "\t" << vStart.size() << "\t" << vEnd.size() << endl;

    //Loop on all runs
    for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
        //Enter the folder
        string runFolder = globalPath + "runs/" +to_string(vRun.at(iRun));
        string runPath = globalPath + "runs/";
        
        cout << runFolder << endl;

        //run file name = path of the folder + run number (runFolder) + fileName
        string runFileName = runFolder+"/"+fileName;
        TFile *fRun = new TFile(runFileName.c_str(),"READ");

        TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

        TH1F *hFiredBothPlanesLB = (TH1F*)d->Get("nFiredBothperBoard");
        TH1F *hFiredBPLB = (TH1F*)d->Get("nFiredBPperBoard");
        TH1F *hFiredNBPLB = (TH1F*)d->Get("nFiredNBPperBoard");
        TH1F *hTotLB = (TH1F*)d->Get("nTotperBoard");
        //Track type
        TH1F *hTrackType = (TH1F*)d->Get("hTrackType");

        THnSparse *hSparseCentFiredTotPerLB = (THnSparse*)d->Get("hSparseCentFiredTotperBoard");
        THnSparse *hSparseCentFiredBothPerLB = (THnSparse*)d->Get("hSparseCentFiredBothperBoard");
        THnSparse *hSparseCentFiredBPPerLB = (THnSparse*)d->Get("hSparseCentFiredBPperBoard");
        THnSparse *hSparseCentFiredNBPPerLB = (THnSparse*)d->Get("hSparseCentFiredNBPperBoard");

        if (merge) { //Merge is enabled
            //continue;
            
            cout << "Merging of tracks is active" << endl;

            //Old 
            //tracks = hTotLB->GetEntries();
            //New 
            tracks = hTrackType->GetEntries();
            cumulativeTracks += tracks;
            
            cout << "Run number " << vRun.at(iRun) << " tot tracks in all LB " << tracks << " cumulative " << cumulativeTracks << endl;

            delete hSparseCentFiredTotPerLB;
            delete hSparseCentFiredBothPerLB;
            delete hSparseCentFiredBPPerLB;
            delete hSparseCentFiredNBPPerLB;

            delete hFiredBothPlanesLB,
            delete hFiredBPLB;
            delete hFiredNBPLB;
            delete hTotLB;
            delete hTrackType;

            mergeRun+=vRun.at(iRun); //sum run number to calculate "average" run number
            avgCalculation++; //increase by one for average calculations

            if (!assigned) {
                //first = vStart.at(iRun);
                startTime = vStart.at(iRun);
                first = vRun.at(iRun);
                assigned = true;
            }

            if (cumulativeTracks < trackGoal && vRun.at(iRun) != vRun.back()) { //If total track number is below the target -> Fill the file with the path of each AnalysisResults.root from each run
                //Open the output file only if it has not been opened (i.e. open == false)
                //meaning that it's the first run to be analyzed
                cout << "Cumulative tracks < track goal and not last run" << endl;
                if (!open) {
                    //hMergeRuns.open("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/merged_files/runs.dat");
                    hMergeRuns.open(runPath+"runs.dat");
                    open = true;
                    cout << "Opening file for merge" << endl;
                }
                hMergeRuns << runFileName << "\n"; //Write to output file
            }

            //Total track number greater than goal 
            // write to file the last run which got above the track goal
            // increase the number of merges
            // set the cumulative number of tracks to 0
            // close the output .dat file and set the open variable to false
            // merge the root files calling hadd.C
            // empty the runs.dat file (close the file and set the bool "open to false")
            // save the merged object inside a folder
            // call efficiency calculator function
            else if (cumulativeTracks >= trackGoal) {
                cout << "Track goal reached!" << endl;
                if (!open) {
                    //hMergeRuns.open("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/merged_files/runs.dat");
                    hMergeRuns.open(runPath+"runs.dat");
                    open = true;
                }
                hMergeRuns << runFileName << "\n";
                mergeCounter++;
                cumulativeTracks = 0;
                hMergeRuns.close();
                open = false;
                assigned = false;
                endTime = vEnd.at(iRun);
                last = vRun.at(iRun);
                //Test - create a sub-folder inside the merged_files directory
                gSystem->mkdir((runPath+"mergedRuns"+to_string(mergeCounter)).c_str());
                //Execute the hadd code (already loaded before the loop)
                //string mergeFilesForHadd = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root"+'"';
                string mergeFilesForHadd = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults_from_" + to_string(first) + "_to_" + to_string(last) + ".root"+'"';
                string mergeRunListForHadd = '"'+runPath+"runs.dat"+'"';
                //string mergeFiles = runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root";
                string mergeFiles = runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults_from_" + to_string(first) + "_to_" + to_string(last) + ".root";
                cout << "mergeFiles for eff calculations: " << mergeFiles << endl;
                cout << "mergeFilesForHadd: " << mergeFilesForHadd << endl;
                //gROOT->ProcessLine(Form("hadd(%s,%s)",mergeFilesForHadd.c_str(),mergeRunListForHadd.c_str()));
                if (last != first) {
                    cout << "Processing hadd" << endl;
                    gROOT->ProcessLine(Form("haddTHNsparse(%s,%s)",mergeFilesForHadd.c_str(),mergeRunListForHadd.c_str()));
                    averageRun.push_back(mergeRun/avgCalculation);
                }

                else {
                    cout << "No need to process hadd since track goal is reached with a single file, copying it directly" << endl;
                    averageRun.push_back(mergeRun/avgCalculation);

                    string cmd = "cp \"" + runFileName + "\" \"" + mergeFiles + "\"";
                    cout << "Executing: " << cmd << endl;
                    
                    int ret = gSystem->Exec(cmd.c_str());

                    if (ret == 0) {
                        std::cout << "File copied successfully!" << std::endl;
                    } else {
                        std::cerr << "Error while copying file!" << std::endl;
                    }
                }
                
                cout << "avg run " << mergeRun/avgCalculation << "\t mergeRun " << mergeRun << "\t avgCalculation " << avgCalculation << endl;

                TFile *fMerged = new TFile(mergeFiles.c_str(),"READ");

                TDirectoryFile *dMerged = (TDirectoryFile*)fMerged->Get("mid-efficiency");

                //TH1F *hFiredBothPlanesLBmerged = (TH1F*)dMerged->Get("nFiredBothperBoard");
                //TH1F *hFiredBPLBmerged = (TH1F*)dMerged->Get("nFiredBPperBoard");
                //TH1F *hFiredNBPLBmerged = (TH1F*)dMerged->Get("nFiredNBPperBoard");
                //TH1F *hTotLBmerged = (TH1F*)dMerged->Get("nTotperBoard");

                THnSparse *hSparseCentFiredTotPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredTotperBoard");
                THnSparse *hSparseCentFiredBothPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredBothperBoard");
                THnSparse *hSparseCentFiredBPPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredBPperBoard");
                THnSparse *hSparseCentFiredNBPPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredNBPperBoard");

                if (createByRun) {
                    //createCCDB(hFiredBPLB,hFiredNBPLB,hFiredBothPlanesLB,hTotLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun));
                    createCCDB(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,startTime,endTime,first,last,period,merge,trackGoal);
                }
            
                //Create .root objects to upload them to CCDB for PbPb
                if (createByRunPbPb) {            
                    createCCDBPbPb(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,startTime,endTime,first,last,period,merge,trackGoal);
                }
                //calculateEfficiencyByRun(mergeFiles,vEffBPLBmerged,vEffNBPLBmerged,vEffBothLBmerged,vEffBPLB_merged,vEffNBPLB_merged,vEffBothLB_merged,
                //vErrEffBPLBmerged, vErrEffNBPLBmerged,vErrEffBothLBmerged,vErrEffBPLB_merged,vErrEffNBPLB_merged,vErrEffBothLB_merged,first,last);
                mergeRun = 0;
                avgCalculation = 0;

                fMerged->Close();
                delete fMerged;
                cout << "Closed merged file" << endl;

                delete hSparseCentFiredTotPerLBmerged;
                delete hSparseCentFiredBothPerLBmerged; 
                delete hSparseCentFiredBPPerLBmerged;
                delete hSparseCentFiredNBPPerLBmerged;

                //gROOT->GetListOfFiles()->Delete();        // forget open files
                //gROOT->GetListOfCleanups()->Delete(); 
            }

            //Edge case, it's the last run of the list and the number of cumulative tracks has not yet reached the goal
            // we still merge those files and that's it so we copied the same code as the merges earlier on
            else if (vRun.at(iRun) == vRun.back() && cumulativeTracks < trackGoal) {
                cout << "Last run of the list and track goal not yet reached, merging anyway" << endl;
                hMergeRuns << runFileName << "\n";
                mergeCounter++;
                //mergeRun+=vRun.at(iRun);
                //avgCalculation++;
                cumulativeTracks = 0;
                hMergeRuns.close();
                open = false; //no need to set it to false but for redudancy we do it
                assigned = false;
                endTime = vEnd.at(iRun);
                last = vRun.at(iRun);
                //Test - create a sub-folder inside the merged_files directory
                gSystem->mkdir((runPath+"mergedRuns"+to_string(mergeCounter)).c_str());
                //Execute the hadd code (already loaded before the loop)
                //string mergeFilesForHadd = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root"+'"';
                string mergeFilesForHadd = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults_from_" + to_string(first) + "_to_" + to_string(last) + ".root"+'"';
                string mergeRunListForHadd = '"'+runPath+"runs.dat"+'"';
                //string mergeFiles = runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root";
                string mergeFiles = runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults_from_" + to_string(first) + "_to_" + to_string(last) + ".root";
                cout << "mergeFilesForHadd: " << mergeFilesForHadd << endl;
                //gROOT->ProcessLine(Form("hadd(%s,%s)",mergeFilesForHadd.c_str(),mergeRunListForHadd.c_str()));
                if (last != first) {
                    cout << "Processing hadd" << endl;
                    gROOT->ProcessLine(Form("haddTHNsparse(%s,%s)",mergeFilesForHadd.c_str(),mergeRunListForHadd.c_str()));
                    averageRun.push_back(mergeRun/avgCalculation);
                }

                else {
                    cout << "No need to process hadd since track goal is reached with a single file, copying it directly" << endl;
                    averageRun.push_back(mergeRun/avgCalculation);

                    string cmd = "cp \"" + runFileName + "\" \"" + mergeFiles + "\"";
                    cout << "Executing: " << cmd << endl;
                    
                    int ret = gSystem->Exec(cmd.c_str());

                    if (ret == 0) {
                        std::cout << "File copied successfully!" << std::endl;
                    } else {
                        std::cerr << "Error while copying file!" << std::endl;
                    }
                }
                
                //gROOT->ProcessLine(Form("haddTHNsparse(%s,%s)",mergeFilesForHadd.c_str(),mergeRunListForHadd.c_str()));
                //averageRun.push_back(mergeRun/avgCalculation);
                cout << "avg run " << mergeRun/avgCalculation << "\t mergeRun " << mergeRun << "\t avgCalculation " << avgCalculation << endl;

                TFile *fMerged = new TFile(mergeFiles.c_str(),"READ");

                TDirectoryFile *dMerged = (TDirectoryFile*)fMerged->Get("mid-efficiency");

                //TH1F *hFiredBothPlanesLBmerged = (TH1F*)dMerged->Get("nFiredBothperBoard");
                //TH1F *hFiredBPLBmerged = (TH1F*)dMerged->Get("nFiredBPperBoard");
                //TH1F *hFiredNBPLBmerged = (TH1F*)dMerged->Get("nFiredNBPperBoard");
                //TH1F *hTotLBmerged = (TH1F*)dMerged->Get("nTotperBoard");

                THnSparse *hSparseCentFiredTotPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredTotperBoard");
                THnSparse *hSparseCentFiredBothPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredBothperBoard");
                THnSparse *hSparseCentFiredBPPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredBPperBoard");
                THnSparse *hSparseCentFiredNBPPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredNBPperBoard");

                if (createByRun) {
                    //createCCDB(hFiredBPLB,hFiredNBPLB,hFiredBothPlanesLB,hTotLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun));
                    createCCDB(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,startTime,endTime,first,last,period,merge,trackGoal);
                }
            
                //Create .root objects to upload them to CCDB for PbPb
                if (createByRunPbPb) {            
                    createCCDBPbPb(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,startTime,endTime,first,last,period,merge,trackGoal);
                }
                //calculateEfficiencyByRun(mergeFiles,vEffBPLBmerged,vEffNBPLBmerged,vEffBothLBmerged,vEffBPLB_merged,vEffNBPLB_merged,vEffBothLB_merged,
                //vErrEffBPLBmerged, vErrEffNBPLBmerged, vErrEffBothLBmerged,vErrEffBPLB_merged, vErrEffNBPLB_merged, vErrEffBothLB_merged,first,last);
                mergeRun = 0;
                avgCalculation = 0;

                delete hSparseCentFiredTotPerLBmerged;
                delete hSparseCentFiredBothPerLBmerged; 
                delete hSparseCentFiredBPPerLBmerged;
                delete hSparseCentFiredNBPPerLBmerged;

                fMerged->Close();
                delete fMerged;
            }
            
        } //End of if on merge
        
        else { //Merge is disabled
            cout << "Merging of tracks is not active" << endl;
            //Create .root objects to upload them to CCDB
            if (createByRun) {
                //createCCDB(hFiredBPLB,hFiredNBPLB,hFiredBothPlanesLB,hTotLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun));
                createCCDB(hSparseCentFiredBPPerLB,hSparseCentFiredNBPPerLB,hSparseCentFiredBothPerLB,hSparseCentFiredTotPerLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun),vRun.at(iRun),period,merge,trackGoal);

            }
            
            //Create .root objects to upload them to CCDB for PbPb
            if (createByRunPbPb) {            
                createCCDBPbPb(hSparseCentFiredBPPerLB,hSparseCentFiredNBPPerLB,hSparseCentFiredBothPerLB,hSparseCentFiredTotPerLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun),vRun.at(iRun),period,merge,trackGoal);
            }
        }
            
        //Calculate efficiency per run -> not needed here but just to print it out for debug purposes
        /*for (int i = 1; i <= nBinsBoard; i++) {

            //if (hTotLB->GetBinContent(i) == 0) {
            //    cout << "Run: " << vRun.at(iRun) << " LB " <<  i << "\t" << hTotLB->GetBinContent(i) << endl;
            //}

            if (hTotLB->GetBinContent(i) != 0) {

                effBothLB = (hFiredBothPlanesLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effBPLB = (hFiredBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effNBPLB = (hFiredNBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;

                errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLB->GetBinContent(i));
                errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLB->GetBinContent(i));
                errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLB->GetBinContent(i));
                //Debug printout
                //cout << "Run: " << vRun.at(iRun) << " LB " <<  i << "\t" << effBothLB << "+-" << errEffBothLB << "\t" << effBPLB << "+-" << errEffBPLB << "\t" << effNBPLB << "+-" << errEffNBPLB << endl;
            }

            //I noticed that sometimes in a given run there is a LB with 0 total counts
            else {
               cout << "For run number " << vRun.at(iRun) << " the LB " << i << " has zero entries" << endl;
            }
        }*/

        fRun->Close();
        delete fRun;
    } //End of loop on all runs
    
    hRun.close();
    //hDate.close();
}