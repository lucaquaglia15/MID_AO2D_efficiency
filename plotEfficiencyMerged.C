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
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

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

//Max and min pt values in GeV/c
float minPt = 0.;
float maxPt = 20.;
int binsPt = 150;
float ptStep = (maxPt-minPt)/binsPt;
int ptMinCut = 1.8/ptStep;

int LBnumber = 35;

//o2::ccdb::CcdbApi api; //CCDB API as global object
o2::mid::Mapping mapping; //MID mapping object to construct ccdb object

//Function to upload to ccdb for PbPb data -> only get events with centrality > 50%
//void plotMergedPbPb(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB, THnSparse *hTotLB, int firstRun, int lastRun, string period, bool merge, int trackGoal) {
tuple<float, float, float, float, float, float, float, float, float, float, float, float> plotMergedPbPb(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB, 
    THnSparse *hTotLB, int firstRun, int lastRun, string period, bool merge, int trackGoal) {
    
    bool debug = false;

    //Variables used in the loop
    int LB936 = 0; //LB 1-> 936
    int plane = 0; //plane
    int counter = 0; //temporary variable to hold number of LB

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

    float effBP_before, effNBP_before, effBoth_before;
    float errEffBP_before, errEffNBP_before, errEffBoth_before;
    
    //vector of struct, push back the struct populated with counts
    //for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
        //for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
            //for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                
                //plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                //LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                float totCounts = totLBCountsProj->GetBinContent(LBnumber); //tot
                float BPcounts = BPLBCountsProj->GetBinContent(LBnumber); //BP
                float NBPcounts = NBPLBCountsProj->GetBinContent(LBnumber); //NBP
                float Bothcounts = BothPlanesLBCountsProj->GetBinContent(LBnumber); //Both

                effBP_before = (BPcounts/totCounts)*100;
                effNBP_before = (NBPcounts/totCounts)*100;
                effBoth_before = (Bothcounts/totCounts)*100;
                errEffBP_before = TMath::Sqrt(effBP_before*(100-effBP_before)/totCounts);
                errEffNBP_before = TMath::Sqrt(effNBP_before*(100-effNBP_before)/totCounts);
                errEffBoth_before = TMath::Sqrt(effBoth_before*(100-effBoth_before)/totCounts);

                hTestTotCounts->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LBnumber));
                hTestBP->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LBnumber));
                hTestNBP->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LBnumber));
                hTestBoth->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LBnumber));

                //counter++;
            //}
        //}
    //}

    counter = 0;

    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    //Only accept muon track if pt > 2 GeV/c
    hTotLB->GetAxis(2)->SetRange(ptMinCut,150);
    hFiredBothPlanesLB->GetAxis(2)->SetRange(ptMinCut,150);
    hFiredBPLB->GetAxis(2)->SetRange(ptMinCut,150);
    hFiredNBPLB->GetAxis(2)->SetRange(ptMinCut,150);

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

    float effBP_after, effNBP_after, effBoth_after;
    float errEffBP_after, errEffNBP_after, errEffBoth_after;

    //vector of struct, push back the struct populated with counts
    //for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
    //   for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
    //       for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                
                //Debug printout
                /*if (debug) {
                    cout << "det ID " << ide << " col " << icol << " line " << iline << " LB " << mapping.getBoardId(iline,icol,ide) << endl;
                    cout << "LB: " << mapping.getBoardId(iline,icol,ide) << "\t unique FEEID:" << o2::mid::detparams::makeUniqueFEEId(ide, icol, iline) << "\t";
                }*/
                
                //plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                //LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                totCounts = totLBCountsProj->GetBinContent(LBnumber); //tot
                BPcounts = BPLBCountsProj->GetBinContent(LBnumber); //BP
                NBPcounts = NBPLBCountsProj->GetBinContent(LBnumber); //NBP
                Bothcounts = BothPlanesLBCountsProj->GetBinContent(LBnumber); //Both

                effBP_after = (BPcounts/totCounts)*100;
                effNBP_after = (NBPcounts/totCounts)*100;
                effBoth_after = (Bothcounts/totCounts)*100;
                errEffBP_after = TMath::Sqrt(effBP_after*(100-effBP_after)/totCounts);
                errEffNBP_after = TMath::Sqrt(effNBP_after*(100-effNBP_after)/totCounts);
                errEffBoth_after = TMath::Sqrt(effBoth_after*(100-effBoth_after)/totCounts);

                hTestTotCountsAfter->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LBnumber));
                hTestBPAfter->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LBnumber));
                hTestNBPAfter->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LBnumber));
                hTestBothAfter->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LBnumber));

                //counter++;

                //Debug printout
                //if (debug) {
                //    cout << LB936 << "\t both: " << BothPlanesLBCountsProj->GetBinContent(LB936) << "\t tot: " << totLBCountsProj->GetBinContent(LB936) << endl;
                //} 
        //    }
      //  }
    //}

    //Delete objects here to avoid saturating RAM
    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    TFile *fTest;

    if (merge) {
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str())) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str());
        }
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str());
        }

        fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/" + to_string(trackGoal) + "/details_run_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");

    }
    else {
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(0)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(0)).c_str());
        }
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str());
        }

        fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/" + to_string(trackGoal) + "details_run_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    }

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

    return  make_tuple(effBP_before, errEffBP_before, effNBP_before, errEffNBP_before, effBoth_before, errEffBoth_before,
                    effBP_after, errEffBP_after, effNBP_after, errEffNBP_after, effBoth_after, errEffBoth_after);
} //end of createCCDB for PbPb

//////////

//Function to upload to ccdb
tuple<float, float, float, float, float, float, float, float, float, float, float, float> plotMerged(THnSparse *hFiredBPLB, THnSparse *hFiredNBPLB, THnSparse *hFiredBothPlanesLB,
    THnSparse *hTotLB, int firstRun, int lastRun, string period, bool merge, int trackGoal) {

    bool debug = false;

    //Variables used in the loop
    int LB936 = 0; //LB 1-> 936
    int plane = 0; //plane

    int counter = 0; //temporary variable to hold number of LB

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

    float effBP_before, effNBP_before, effBoth_before;
    float errEffBP_before, errEffNBP_before, errEffBoth_before;

    //vector of struct, push back the struct populated with counts
    //for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
    //   for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
    //       for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
    //            plane = o2::mid::detparams::getChamber(ide); //Get detection plane

    //            LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                float totCounts = totLBCountsProj->GetBinContent(LBnumber); //tot
                float BPcounts = BPLBCountsProj->GetBinContent(LBnumber); //BP
                float NBPcounts = NBPLBCountsProj->GetBinContent(LBnumber); //NBP
                float Bothcounts = BothPlanesLBCountsProj->GetBinContent(LBnumber); //Both

                effBP_before = (BPcounts/totCounts)*100;
                effNBP_before = (NBPcounts/totCounts)*100;
                effBoth_before = (Bothcounts/totCounts)*100;
                errEffBP_before = TMath::Sqrt(effBP_before*(100-effBP_before)/totCounts);
                errEffNBP_before = TMath::Sqrt(effNBP_before*(100-effNBP_before)/totCounts);
                errEffBoth_before = TMath::Sqrt(effBoth_before*(100-effBoth_before)/totCounts);

                hTestTotCounts->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LBnumber));
                hTestBP->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LBnumber));
                hTestNBP->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LBnumber));
                hTestBoth->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LBnumber));

                //counter++;
    //        }
    //    }
    //}

    counter = 0;

    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    //Only accept muon track if pt > 2 GeV/c
    hTotLB->GetAxis(1)->SetRange(ptMinCut,150);
    hFiredBothPlanesLB->GetAxis(1)->SetRange(ptMinCut,150);
    hFiredBPLB->GetAxis(1)->SetRange(ptMinCut,150);
    hFiredNBPLB->GetAxis(1)->SetRange(ptMinCut,150);

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

    float effBP_after, effNBP_after, effBoth_after;
    float errEffBP_after, errEffNBP_after, errEffBoth_after;

    //vector of struct, push back the struct populated with counts
    //for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
    //    for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
    //        for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                
                //Debug printout
                /*if (debug) {
                    cout << "det ID " << ide << " col " << icol << " line " << iline << " LB " << mapping.getBoardId(iline,icol,ide);
                    cout << "LB: " << mapping.getBoardId(iline,icol,ide) << "\t unique FEEID:" << o2::mid::detparams::makeUniqueFEEId(ide, icol, iline) << "\t";
                }*/
                
                //plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                //LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                totCounts = totLBCountsProj->GetBinContent(LBnumber); //tot
                BPcounts = BPLBCountsProj->GetBinContent(LBnumber); //BP
                NBPcounts = NBPLBCountsProj->GetBinContent(LBnumber); //NBP
                Bothcounts = BothPlanesLBCountsProj->GetBinContent(LBnumber); //Both

                effBP_after = (BPcounts/totCounts)*100;
                effNBP_after = (NBPcounts/totCounts)*100;
                effBoth_after = (Bothcounts/totCounts)*100;
                errEffBP_after = TMath::Sqrt(effBP_after*(100-effBP_after)/totCounts);
                errEffNBP_after = TMath::Sqrt(effNBP_after*(100-effNBP_after)/totCounts);
                errEffBoth_after = TMath::Sqrt(effBoth_after*(100-effBoth_after)/totCounts);

                hTestTotCountsAfter->SetBinContent(counter+1,totLBCountsProj->GetBinContent(LBnumber));
                hTestBPAfter->SetBinContent(counter+1,BPLBCountsProj->GetBinContent(LBnumber));
                hTestNBPAfter->SetBinContent(counter+1,NBPLBCountsProj->GetBinContent(LBnumber));
                hTestBothAfter->SetBinContent(counter+1,BothPlanesLBCountsProj->GetBinContent(LBnumber));

                //counter++;
                
                //Debug printout
                /*if (debug) {
                    cout << LB936 << "\t both: " << hFiredBothPlanesLB->GetBinContent(LB936) << "\t tot: " << hTotLB->GetBinContent(LB936) << endl;
                }*/
    //        }
    //    }
    //}

    //Delete objects here to avoid saturating RAM
    delete totLBCountsProj;
    delete BothPlanesLBCountsProj;
    delete BPLBCountsProj;
    delete NBPLBCountsProj;

    cout << "Opening file for CCDB" << endl;
    TFile *fTest;
    
    if (merge) {
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(trackGoal)).c_str());
        }

        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str());
        }

        fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/" +to_string(trackGoal) + "/details_run_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");
    }

    else {
        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(0)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/merged/"+to_string(0)).c_str());
        }

        if (!fs::exists(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str() )) {
            fs::create_directories(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/"+to_string(trackGoal)).c_str());
        }

        fTest = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/detailedOutput/" + to_string(trackGoal) + "/details_run_"+to_string(firstRun)+"_"+to_string(lastRun)+".root").c_str(),"RECREATE");

    }

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

    return  make_tuple(effBP_before, errEffBP_before, effNBP_before, errEffNBP_before, effBoth_before, errEffBoth_before,
                    effBP_after, errEffBP_after, effNBP_after, errEffNBP_after, effBoth_after, errEffBoth_after);

} //end of createCCDB


//---------------//
//               //
// Main function //
//               //
//---------------//

void plotEfficiencyMerged() { //Main function

    cout << "ptMinCut: " << ptMinCut << endl;

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    vector<double> v_effBP_before, v_effNBP_before, v_effBoth_before, v_errEffBP_before, v_errEffNBP_before, v_errEffBoth_before;
    vector<double> v_effBP_after, v_effNBP_after, v_effBoth_after, v_errEffBP_after, v_errEffNBP_after, v_errEffBoth_after;
    vector<double> vRun_average;
    vector<double> vRunLow, vRunHigh;

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
    int trackGoal = 1e+6;

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    //string period = "LHC22_pass7_skimmed";
    string period = "LHC23_pass4_skimmed";
    //string period = "LHC23_PbPb_pass4";
    //string period = "LHC24_PbPb_pass2";
    //string period = "LHC24_pass1_skimmed";
    //string period = "LHC25ad_pass2"; //pO
    //string period = "LHC25ae_pass2"; //O-O
    //string period = "LHC25af_pass2"; //Ne-Ne
    //string period = "LHC24_ppref_pass1"; //pp ref 2024

    string globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/";

    string fileName = "AnalysisResults.root";

    //Check if folder exists

    //Check period name and decide if it's Pb-Pb or pp
    if (period == "LHC23_pass4_skimmed" || period == "LHC24_pass1_skimmed" || period == "LHC22_pass7_skimmed" || period == "LHC24_ppref_pass1") { //pp data
        createByRun = true;
        createByRunPbPb = false;
    }

    else if (period == "LHC23_PbPb_pass3_I-A11" || period == "LHC23_PbPb_pass4" || period == "LHC25ad_pass2"
          || period == "LHC25ae_pass1" || period == "LHC25ae_pass2" || period == "LHC25af_pass2") { //Pb-Pb/p-O/O-O/Ne-Ne data
        createByRun = false;
        createByRunPbPb = true;
    }

    else {
        cout << "Wrong period name inserted, please check it" << endl;
        return;
    }

    //Load hadd.C macro to merge the root files from different runs
    gROOT->ProcessLine(".L /home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/haddTHNsparse.C");

    //Output file for the merge of root files if the number of tracks reaches the desired goal
    ofstream hMergeRuns;
    
    string runNumbers = globalPath+"run_IR_Bfield.txt";

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
                    cout << "Opening file for merge " << runPath << "runs.dat" << endl;
                    cout << "Run file name: " << runFileName << endl;
                }
                hMergeRuns << runFileName << "\n"; //Write to output file
                cout << "Written" << endl;
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

                THnSparse *hSparseCentFiredTotPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredTotperBoard");
                THnSparse *hSparseCentFiredBothPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredBothperBoard");
                THnSparse *hSparseCentFiredBPPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredBPperBoard");
                THnSparse *hSparseCentFiredNBPPerLBmerged = (THnSparse*)dMerged->Get("hSparseCentFiredNBPperBoard");

                float effBP_before, effNBP_before, effBoth_before, errEffBP_before, errEffNBP_before, errEffBoth_before;
                float effBP_after, effNBP_after, effBoth_after, errEffBP_after, errEffNBP_after, errEffBoth_after;

                tuple<float,float,float,float,float,float,float,float,float,float,float,float> efficiency;

                if (createByRun) {
                    //plotMerged(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,first,last,period,merge,trackGoal);
                    efficiency = plotMerged(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,first,last,period,merge,trackGoal);
                    
                    effBP_before = get<0>(efficiency);
                    errEffBP_before = get<1>(efficiency);
                    effNBP_before = get<2>(efficiency);
                    errEffNBP_before = get<3>(efficiency);
                    effBoth_before = get<4>(efficiency);
                    errEffBoth_before = get<5>(efficiency);
                    effBP_after = get<6>(efficiency);
                    errEffBP_after = get<7>(efficiency);
                    effNBP_after = get<8>(efficiency);
                    errEffNBP_after = get<9>(efficiency);
                    effBoth_after = get<10>(efficiency);
                    errEffBoth_after = get<11>(efficiency);

                    v_effBP_before.push_back(effBP_before);
                    v_effNBP_before.push_back(effNBP_before);
                    v_effBoth_before.push_back(effBoth_before);
                    v_errEffBP_before.push_back(errEffBP_before);
                    v_errEffNBP_before.push_back(errEffNBP_before);
                    v_errEffBoth_before.push_back(errEffBoth_before);

                    v_effBP_after.push_back(effBP_after);
                    v_effNBP_after.push_back(effNBP_after);
                    v_effBoth_after.push_back(effBoth_after);
                    v_errEffBP_after.push_back(errEffBP_after);
                    v_errEffNBP_after.push_back(errEffNBP_after);
                    v_errEffBoth_after.push_back(errEffBoth_after);
                    
                    vRun_average.push_back(mergeRun/avgCalculation);
                    vRunLow.push_back((mergeRun/avgCalculation) - first);
                    vRunHigh.push_back(last - (mergeRun/avgCalculation));
                    cout << "first " << first << " last " << last << " avg " << mergeRun/avgCalculation << endl;
                }
            
                //Create .root objects to upload them to CCDB for PbPb
                if (createByRunPbPb) {            
                    //plotMergedPbPb(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,first,last,period,merge,trackGoal);
                    efficiency = plotMergedPbPb(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,first,last,period,merge,trackGoal);
                    
                    effBP_before = get<0>(efficiency);
                    errEffBP_before = get<1>(efficiency);
                    effNBP_before = get<2>(efficiency);
                    errEffNBP_before = get<3>(efficiency);
                    effBoth_before = get<4>(efficiency);
                    errEffBoth_before = get<5>(efficiency);
                    effBP_after = get<6>(efficiency);
                    errEffBP_after = get<7>(efficiency);
                    effNBP_after = get<8>(efficiency);
                    errEffNBP_after = get<9>(efficiency);
                    effBoth_after = get<10>(efficiency);
                    errEffBoth_after = get<11>(efficiency);

                    v_effBP_before.push_back(effBP_before);
                    v_effNBP_before.push_back(effNBP_before);
                    v_effBoth_before.push_back(effBoth_before);
                    v_errEffBP_before.push_back(errEffBP_before);
                    v_errEffNBP_before.push_back(errEffNBP_before);
                    v_errEffBoth_before.push_back(errEffBoth_before);

                    v_effBP_after.push_back(effBP_after);
                    v_effNBP_after.push_back(effNBP_after);
                    v_effBoth_after.push_back(effBoth_after);
                    v_errEffBP_after.push_back(errEffBP_after);
                    v_errEffNBP_after.push_back(errEffNBP_after);
                    v_errEffBoth_after.push_back(errEffBoth_after);
                    
                    vRun_average.push_back(mergeRun/avgCalculation);
                    vRunLow.push_back((mergeRun/avgCalculation) - first);
                    vRunHigh.push_back(last - (mergeRun/avgCalculation));
                    cout << "first " << first << " last " << last << " avg " << mergeRun/avgCalculation << endl;
                }
                //calculateEfficiencyByRun(mergeFiles,vEffBPLBmerged,vEffNBPLBmerged,vEffBothLBmerged,vEffBPLB_merged,vEffNBPLB_merged,vEffBothLB_merged,
                //vErrEffBPLBmerged, vErrEffNBPLBmerged,vErrEffBothLBmerged,vErrEffBPLB_merged,vErrEffNBPLB_merged,vErrEffBothLB_merged,first,last);
                mergeRun = 0;
                avgCalculation = 0;

                /*delete hSparseCentFiredTotPerLBmerged;
                delete hSparseCentFiredBothPerLBmerged; 
                delete hSparseCentFiredBPPerLBmerged;
                delete hSparseCentFiredNBPPerLBmerged;*/

                fMerged->Close();
                delete fMerged;

                cout << "Closed merged file" << endl;

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

                float effBP_before, effNBP_before, effBoth_before, errEffBP_before, errEffNBP_before, errEffBoth_before;
                float effBP_after, effNBP_after, effBoth_after, errEffBP_after, errEffNBP_after, errEffBoth_after;

                tuple<float,float,float,float,float,float,float,float,float,float,float,float> efficiency;

                if (createByRun) {
                    //createCCDB(hFiredBPLB,hFiredNBPLB,hFiredBothPlanesLB,hTotLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun));
                    efficiency = plotMerged(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,first,last,period,merge,trackGoal);
                    
                    effBP_before = get<0>(efficiency);
                    errEffBP_before = get<1>(efficiency);
                    effNBP_before = get<2>(efficiency);
                    errEffNBP_before = get<3>(efficiency);
                    effBoth_before = get<4>(efficiency);
                    errEffBoth_before = get<5>(efficiency);
                    effBP_after = get<6>(efficiency);
                    errEffBP_after = get<7>(efficiency);
                    effNBP_after = get<8>(efficiency);
                    errEffNBP_after = get<9>(efficiency);
                    effBoth_after = get<10>(efficiency);
                    errEffBoth_after = get<11>(efficiency);

                    v_effBP_before.push_back(effBP_before);
                    v_effNBP_before.push_back(effNBP_before);
                    v_effBoth_before.push_back(effBoth_before);
                    v_errEffBP_before.push_back(errEffBP_before);
                    v_errEffNBP_before.push_back(errEffNBP_before);
                    v_errEffBoth_before.push_back(errEffBoth_before);

                    v_effBP_after.push_back(effBP_after);
                    v_effNBP_after.push_back(effNBP_after);
                    v_effBoth_after.push_back(effBoth_after);
                    v_errEffBP_after.push_back(errEffBP_after);
                    v_errEffNBP_after.push_back(errEffNBP_after);
                    v_errEffBoth_after.push_back(errEffBoth_after);
                    
                    vRun_average.push_back(mergeRun/avgCalculation);
                    vRunLow.push_back((mergeRun/avgCalculation) - first);
                    vRunHigh.push_back(last - (mergeRun/avgCalculation));
                    cout << "first " << first << " last " << last << " avg " << mergeRun/avgCalculation << endl;
                }
            
                //Create .root objects to upload them to CCDB for PbPb
                if (createByRunPbPb) {            
                    efficiency = plotMergedPbPb(hSparseCentFiredBPPerLBmerged,hSparseCentFiredNBPPerLBmerged,hSparseCentFiredBothPerLBmerged,hSparseCentFiredTotPerLBmerged,first,last,period,merge,trackGoal);
                    
                    effBP_before = get<0>(efficiency);
                    errEffBP_before = get<1>(efficiency);
                    effNBP_before = get<2>(efficiency);
                    errEffNBP_before = get<3>(efficiency);
                    effBoth_before = get<4>(efficiency);
                    errEffBoth_before = get<5>(efficiency);
                    effBP_after = get<6>(efficiency);
                    errEffBP_after = get<7>(efficiency);
                    effNBP_after = get<8>(efficiency);
                    errEffNBP_after = get<9>(efficiency);
                    effBoth_after = get<10>(efficiency);
                    errEffBoth_after = get<11>(efficiency);

                    v_effBP_before.push_back(effBP_before);
                    v_effNBP_before.push_back(effNBP_before);
                    v_effBoth_before.push_back(effBoth_before);
                    v_errEffBP_before.push_back(errEffBP_before);
                    v_errEffNBP_before.push_back(errEffNBP_before);
                    v_errEffBoth_before.push_back(errEffBoth_before);

                    v_effBP_after.push_back(effBP_after);
                    v_effNBP_after.push_back(effNBP_after);
                    v_effBoth_after.push_back(effBoth_after);
                    v_errEffBP_after.push_back(errEffBP_after);
                    v_errEffNBP_after.push_back(errEffNBP_after);
                    v_errEffBoth_after.push_back(errEffBoth_after);
                    
                    vRun_average.push_back(mergeRun/avgCalculation);
                    vRunLow.push_back((mergeRun/avgCalculation) - first);
                    vRunHigh.push_back(last - (mergeRun/avgCalculation));
                    cout << "first " << first << " last " << last << " avg " << mergeRun/avgCalculation << endl;
                }

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
            float effBP_before, effNBP_before, effBoth_before, errEffBP_before, errEffNBP_before, errEffBoth_before;
            float effBP_after, effNBP_after, effBoth_after, errEffBP_after, errEffNBP_after, errEffBoth_after;

            tuple<float,float,float,float,float,float,float,float,float,float,float,float> efficiency;
            //Create .root objects to upload them to CCDB
            if (createByRun) {
                //plotMerged(hSparseCentFiredBPPerLB,hSparseCentFiredNBPPerLB,hSparseCentFiredBothPerLB,hSparseCentFiredTotPerLB,vRun.at(iRun),vRun.at(iRun),period,merge,trackGoal);
                efficiency = plotMergedPbPb(hSparseCentFiredBPPerLB,hSparseCentFiredNBPPerLB,hSparseCentFiredBothPerLB,hSparseCentFiredTotPerLB,vRun.at(iRun),vRun.at(iRun),period,merge,trackGoal);
                    
                effBP_before = get<0>(efficiency);
                errEffBP_before = get<1>(efficiency);
                effNBP_before = get<2>(efficiency);
                errEffNBP_before = get<3>(efficiency);
                effBoth_before = get<4>(efficiency);
                errEffBoth_before = get<5>(efficiency);
                effBP_after = get<6>(efficiency);
                errEffBP_after = get<7>(efficiency);
                effNBP_after = get<8>(efficiency);
                errEffNBP_after = get<9>(efficiency);
                effBoth_after = get<10>(efficiency);
                errEffBoth_after = get<11>(efficiency);

                v_effBP_before.push_back(effBP_before);
                v_effNBP_before.push_back(effNBP_before);
                v_effBoth_before.push_back(effBoth_before);
                v_errEffBP_before.push_back(errEffBP_before);
                v_errEffNBP_before.push_back(errEffNBP_before);
                v_errEffBoth_before.push_back(errEffBoth_before);

                v_effBP_after.push_back(effBP_after);
                v_effNBP_after.push_back(effNBP_after);
                v_effBoth_after.push_back(effBoth_after);
                v_errEffBP_after.push_back(errEffBP_after);
                v_errEffNBP_after.push_back(errEffNBP_after);
                v_errEffBoth_after.push_back(errEffBoth_after);
                
                vRun_average.push_back(vRun.at(iRun));
                vRunLow.push_back(0);
                vRunHigh.push_back(0);
                cout << "first " << vRun.at(iRun) << " last " << vRun.at(iRun) << " avg " << vRun.at(iRun) << endl;
            
            }
            
            //Create .root objects to upload them to CCDB for PbPb
            if (createByRunPbPb) {            
                //plotMergedPbPb(hSparseCentFiredBPPerLB,hSparseCentFiredNBPPerLB,hSparseCentFiredBothPerLB,hSparseCentFiredTotPerLB,vRun.at(iRun),vRun.at(iRun),period,merge,trackGoal);
                efficiency = plotMergedPbPb(hSparseCentFiredBPPerLB,hSparseCentFiredNBPPerLB,hSparseCentFiredBothPerLB,hSparseCentFiredTotPerLB,vRun.at(iRun),vRun.at(iRun),period,merge,trackGoal);
                    
                    effBP_before = get<0>(efficiency);
                    errEffBP_before = get<1>(efficiency);
                    effNBP_before = get<2>(efficiency);
                    errEffNBP_before = get<3>(efficiency);
                    effBoth_before = get<4>(efficiency);
                    errEffBoth_before = get<5>(efficiency);
                    effBP_after = get<6>(efficiency);
                    errEffBP_after = get<7>(efficiency);
                    effNBP_after = get<8>(efficiency);
                    errEffNBP_after = get<9>(efficiency);
                    effBoth_after = get<10>(efficiency);
                    errEffBoth_after = get<11>(efficiency);

                    v_effBP_before.push_back(effBP_before);
                    v_effNBP_before.push_back(effNBP_before);
                    v_effBoth_before.push_back(effBoth_before);
                    v_errEffBP_before.push_back(errEffBP_before);
                    v_errEffNBP_before.push_back(errEffNBP_before);
                    v_errEffBoth_before.push_back(errEffBoth_before);

                    v_effBP_after.push_back(effBP_after);
                    v_effNBP_after.push_back(effNBP_after);
                    v_effBoth_after.push_back(effBoth_after);
                    v_errEffBP_after.push_back(errEffBP_after);
                    v_errEffNBP_after.push_back(errEffNBP_after);
                    v_errEffBoth_after.push_back(errEffBoth_after);
                    
                    vRun_average.push_back(vRun.at(iRun));
                    vRunLow.push_back(0);
                    vRunHigh.push_back(0);
                    cout << "first " << vRun.at(iRun) << " last " << vRun.at(iRun) << " avg " << vRun.at(iRun) << endl;
            }

            delete hSparseCentFiredBPPerLB;
            delete hSparseCentFiredNBPPerLB;
            delete hSparseCentFiredBothPerLB;
            delete hSparseCentFiredTotPerLB;
        }
            
        /*delete hSparseCentFiredBPPerLB;
        delete hSparseCentFiredNBPPerLB;
        delete hSparseCentFiredBothPerLB;
        delete hSparseCentFiredTotPerLB;*/

        fRun->Close();
        delete fRun;
    } //End of loop on all runs

    TGraphAsymmErrors *gBPEff_before = new TGraphAsymmErrors(vRun_average.size(),&vRun_average[0],&v_effBP_before[0],&vRunLow[0],&vRunHigh[0],&v_errEffBP_before[0],&v_errEffBP_before[0]);
    gBPEff_before->SetMarkerSize(1.4);
    gBPEff_before->SetMarkerColor(kBlack);

    TGraphAsymmErrors *gBPEff_after = new TGraphAsymmErrors(vRun_average.size(),&vRun_average[0],&v_effBP_after[0],&vRunLow[0],&vRunHigh[0],&v_errEffBP_after[0],&v_errEffBP_after[0]);
    gBPEff_after->SetMarkerSize(1.4);
    gBPEff_after->SetMarkerColor(kRed);

    TGraphAsymmErrors *gNBPEff_before = new TGraphAsymmErrors(vRun_average.size(),&vRun_average[0],&v_effNBP_before[0],&vRunLow[0],&vRunHigh[0],&v_errEffNBP_before[0],&v_errEffNBP_before[0]);
    gNBPEff_before->SetMarkerSize(1.4);
    gNBPEff_before->SetMarkerColor(kBlack);

    TGraphAsymmErrors *gNBPEff_after = new TGraphAsymmErrors(vRun_average.size(),&vRun_average[0],&v_effNBP_after[0],&vRunLow[0],&vRunHigh[0],&v_errEffNBP_after[0],&v_errEffNBP_after[0]);
    gNBPEff_after->SetMarkerSize(1.4);
    gNBPEff_after->SetMarkerColor(kRed);

    TGraphAsymmErrors *gBothEff_before = new TGraphAsymmErrors(vRun_average.size(),&vRun_average[0],&v_effBoth_before[0],&vRunLow[0],&vRunHigh[0],&v_errEffBoth_before[0],&v_errEffBoth_before[0]);
    gBothEff_before->SetMarkerSize(1.4);
    gBothEff_before->SetMarkerColor(kBlack);

    TGraphAsymmErrors *gBothEff_after = new TGraphAsymmErrors(vRun_average.size(),&vRun_average[0],&v_effBoth_after[0],&vRunLow[0],&vRunHigh[0],&v_errEffBoth_after[0],&v_errEffBoth_after[0]);
    gBothEff_after->SetMarkerSize(1.4);
    gBothEff_after->SetMarkerColor(kRed);

    if (merge) {
        gBPEff_before->SetMarkerStyle(8);
        gBPEff_after->SetMarkerStyle(8);
        gNBPEff_before->SetMarkerStyle(8);
        gNBPEff_after->SetMarkerStyle(8);
        gBothEff_before->SetMarkerStyle(8);
        gBothEff_after->SetMarkerStyle(8);
    }

    else {
        gBPEff_before->SetMarkerStyle(22);
        gBPEff_after->SetMarkerStyle(22);
        gNBPEff_before->SetMarkerStyle(22);
        gNBPEff_after->SetMarkerStyle(22);
        gBothEff_before->SetMarkerStyle(22);
        gBothEff_after->SetMarkerStyle(22);
    }

    TMultiGraph *mBP = new TMultiGraph();
    mBP->Add(gBPEff_before);
    mBP->Add(gBPEff_after);
    
    TMultiGraph *mNBP = new TMultiGraph();
    mNBP->Add(gNBPEff_before);
    mNBP->Add(gNBPEff_after);

    TMultiGraph *mBoth = new TMultiGraph();
    mBoth->Add(gBothEff_before);
    mBoth->Add(gBothEff_after);

    hRun.close();
    //hDate.close();

    TCanvas *cEff_BP = new TCanvas();
    cEff_BP->cd();
    cEff_BP->SetGridx();
    cEff_BP->SetGridy();
    mBP->Draw("AP");
    mBP->GetXaxis()->SetNoExponent(1);
    mBP->GetXaxis()->SetTitle("Run #");
    mBP->GetXaxis()->CenterTitle(true);
    mBP->GetYaxis()->SetTitle("Efficiency [%]");
    mBP->GetYaxis()->CenterTitle(true);
    mBP->GetXaxis()->SetLabelFont(62); 
    mBP->GetYaxis()->SetLabelFont(62); 
    mBP->GetYaxis()->SetTitleFont(62); 
    mBP->GetYaxis()->SetTitleFont(62); 

    TCanvas *cEff_NBP = new TCanvas();
    cEff_NBP->cd();
    cEff_NBP->SetGridx();
    cEff_NBP->SetGridy();
    mNBP->Draw("AP");
    mNBP->GetXaxis()->SetNoExponent(1);
    mNBP->GetXaxis()->SetTitle("Run #");
    mNBP->GetXaxis()->CenterTitle(true);
    mNBP->GetYaxis()->SetTitle("Efficiency [%]");
    mNBP->GetYaxis()->CenterTitle(true);
    mNBP->GetXaxis()->SetLabelFont(62); 
    mNBP->GetYaxis()->SetLabelFont(62); 
    mNBP->GetYaxis()->SetTitleFont(62); 
    mNBP->GetYaxis()->SetTitleFont(62);

    TCanvas *cEff_Both = new TCanvas();
    cEff_Both->cd();
    cEff_Both->SetGridx();
    cEff_Both->SetGridy();
    mBoth->Draw("AP");
    mBoth->GetXaxis()->SetNoExponent(1);
    mBoth->GetXaxis()->SetTitle("Run #");
    mBoth->GetXaxis()->CenterTitle(true);
    mBoth->GetYaxis()->SetTitle("Efficiency [%]");
    mBoth->GetYaxis()->CenterTitle(true);
    mBoth->GetXaxis()->SetLabelFont(62); 
    mBoth->GetYaxis()->SetLabelFont(62); 
    mBoth->GetYaxis()->SetTitleFont(62); 
    mBoth->GetYaxis()->SetTitleFont(62);

    TFile *fSummary;

    if (merge) {
        fSummary = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/outFile_" + period + "_merged.root").c_str(),"RECREATE");
    }
    else {
        fSummary = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/outFile_" + period + ".root").c_str(),"RECREATE");
    }

    fSummary->cd();
    cEff_BP->Write("BP");
    cEff_NBP->Write("NBP");
    cEff_Both->Write("Both");

    fSummary->Close();
    //gBPEff_after->Draw("SAME")
}