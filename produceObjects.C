#include "TH1F.h"
#include "TH2F.h"
//#include "TROOT.h"
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

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameters
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

int nBinsPlane = 4; //Number of planes
int nBinsRPC = 72; //Number of RPCs
int nBinsBoard = 936; //Number of LBs

o2::ccdb::CcdbApi api; //CCDB API as global object
o2::mid::Mapping mapping; //MID mapping object to construct ccdb object

//Function to upload to ccdb
void createCCDB(TH1F *hFiredBPLB, TH1F *hFiredNBPLB, TH1F *hFiredBothPlanesLB, TH1F *hTotLB, long int startValidity, long int endValidity, int runNumber) {

    bool debug = false;

    //Variables used in the loop
    int LB936 = 0; //LB 1-> 936
    int plane = 0; //plane

    std::vector<o2::mid::ChEffCounter> counterVector; //Vector of efficiency counter structs
    struct o2::mid::ChEffCounter entry; //Struct for efficiency counters
    o2::mid::ChamberEfficiency effMap; //Chamber efficiency map

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

                entry.counts[0] = hFiredBPLB->GetBinContent(LB936); //BP
                entry.counts[1] = hFiredNBPLB->GetBinContent(LB936); //NBP
                entry.counts[2] = hFiredBothPlanesLB->GetBinContent(LB936); //Both
                entry.counts[3] = hTotLB->GetBinContent(LB936); //Total

                counterVector.push_back(entry);
                
                //Debug printout
                if (debug) {
                    cout << LB936 << "\t both: " << hFiredBothPlanesLB->GetBinContent(LB936) << "\t tot: " << hTotLB->GetBinContent(LB936) << endl;
                } 
            }
        }
    }

    //TFile *fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/LHC22_pass7_skimmed/ccdb/o2-mid-ChEffCounter_"+to_string(runNumber)+".root").c_str(),"RECREATE");
    TFile *fOut = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/LHC23_pass4_skimmed/ccdb/o2-mid-ChEffCounter_"+to_string(runNumber)+".root").c_str(),"RECREATE");
    fOut->WriteObjectAny(&counterVector, "std::vector<o2::mid::ChEffCounter>","ccdb-object");
    fOut->Close();
    
} //end of createCCDB


void produceObjects() { //Main function

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    //string period = "LHC22_pass7_skimmed";
    string period = "LHC23_pass4_skimmed";
    string globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/";

    string fileName = "AnalysisResults.root";
    
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

    cout << vRun.size() << vStart.size() << vEnd.size() << endl;

    //Loop on all runs
    for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
        //Enter the folder
        string runFolder = globalPath + "runs/" +to_string(vRun.at(iRun));
        cout << runFolder << endl;

        //run file name = path of the folder + run number (runFolder) + fileName
        string runFileName = runFolder+"/"+fileName;
        TFile *fRun = new TFile(runFileName.c_str(),"READ");

        TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

        TH1F *hFiredBothPlanesLB = (TH1F*)d->Get("nFiredBothperBoard");
        TH1F *hFiredBPLB = (TH1F*)d->Get("nFiredBPperBoard");
        TH1F *hFiredNBPLB = (TH1F*)d->Get("nFiredNBPperBoard");
        TH1F *hTotLB = (TH1F*)d->Get("nTotperBoard");

        //Create .root objects to upload them to CCDB
        bool createByRun = false;
        if (createByRun) {
            createCCDB(hFiredBPLB,hFiredNBPLB,hFiredBothPlanesLB,hTotLB,vStart.at(iRun),vEnd.at(iRun),vRun.at(iRun));
        }

        //Calculate efficiency per run -> not needed here but just to print it out for debug purposes
        for (int i = 1; i <= nBinsBoard; i++) {

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
        }
    } //End of loop on all runs
    
    hRun.close();
    //hDate.close();
}