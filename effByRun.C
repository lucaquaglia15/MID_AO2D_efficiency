#include "TH1F.h"
#include "TH2F.h"
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
#include "TKey.h"

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameters
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

//#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

bool debug = false;

void effByRun() {

    int nBinsPlane = 4; //Number of planes
    int nBinsRPC = 72; //Number of RPCs
    int nBinsBoard = 936; //Number of LBs

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    vector<float> vEffBothLB, vEffBPLB, vEffNBPLB;
    vector<float> vErrEffBothLB, vErrEffBPLB, vErrEffNBPLB;

    vector<vector<float>> vEffBothLB_runs, vEffBPLB_runs, vEffNBPLB_runs;
    vector<vector<float>> vErrEffBothLB_runs, vErrEffBPLB_runs, vErrEffNBPLB_runs;

    /*TH1F *hEffLBplanes1D_both[4]; //4 1D histograms for LB (one per plane eff on both planes)
    TH1F *hEffLBplanes1D_BP[4]; //4 1D histograms for LB (one per plane eff on BP)
    TH1F *hEffLBplanes1D_NBP[4]; //4 1D histograms for LB (one per plane eff on NBP)
    //Initialize them
    for (int i = 0; i < 4; i++) {
        hEffLBplanes1D_both[i] = new TH1F(("LB efficiency both "+planeName[i]).c_str(),("LB efficiency both "+planeName[i]).c_str(),234,0.5,234.5);
        hEffLBplanes1D_BP[i] = new TH1F(("LB efficiency BP "+planeName[i]).c_str(),("LB efficiency BP "+planeName[i]).c_str(),234,0.5,234.5);
        hEffLBplanes1D_NBP[i] = new TH1F(("LB efficiency NBP "+planeName[i]).c_str(),("LB efficiency NBP "+planeName[i]).c_str(),234,0.5,234.5);
    }*/

    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //Path of the merged file, run-by-run
    string runPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/LHC23_pass4_skimmed_QC1/merged_files/"; 

    //Path for the .txt file of the run list of the period
    string runNumbers = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/LHC23_pass4_skimmed_QC1/run_list.txt"; 

    //Open txt file of runs
    ifstream hRun;
    hRun.open(runNumbers.c_str());

    //Push back to a vector of int (no need to care about size)
    float run;
    vector<float> vRun;

    while(hRun >> run) {
        vRun.push_back(run);
    }

    //General string name
    string fileName = "AnalysisResults.root";

    //Loop on all runs
    for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
        //Enter the folder
        string runFolder = runPath+to_string((int)vRun.at(iRun));
        gSystem->cd(runFolder.c_str());

        //run file name = path of the folder + run number (runFolder) + fileName
        string runFileName = runFolder+"/"+fileName;
        TFile *fRun = new TFile(runFileName.c_str(),"READ");

        TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

        TH1F *hFiredBothPlanesLB = (TH1F*)d->Get("nFiredBothperBoard");
        TH1F *hFiredBPLB = (TH1F*)d->Get("nFiredBPperBoard");
        TH1F *hFiredNBPLB = (TH1F*)d->Get("nFiredNBPperBoard");
        TH1F *hTotLB = (TH1F*)d->Get("nTotperBoard");

        for (int i = 1; i <= nBinsBoard; i++) {
            cout << "LB Both planes " <<  i << "\t" << hFiredBothPlanesLB->GetBinContent(i) << "\t" << hTotLB->GetBinContent(i) << endl;

            if (hTotLB->GetBinContent(i) != 0) {

                effBothLB = (hFiredBothPlanesLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effBPLB = (hFiredBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effNBPLB = (hFiredNBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;

                errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLB->GetBinContent(i));
                errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLB->GetBinContent(i));
                errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLB->GetBinContent(i));

                vEffBothLB.push_back(effBothLB);
                vEffBPLB.push_back(effBPLB);
                vEffNBPLB.push_back(effNBPLB);
                
                vErrEffBothLB.push_back(errEffBothLB);
                vErrEffBPLB.push_back(errEffBPLB);
                vErrEffNBPLB.push_back(errEffNBPLB);

                /*hEffLB_both->Fill(i,effBothLB);
                hEffLB_both->SetBinError(i,errEffBothLB);
                hEffLB_BP->Fill(i,effBPLB);
                hEffLB_BP->SetBinError(i,errEffBPLB);
                hEffLB_NBP->Fill(i,effNBPLB);
                hEffLB_NBP->SetBinError(i,errEffNBPLB);

                if (i <= 234) {
                    hEffLBplanes1D_both[0]->SetBinContent(i,effBothLB);
                    hEffLBplanes1D_both[0]->SetBinError(i,errEffBothLB);
                    hEffLBplanes1D_BP[0]->SetBinContent(i,effBPLB);
                    hEffLBplanes1D_BP[0]->SetBinError(i,errEffBPLB);
                    hEffLBplanes1D_NBP[0]->SetBinContent(i,effNBPLB);
                    hEffLBplanes1D_NBP[0]->SetBinError(i,errEffNBPLB);
                }
                else if (i >= 235 && i <= 468) {
                    hEffLBplanes1D_both[1]->SetBinContent(i-234,effBothLB);
                    hEffLBplanes1D_both[1]->SetBinError(i-234,errEffBothLB);
                    hEffLBplanes1D_BP[1]->SetBinContent(i-234,effBPLB);
                    hEffLBplanes1D_BP[1]->SetBinError(i-234,errEffBPLB);
                    hEffLBplanes1D_NBP[1]->SetBinContent(i-234,effNBPLB);
                    hEffLBplanes1D_NBP[1]->SetBinError(i-234,errEffNBPLB);
                }
                else if (i>= 469 && i <= 702) {
                    hEffLBplanes1D_both[2]->SetBinContent(i-468,effBothLB);
                    hEffLBplanes1D_both[2]->SetBinError(i-468,errEffBothLB);
                    hEffLBplanes1D_BP[2]->SetBinContent(i-468,effBPLB);
                    hEffLBplanes1D_BP[2]->SetBinError(i-468,errEffBPLB);
                    hEffLBplanes1D_NBP[2]->SetBinContent(i-468,effNBPLB);
                    hEffLBplanes1D_NBP[2]->SetBinError(i-468,errEffNBPLB);
                }

                else {
                    hEffLBplanes1D_both[3]->SetBinContent(i-702,effBothLB);
                    hEffLBplanes1D_both[3]->SetBinError(i-702,errEffBothLB);
                    hEffLBplanes1D_BP[3]->SetBinContent(i-702,effBPLB);
                    hEffLBplanes1D_BP[3]->SetBinError(i-702,errEffBPLB);
                    hEffLBplanes1D_NBP[3]->SetBinContent(i-702,effNBPLB);
                    hEffLBplanes1D_NBP[3]->SetBinError(i-702,errEffNBPLB);
                }*/
            }
        }

        vEffBothLB_runs.push_back(vEffBothLB);
        vEffBPLB_runs.push_back(vEffBPLB);
        vEffNBPLB_runs.push_back(vEffNBPLB);
         
        vErrEffBothLB_runs.push_back(vErrEffBothLB);
        vErrEffBPLB_runs.push_back(vErrEffBPLB); 
        vErrEffNBPLB_runs.push_back(vErrEffNBPLB);

        vEffBothLB.clear();
        vEffBPLB.clear();
        vEffNBPLB.clear();
            
        vErrEffBothLB.clear();
        vErrEffBPLB.clear();
        vErrEffNBPLB.clear();

        //cout << runFileName << endl;

    } //End of loop on all runs

    cout << vEffBothLB_runs[0].size() << "\t" << vEffBPLB_runs[0].size() << "\t" << vEffNBPLB_runs[0].size() << endl;

    //First [xx] is the run number in the period and the second [xx] is the LB number
    TGraphErrors *gExample = new TGraphErrors(vRun.size(),&vRun[0],&vEffBPLB_runs[0][30],NULL,&vErrEffBPLB_runs[0][30]);
    TGraphErrors *gExample2 = new TGraphErrors(vRun.size(),&vRun[0],&vEffBPLB_runs[22][30],NULL,&vErrEffBPLB_runs[22][30]);
    gExample->SetMarkerStyle(8);
    gExample2->SetMarkerColor(kGreen);
    gExample2->SetMarkerStyle(8);
    gExample2->SetMarkerColor(kRed);

    TMultiGraph *m = new TMultiGraph();
    m->Add(gExample);
    m->Add(gExample2);

    TCanvas *cExample = new TCanvas();
    cExample->cd();
    m->Draw("AP");
}