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
#include "THnSparse.h"
#include "TKey.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include <bits/stdc++.h>

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameter
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

bool debug = false;
int nBinsPlane = 4; //Number of planes
int nBinsRPC = 72; //Number of RPCs
int nBinsBoard = 936; //Number of LBs

//Max and min pt values in GeV/c
float minPt = 0.;
float maxPt = 20.;
int binsPt = 150;
float ptStep = (maxPt-minPt)/binsPt;
int ptMinCut = 1.8/ptStep;

o2::ccdb::CcdbApi api; //CCDB API as global object
o2::mid::Mapping mapping; //MID mapping object to construct ccdb object

void plotEfficiency_MT22IN3() { //Main function

    //Per plane
    double effBothPlane = 0, effBPPlane = 0, effNBPPlane =0;
    double errEffBothPlane = 0, errEffBPPlane = 0, errEffNBPPlane = 0;

    //Per RPC
    double effBothRPC = 0, effBPRPC = 0, effNBPPRPC =0;
    double errEffBothRPC = 0, errEffBPRPC = 0, errEffNBPRPC = 0;

    //Run by run - Plane
    vector<double> vEffBoth_Planes, vEffBP_Planes, vEffNBP_Planes;
    vector<double> vErrEffBoth_Planes, vErrEffBP_Planes, vErrEffNBP_Planes;
    
    vector<vector<double>> vEffBoth_Planes_runs, vEffBP_Planes_runs, vEffNBP_Planes_runs;
    vector<vector<double>> vErrEffBoth_Planes_runs, vErrEffBP_Planes_runs, vErrEffNBP_Planes_runs;
    
    //Run by run - RPC
    vector<double> vEffBoth_RPC, vEffBP_RPC, vEffNBP_RPC;
    vector<double> vErrEffBoth_RPC, vErrEffBP_RPC, vErrEffNBP_RPC;

    vector<double> vEffBoth_RPC_tot, vEffBP_RPC_tot, vEffNBP_RPC_tot;
    vector<double> vErrEffBoth_RPC_tot, vErrEffBP_RPC_tot, vErrEffNBP_RPC_tot;

    vector<vector<double>> vEffBoth_RPC_runs, vEffBP_RPC_runs, vEffNBP_RPC_runs;
    vector<vector<double>> vErrEffBoth_RPC_runs, vErrEffBP_RPC_runs, vErrEffNBP_RPC_runs;

    //start/end/IR/B-field
    double start, end, run, IR;
    vector<double> vStart, vEnd, vRun, vIR;  
    vector<double> vRunTot; //for all runs, not reset after one period, to keep track of all runs and change marker color when B-field polarity changes 
    
    string bField;
    vector<string> vBField;
    vector<string> vBFieldTot; //for all B-field values, not reset after one period, to keep track of all B-fields and change marker color when B-field polarity changes

    bool isIn;
    bool isTime = false;

    //Legend for all periods
    TLegend *lPeriods = new TLegend(0.6,0.6,0.9,0.9,"","rNDC"); 

    //Multigraph for different periods with different colors
    TMultiGraph *mEffBothPeriod = new TMultiGraph();
    TMultiGraph *mEffBPPeriod = new TMultiGraph();
    TMultiGraph *mEffNBPPeriod = new TMultiGraph();
    
    //General output file name
    string fileName = "AnalysisResults.root";

    bool ptCut = true; //If true, a pt cut is also applied in the computation of the efficiency vs eta

    vector<string> vPeriods = {
        "LHC22_pass7_skimmed",
        "LHC23_pass4_skimmed",
        "LHC23_PbPb_pass4",
        "LHC24_pass1_skimmed",
        "LHC24_PbPb_pass1",
        "LHC25ac_pass1_skimmed",
        "LHC25ad_pass2",
        "LHC25ae_pass2",
        "LHC25af_pass2",
        "LHC25ah_pass1_skimmed_small",
        "LHC25ai_pass1_skimmed",
        "LHC25an_cpass0_QC1_sampling"
    };

    //For now we are interested in RPC MT 22 IN 3
    //-> RPC numbering scheme (0-71) -> 29 (https://github.com/ariffero/aQC-studies/blob/master/MID%20geometry.pdf)
    int binNumber = 30; //bin 1 is RPC 0

    vector<bool> isAnalyzed = {1,1,1,1,1,1,1,1,1,1,1,1};
    //vector<bool> isAnalyzed = {0,0,0,0,1,0,0,0,0,0,0,0};
    vector<bool> isPbPb = {false,false,true,false,true,false,true,true,true,false,false,true};
    vector<bool> isAnalyzedPbPb;
    vector<int> color = {1,2,3,4,5,6,7,8,9,12,16,20};

    //vector<int> markerStyle = {89,90,91,92,93,94,95,96};
    //vector<int> markerStyle = {71,72,73,74,75,76,77,78};
    //vector<int> markerStyle = {53,54,55,56,57,58,59,60};

    //Check that vPeriods and isAnalyzed have the same number of elements
    if (vPeriods.size() != isAnalyzed.size()) {
        cout << "Careful,  vPeriods and isAnalyzed don't have the same number of elements!" << endl;
        //break;
    }

    vector<string> analyzedPeriod;

    for (unsigned int i = 0; i < vPeriods.size(); i++) {
        if (isAnalyzed.at(i)) {
            analyzedPeriod.push_back(vPeriods.at(i));
            isAnalyzedPbPb.push_back(isPbPb.at(i));
        }
    }

    //Vector of TGraphErrors for all 4 planes on both, BP and NBP - plane eff vs run #
    vector<TGraphErrors*> gEffPerRunPlaneBoth, gEffPerRunPlaneBP, gEffPerRunPlaneNBP;

    for (unsigned int per = 0; per < analyzedPeriod.size(); per++) {
        //base path for all periods
        string globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+analyzedPeriod.at(per)+"/";

        //Path of the merged file, run-by-run
        string runPath = globalPath+"runs/";
        
        //runs + info on magnetic field + IR (for PbPb runs)
        string fIR = globalPath+"run_IR_Bfield.txt"; 
        ifstream hIr;
        hIr.open(fIR.c_str()); 

        while (hIr >> isIn >> run >> IR >> bField >> start >> end) {
            if (isIn) {
                vRun.push_back(run);
                vRunTot.push_back(run);
                vIR.push_back(IR/1000);
                vBField.push_back(bField);
                vBFieldTot.push_back(bField);
                vStart.push_back(start);
                vEnd.push_back(end);
            }
        }

        cout << "Analyzing period: " << analyzedPeriod.at(per) << endl;

        /*for (unsigned int zz = 0; zz < vRun.size(); zz++) {
            if (zz < vRun.size() - 1) {
                cout << vRun.at(zz) << ",";
            }
            else {
                cout << vRun.at(zz) << endl;
            }
        }*/

        //Loop on all runs
        for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
            cout << "Processing run # " << vRun.at(iRun) << endl;
            
            //Enter the folder
            string runFolder = runPath+to_string((int)vRun.at(iRun));

            //run file name = path of the folder + run number (runFolder) + fileName
            string runFileName = runFolder+"/"+fileName;
            TFile *fRun = new TFile(runFileName.c_str(),"READ");

            TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

            //Histo of counts per RPC
            THnSparse *hSparseCentFiredTotPerRPC = (THnSparse*)d->Get("hSparseCentFiredTotperRPC");
            THnSparse *hSparseCentFiredBothPerRPC = (THnSparse*)d->Get("hSparseCentFiredBothperRPC");
            THnSparse *hSparseCentFiredBPPerRPC = (THnSparse*)d->Get("hSparseCentFiredBPperRPC");
            THnSparse *hSparseCentFiredNBPPerRPC = (THnSparse*)d->Get("hSparseCentFiredNBPperRPC");

            //Add a cut in pt in the efficiency computation (pt > 2 GeV/c)
            if (ptCut) {

                int ptAxis;

                if (isAnalyzedPbPb.at(per) == true) {
                    ptAxis = 2;
                }
                else {
                    ptAxis = 1;
                }
                
                //cout << "Cutting on pt axis (2 for PbPb, 1 for pp): " << ptAxis << endl;
                hSparseCentFiredTotPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
                hSparseCentFiredBothPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
                hSparseCentFiredBPPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
                hSparseCentFiredNBPPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
            }

            //Project THNsparse            
            TH1D* totRPCCountsProj = hSparseCentFiredTotPerRPC->Projection(0); 
            TH1D* BothRPCCountsProj = hSparseCentFiredBothPerRPC->Projection(0);
            TH1D* BPRPCCountsProj = hSparseCentFiredBPPerRPC->Projection(0);
            TH1D* NBPRPCCountsProj = hSparseCentFiredNBPPerRPC->Projection(0);

            effBothRPC = (BothRPCCountsProj->GetBinContent(binNumber)/totRPCCountsProj->GetBinContent(binNumber))*100;
            effBPRPC = (BPRPCCountsProj->GetBinContent(binNumber)/totRPCCountsProj->GetBinContent(binNumber))*100;
            effNBPPRPC = (NBPRPCCountsProj->GetBinContent(binNumber)/totRPCCountsProj->GetBinContent(binNumber))*100;

            errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/totRPCCountsProj->GetBinContent(binNumber));
            errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/totRPCCountsProj->GetBinContent(binNumber));
            errEffNBPRPC = TMath::Sqrt(effNBPPRPC*(100-effNBPPRPC)/totRPCCountsProj->GetBinContent(binNumber));

            //Fill vector for efficiency per LB in the run
            vEffBoth_RPC.push_back(effBothRPC);
            vEffBP_RPC.push_back(effBPRPC);
            vEffNBP_RPC.push_back(effNBPPRPC);
            
            //Fill vector for error on efficiency per LB in the run
            vErrEffBoth_RPC.push_back(errEffBothRPC);
            vErrEffBP_RPC.push_back(errEffBPRPC);
            vErrEffNBP_RPC.push_back(errEffNBPRPC);
            
            //-----//
            delete totRPCCountsProj;
            delete BothRPCCountsProj;
            delete BPRPCCountsProj;
            delete NBPRPCCountsProj; 

            //-----//
            delete hSparseCentFiredTotPerRPC; 
            delete hSparseCentFiredBothPerRPC;
            delete hSparseCentFiredBPPerRPC; 
            delete hSparseCentFiredNBPPerRPC;

            fRun->Close();
            delete fRun;

        } //End of loop on all runs in a given period

        cout << "Size of vStart in period: " << analyzedPeriod.at(per) << ": " << vStart.size() << endl;
        cout << "Size of vRun " << analyzedPeriod.at(per) << ": " << vRun.size() << endl;
        cout << "Size of vEffBoth_Planes " << analyzedPeriod.at(per) << ": " << vEffBoth_Planes_runs.size() << endl;
        cout << "size of B field vector " << analyzedPeriod.at(per) << ": " << vBField.size() << endl;

        TGraphErrors *gBoth = new TGraphErrors(vRun.size(),&vRun[0],&vEffBoth_RPC[0],NULL,&vErrEffBoth_RPC[0]);
        TGraphErrors *gBP = new TGraphErrors(vRun.size(),&vRun[0],&vEffBP_RPC[0],NULL,&vErrEffBP_RPC[0]);
        TGraphErrors *gNBP = new TGraphErrors(vRun.size(),&vRun[0],&vEffNBP_RPC[0],NULL,&vErrEffNBP_RPC[0]);

        gBoth->SetMarkerStyle(8);
        gBoth->SetMarkerColor(color.at(per));

        gBP->SetMarkerStyle(8);
        gBP->SetMarkerColor(color.at(per));

        gNBP->SetMarkerStyle(8);
        gNBP->SetMarkerColor(color.at(per));
        
        mEffBothPeriod->Add(gBoth);
        mEffBPPeriod->Add(gBP);
        mEffNBPPeriod->Add(gNBP);

        TMarker *m = new TMarker();
        m->SetMarkerStyle(8);
        m->SetMarkerSize(1.5);
        m->SetMarkerColor(color.at(per));
        lPeriods->AddEntry(m,analyzedPeriod.at(per).c_str(),"p");

        vRun.clear();
        vEffBoth_RPC.clear();
        vEffBP_RPC.clear();
        vEffNBP_RPC.clear();
        vErrEffBoth_RPC.clear();
        vErrEffBP_RPC.clear();
        vErrEffNBP_RPC.clear();
        //vRunTot.push_back(run);
        vIR.clear();
        vBField.clear();
        //vBFieldTot.push_back(bField);
        vStart.clear();
        vEnd.clear();
    } //End of loop on all periods
    
    //cout << vEffPerPlanePerRunBoth_PbPb[0].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[1].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[2].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[3].size() << endl;
    //cout << vEffPerPlanePerRunBP_PbPb[0].size() << "\t" << vEffPerPlanePerRunBP_PbPb[1].size() << "\t" << vEffPerPlanePerRunBP_PbPb[2].size() << "\t" << vEffPerPlanePerRunBP_PbPb[3].size() << endl;
    //cout << vEffPerPlanePerRunNBP_PbPb[0].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[1].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[2].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[3].size() << endl;

    TGraphErrors *gEffRPCrunBoth = new TGraphErrors(vRunTot.size(),&vRunTot[0],&vEffBoth_RPC_tot[0],NULL,&vErrEffBoth_RPC_tot[0]);
    TGraphErrors *gEffRPCrunBP = new TGraphErrors(vRunTot.size(),&vRunTot[0],&vEffBP_RPC_tot[0],NULL,&vErrEffBP_RPC_tot[0]);
    TGraphErrors *gEffRPCrunNBP = new TGraphErrors(vRunTot.size(),&vRunTot[0],&vEffNBP_RPC_tot[0],NULL,&vErrEffNBP_RPC_tot[0]);

    //Both
    //TCanvas *cEffPerPlaneBoth = new TCanvas();
    //cEffPerPlaneBoth->Divide(1,4);
    TCanvas *cEffPerPlaneBoth = new TCanvas();
    TCanvas *cEffPerPlaneBP = new TCanvas();
    TCanvas *cEffPerPlaneNBP = new TCanvas();
    
        
    cEffPerPlaneBoth->cd();
    cEffPerPlaneBoth->SetGridx();
    cEffPerPlaneBoth->SetGridy();
    gEffRPCrunBoth->SetTitle("Both planes efficiency");
    gEffRPCrunBoth->GetXaxis()->SetTitle("Run #");
    gEffRPCrunBoth->GetXaxis()->SetTitleOffset(0.9);
    gEffRPCrunBoth->GetXaxis()->SetTitleSize(0.04);
    gEffRPCrunBoth->GetXaxis()->SetTitleFont(62);
    gEffRPCrunBoth->GetXaxis()->SetLabelSize(0.04);
    gEffRPCrunBoth->GetXaxis()->SetLabelFont(62);
    gEffRPCrunBoth->GetXaxis()->SetNoExponent(1);
    gEffRPCrunBoth->GetYaxis()->SetTitle("Efficiency [%]");
    gEffRPCrunBoth->GetYaxis()->SetTitleOffset(0.85);
    gEffRPCrunBoth->GetYaxis()->SetTitleSize(0.04);
    gEffRPCrunBoth->GetYaxis()->SetTitleFont(62);
    gEffRPCrunBoth->GetYaxis()->SetLabelSize(0.04);
    gEffRPCrunBoth->GetYaxis()->SetLabelFont(62); 
    gEffRPCrunBoth->GetYaxis()->SetRangeUser(80,100); 
    gEffRPCrunBoth->Draw("AP");
        
        
    cEffPerPlaneBP->cd();
    cEffPerPlaneBP->SetGridx();
    cEffPerPlaneBP->SetGridy();
    gEffRPCrunBP->SetTitle("BP efficiency");
    gEffRPCrunBP->GetXaxis()->SetTitle("Run #");
    gEffRPCrunBP->GetXaxis()->SetTitleOffset(0.9);
    gEffRPCrunBP->GetXaxis()->SetTitleSize(0.04);
    gEffRPCrunBP->GetXaxis()->SetTitleFont(62);
    gEffRPCrunBP->GetXaxis()->SetLabelSize(0.04);
    gEffRPCrunBP->GetXaxis()->SetLabelFont(62);
    gEffRPCrunBP->GetXaxis()->SetNoExponent(1);
    gEffRPCrunBP->GetYaxis()->SetTitle("Efficiency [%]");
    gEffRPCrunBP->GetYaxis()->SetTitleOffset(0.85);
    gEffRPCrunBP->GetYaxis()->SetTitleSize(0.04);
    gEffRPCrunBP->GetYaxis()->SetTitleFont(62);
    gEffRPCrunBP->GetYaxis()->SetLabelSize(0.04);
    gEffRPCrunBP->GetYaxis()->SetLabelFont(62); 
    gEffRPCrunBP->GetYaxis()->SetRangeUser(80,100);
    gEffRPCrunBP->Draw("AP");
        
    
    cEffPerPlaneNBP->cd();
    cEffPerPlaneNBP->SetGridx();
    cEffPerPlaneNBP->SetGridy();
    gEffRPCrunNBP->SetTitle("NBP efficiency");
    gEffRPCrunNBP->GetXaxis()->SetTitle("Run #");
    gEffRPCrunNBP->GetXaxis()->SetTitleOffset(0.9);
    gEffRPCrunNBP->GetXaxis()->SetTitleSize(0.04);
    gEffRPCrunNBP->GetXaxis()->SetTitleFont(62);
    gEffRPCrunNBP->GetXaxis()->SetLabelSize(0.04);
    gEffRPCrunNBP->GetXaxis()->SetLabelFont(62);
    gEffRPCrunNBP->GetXaxis()->SetNoExponent(1);
    gEffRPCrunNBP->GetYaxis()->SetTitle("Efficiency [%]");
    gEffRPCrunNBP->GetYaxis()->SetTitleOffset(0.85);
    gEffRPCrunNBP->GetYaxis()->SetTitleSize(0.04);
    gEffRPCrunNBP->GetYaxis()->SetTitleFont(62);
    gEffRPCrunNBP->GetYaxis()->SetLabelSize(0.04);
    gEffRPCrunNBP->GetYaxis()->SetLabelFont(62);
    gEffRPCrunNBP->GetYaxis()->SetRangeUser(80,100); 
    gEffRPCrunNBP->Draw("AP"); 

    TCanvas *cEffBoth_period = new TCanvas();
    cEffBoth_period->SetGridx();
    cEffBoth_period->SetGridy();
    mEffBothPeriod->Draw("AP");
    mEffBothPeriod->SetTitle("Both planes");
    mEffBothPeriod->GetXaxis()->SetTitle("Run #");
    mEffBothPeriod->GetXaxis()->SetTitleOffset(0.9);
    mEffBothPeriod->GetXaxis()->SetTitleSize(0.04);
    mEffBothPeriod->GetXaxis()->SetTitleFont(62);
    mEffBothPeriod->GetXaxis()->SetLabelSize(0.04);
    mEffBothPeriod->GetXaxis()->SetLabelFont(62);
    mEffBothPeriod->GetXaxis()->SetNoExponent(1);
    mEffBothPeriod->GetYaxis()->SetTitle("Efficiency [%]");
    mEffBothPeriod->GetYaxis()->SetTitleOffset(0.85);
    mEffBothPeriod->GetYaxis()->SetTitleSize(0.04);
    mEffBothPeriod->GetYaxis()->SetTitleFont(62);
    mEffBothPeriod->GetYaxis()->SetLabelSize(0.04);
    mEffBothPeriod->GetYaxis()->SetLabelFont(62);
    mEffBothPeriod->GetYaxis()->SetRangeUser(80,100); 
    lPeriods->Draw("SAME");

    TCanvas *cEffBP_period = new TCanvas();
    cEffBP_period->SetGridx();
    cEffBP_period->SetGridy();
    mEffBPPeriod->Draw("AP");
    mEffBPPeriod->SetTitle("BP");
    mEffBPPeriod->GetXaxis()->SetTitle("Run #");
    mEffBPPeriod->GetXaxis()->SetTitleOffset(0.9);
    mEffBPPeriod->GetXaxis()->SetTitleSize(0.04);
    mEffBPPeriod->GetXaxis()->SetTitleFont(62);
    mEffBPPeriod->GetXaxis()->SetLabelSize(0.04);
    mEffBPPeriod->GetXaxis()->SetLabelFont(62);
    mEffBPPeriod->GetXaxis()->SetNoExponent(1);
    mEffBPPeriod->GetYaxis()->SetTitle("Efficiency [%]");
    mEffBPPeriod->GetYaxis()->SetTitleOffset(0.85);
    mEffBPPeriod->GetYaxis()->SetTitleSize(0.04);
    mEffBPPeriod->GetYaxis()->SetTitleFont(62);
    mEffBPPeriod->GetYaxis()->SetLabelSize(0.04);
    mEffBPPeriod->GetYaxis()->SetLabelFont(62);
    mEffBPPeriod->GetYaxis()->SetRangeUser(80,100); 
    lPeriods->Draw("SAME");

    TCanvas *cEffNBP_period = new TCanvas();
    cEffNBP_period->SetGridx();
    cEffNBP_period->SetGridy();
    mEffNBPPeriod->Draw("AP");
    mEffNBPPeriod->SetTitle("NBP");
    mEffNBPPeriod->GetXaxis()->SetTitle("Run #");
    mEffNBPPeriod->GetXaxis()->SetTitleOffset(0.9);
    mEffNBPPeriod->GetXaxis()->SetTitleSize(0.04);
    mEffNBPPeriod->GetXaxis()->SetTitleFont(62);
    mEffNBPPeriod->GetXaxis()->SetLabelSize(0.04);
    mEffNBPPeriod->GetXaxis()->SetLabelFont(62);
    mEffNBPPeriod->GetXaxis()->SetNoExponent(1);
    mEffNBPPeriod->GetYaxis()->SetTitle("Efficiency [%]");
    mEffNBPPeriod->GetYaxis()->SetTitleOffset(0.85);
    mEffNBPPeriod->GetYaxis()->SetTitleSize(0.04);
    mEffNBPPeriod->GetYaxis()->SetTitleFont(62);
    mEffNBPPeriod->GetYaxis()->SetLabelSize(0.04);
    mEffNBPPeriod->GetYaxis()->SetLabelFont(62);
    mEffNBPPeriod->GetYaxis()->SetRangeUser(80,100);
    lPeriods->Draw("SAME"); 
    

    TFile *fOutTot = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/MT22IN3/outFile_binNumber_" + to_string(binNumber) + ".root").c_str(),"RECREATE");
    fOutTot->cd();

    cEffPerPlaneBoth->Write(("efficiency_binNumber_" + to_string(binNumber) + "_bothPlanes").c_str());
    cEffPerPlaneBP->Write(("efficiency_binNumber_" + to_string(binNumber) + "_BP").c_str());
    cEffPerPlaneNBP->Write(("efficiency_binNumber_" + to_string(binNumber) + "_NBP").c_str());
    
    //Different color
    cEffBoth_period->Write(("efficiency_binNumber_" + to_string(binNumber) + "_bothPlanes_color").c_str());
    cEffBP_period->Write(("efficiency_binNumber_" + to_string(binNumber) + "BP_colors").c_str());
    cEffNBP_period->Write(("efficiency_binNumber_" + to_string(binNumber) + "_NBP_colors").c_str());

    fOutTot->Close();
    
}