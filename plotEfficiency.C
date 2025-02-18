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
#include "TMarker.h"
#include "TLegend.h"
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

o2::ccdb::CcdbApi api; //CCDB API as global object
o2::mid::Mapping mapping; //MID mapping object to construct ccdb object

void plotEfficiency() { //Main function

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
    
    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //General output file name
    string fileName = "AnalysisResults.root";

    //Legends for eff vs run #
    TLegend *lAllAnalyzedPeriodsBoth = new TLegend(0.13,0.38,0.2,0.46,"","rNDC");
    lAllAnalyzedPeriodsBoth->SetBorderSize(0);	//No borders in legend
    lAllAnalyzedPeriodsBoth->SetFillStyle(0);		//Transparent background
    lAllAnalyzedPeriodsBoth->SetTextFont(62);		//Bold legend
    lAllAnalyzedPeriodsBoth->SetTextSize(0.025);
    TLegend *lAllAnalyzedPeriodsBP = new TLegend(0.13,0.38,0.2,0.46,"","rNDC");
    lAllAnalyzedPeriodsBP->SetBorderSize(0);	//No borders in legend
    lAllAnalyzedPeriodsBP->SetFillStyle(0);		//Transparent background
    lAllAnalyzedPeriodsBP->SetTextFont(62);		//Bold legend
    lAllAnalyzedPeriodsBP->SetTextSize(0.025);
    TLegend *lAllAnalyzedPeriodsNBP = new TLegend(0.13,0.38,0.2,0.46,"","rNDC");
    lAllAnalyzedPeriodsNBP->SetBorderSize(0);	//No borders in legend
    lAllAnalyzedPeriodsNBP->SetFillStyle(0);		//Transparent background
    lAllAnalyzedPeriodsNBP->SetTextFont(62);		//Bold legend
    lAllAnalyzedPeriodsNBP->SetTextSize(0.025);
    //General path to add flexibility to the code + period name
    string period_DQ_LHC23_PbPb_pass4 = "DQ_LHC23_PbPb_pass4";
    string period_LHC22o_pass7_minBias = "LHC22o_pass7_minBias";
    string period_LHC22_pass7_skimmed = "LHC22_pass7_skimmed";
    string period_LHC23_pass4_skimmed = "LHC23_pass4_skimmed";
    string period_LHC23_PbPb_pass3_I_A11 = "LHC23_PbPb_pass3_I-A11";
    string period_LHC23_PbPb_pass4 = "LHC23_PbPb_pass4";
    string period_LHC24_pass1_skimmed = "LHC24_pass1_skimmed";
    string period_LHC24_PbPb_pass1 = "LHC24_PbPb_pass1";

    vector<string> vPeriods = {"DQ_LHC23_PbPb_pass4_muon","LHC22o_pass7_minBias","LHC22_pass7_skimmed","LHC23_pass4_skimmed","LHC23_PbPb_pass3_I-A11","LHC23_PbPb_pass4","LHC24_pass1_skimmed","LHC24_PbPb_pass1"};

    vector<bool> isAnalyzed = {0,0,1,1,0,1,1,0};

    //vector<int> markerStyle = {89,90,91,92,93,94,95,96};
    //vector<int> markerStyle = {71,72,73,74,75,76,77,78};
    vector<int> markerStyle = {53,54,55,56,57,58,59,60};
    
    //Each element of this vector is a plane
    vector<TMultiGraph*> mEffPlaneBoth, mEffPlaneBP, mEffPlaneNBP;

    for (unsigned int i = 0; i < nBinsPlane; i++) {
        TMultiGraph *m1 = new TMultiGraph();
        TMultiGraph *m2 = new TMultiGraph();
        TMultiGraph *m3 = new TMultiGraph();  

        mEffPlaneBoth.push_back(m1);
        mEffPlaneBP.push_back(m2);
        mEffPlaneNBP.push_back(m3);
    }

    //Check that vPeriods and isAnalyzed have the same number of elements
    if (vPeriods.size() != isAnalyzed.size()) {
        cout << "Careful,  vPeriods and isAnalyzed don't have the same number of elements!" << endl;
        //break;
    }

    vector<string> analyzedPeriod;

    for (unsigned int i = 0; i < vPeriods.size(); i++) {
        if (isAnalyzed.at(i)) {
            analyzedPeriod.push_back(vPeriods.at(i));
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

        for (unsigned int zz = 0; zz < vRun.size(); zz++) {
            cout << vRun.at(zz) << endl;
        }

        //Loop on all runs
        for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
            //Enter the folder
            string runFolder = runPath+to_string((int)vRun.at(iRun));

            //run file name = path of the folder + run number (runFolder) + fileName
            string runFileName = runFolder+"/"+fileName;
            TFile *fRun = new TFile(runFileName.c_str(),"READ");

            TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

            //Histo of counts per plane
            TH1F *hFiredBoth_Planes = (TH1F*)d->Get("nFiredBothperPlane");
            TH1F *hFiredBP_Planes = (TH1F*)d->Get("nFiredBPperPlane");
            TH1F *hFiredNBP_Planes = (TH1F*)d->Get("nFiredNBPperPlane");
            TH1F *hTotPlanes = (TH1F*)d->Get("nTotperPlane");    
            //Histo of counts per RPC
            TH1F *hFiredBoth_RPC = (TH1F*)d->Get("nFiredBothperRPC");
            TH1F *hFiredBP_RPC = (TH1F*)d->Get("nFiredBPperRPC");
            TH1F *hFiredNBP_RPC = (TH1F*)d->Get("nFiredNBPperRPC");
            TH1F *hTotRPC = (TH1F*)d->Get("nTotperRPC"); 

            //Calculate eff per run per plane
            for (int i = 1; i <= nBinsPlane; i++) {

                effBothPlane = (hFiredBoth_Planes->GetBinContent(i)/hTotPlanes->GetBinContent(i))*100;
                effBPPlane = (hFiredBP_Planes->GetBinContent(i)/hTotPlanes->GetBinContent(i))*100;
                effNBPPlane = (hFiredNBP_Planes->GetBinContent(i)/hTotPlanes->GetBinContent(i))*100;

                errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/hTotPlanes->GetBinContent(i));
                errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/hTotPlanes->GetBinContent(i));
                errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/hTotPlanes->GetBinContent(i));

                //Fill vector for efficiency per LB in the run
                vEffBoth_Planes.push_back(effBothPlane);
                vEffBP_Planes.push_back(effBPPlane);
                vEffNBP_Planes.push_back(effNBPPlane);
                
                //Fill vector for error on efficiency per LB in the run
                vErrEffBoth_Planes.push_back(errEffBothPlane);
                vErrEffBP_Planes.push_back(errEffBPPlane);
                vErrEffNBP_Planes.push_back(errEffNBPPlane);
            
            }

            //Calculate eff per run per RPC
            for (int i = 1; i <= nBinsRPC; i++) {

                effBothRPC = (hFiredBoth_RPC->GetBinContent(i)/hTotRPC->GetBinContent(i))*100;
                effBPRPC = (hFiredBP_RPC->GetBinContent(i)/hTotRPC->GetBinContent(i))*100;
                effNBPPRPC = (hFiredNBP_RPC->GetBinContent(i)/hTotRPC->GetBinContent(i))*100;

                errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/hTotRPC->GetBinContent(i));
                errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/hTotRPC->GetBinContent(i));
                errEffNBPRPC = TMath::Sqrt(effNBPPRPC*(100-effNBPPRPC)/hTotRPC->GetBinContent(i));

                //Fill vector for efficiency per LB in the run
                vEffBoth_RPC.push_back(effBothRPC);
                vEffBP_RPC.push_back(effBPRPC);
                vEffNBP_RPC.push_back(effNBPPRPC);
                
                //Fill vector for error on efficiency per LB in the run
                vErrEffBoth_RPC.push_back(errEffBothRPC);
                vErrEffBP_RPC.push_back(errEffBPRPC);
                vErrEffNBP_RPC.push_back(errEffNBPRPC);
            
            }

            //Push back the vector with the eff of LB to a larger vector of vectors (one element of this = one run) - per plane
            vEffBoth_Planes_runs.push_back(vEffBoth_Planes);
            vEffBP_Planes_runs.push_back(vEffBP_Planes);
            vEffNBP_Planes_runs.push_back(vEffNBP_Planes);
            
            //Push back the vector with the error on eff of LB to a larger vector of vectors (one element of this = one run) - per plane
            vErrEffBoth_Planes_runs.push_back(vErrEffBoth_Planes);
            vErrEffBP_Planes_runs.push_back(vErrEffBP_Planes); 
            vErrEffNBP_Planes_runs.push_back(vErrEffNBP_Planes);

            //Push back the vector with the eff of LB to a larger vector of vectors (one element of this = one run) - per RPC
            vEffBoth_RPC_runs.push_back(vEffBoth_RPC);
            vEffBP_RPC_runs.push_back(vEffBP_RPC);
            vEffNBP_RPC_runs.push_back(vEffNBP_RPC);
            
            //Push back the vector with the error on eff of LB to a larger vector of vectors (one element of this = one run) - per RPC
            vErrEffBoth_RPC_runs.push_back(vErrEffBoth_RPC);
            vErrEffBP_RPC_runs.push_back(vErrEffBP_RPC); 
            vErrEffNBP_RPC_runs.push_back(vErrEffNBP_RPC);

            //Clear vector of eff per Plane in a run 
            vEffBoth_Planes.clear();
            vEffBP_Planes.clear();
            vEffNBP_Planes.clear();
            //Clear vector of error on eff per Plane in a run
            vErrEffBoth_Planes.clear();
            vErrEffBP_Planes.clear();
            vErrEffNBP_Planes.clear();

            //Clear vector of eff per RPC in a run 
            vEffBoth_RPC.clear();
            vEffBP_RPC.clear();
            vEffNBP_RPC.clear();
            //Clear vector of error on eff per RPC in a run
            vErrEffBoth_RPC.clear();
            vErrEffBP_RPC.clear();
            vErrEffNBP_RPC.clear();

        } //End of loop on all runs in a given period

        cout << "Size of vStart in period: " << analyzedPeriod.at(per) << ": " << vStart.size() << endl;
        cout << "Size of vRun " << analyzedPeriod.at(per) << ": " << vRun.size() << endl;
        cout << "Size of vEffBoth_Planes " << analyzedPeriod.at(per) << ": " << vEffBoth_Planes_runs.size() << endl;
        cout << "size of B field vector " << analyzedPeriod.at(per) << ": " << vBField.size() << endl;

        //Convert in vectors with efficiency per plane and per run
        vector<double> vEffPerPlaneBoth, vEffPerPlaneBP, vEffPerPlaneNBP;
        vector<double> vErrEffPerPlaneBoth, vErrEffPerPlaneBP, vErrEffPerPlaneNBP;

        vector<vector<double>> vEffPerPlanePerRunBoth, vEffPerPlanePerRunBP, vEffPerPlanePerRunNBP;
        vector<vector<double>> vErrEffPerPlanePerRunBoth, vErrEffPerPlanePerRunBP, vErrEffPerPlanePerRunNBP;
        
        int plane = 0; //Variable to be used in the following loop to keep track of the plane (MT11, MT12, MT21, MT22)
        //----------------------------------------------------------//
        //eff of both planes per plane
        //vEffBoth_Planes_runs = [[11,12,21,22]......[11,12,21,22]]
        //                           run #1             run #n 
        //
        //We want to get efficiency per plane and per run so the following loop creates this
        //
        //vEffPerPlanePerRunBoth = [[11_run#1...11_run#n],[12_run#1...12_run#n],[21_run#1...21_run#n],[22_run#1...22_run#n]]
        //
        //So that we can plot it later on and also the sae for BP and NBP alone

        for (int i = 0; i <= (vEffBoth_Planes_runs.size())*4; i++) {
            
            if ((i % vEffBoth_Planes_runs.size() == 0) && (i != 0)) {
                //cout << "Pushing back vector" << endl;
                //cout << i-(plane*vEffBoth_Planes_runs.size()) << "\t" << plane << endl;
                if (plane > 0) {
                    vEffPerPlaneBoth.push_back(vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())-1][plane]);
                    vEffPerPlaneBP.push_back(vEffBP_Planes_runs[i-(plane*vEffBP_Planes_runs.size())-1][plane]);
                    vEffPerPlaneNBP.push_back(vEffNBP_Planes_runs[i-(plane*vEffNBP_Planes_runs.size())-1][plane]);
                    vErrEffPerPlaneBoth.push_back(vErrEffBoth_Planes_runs[i-(plane*vErrEffBoth_Planes_runs.size())-1][plane]);
                    vErrEffPerPlaneBP.push_back(vErrEffBP_Planes_runs[i-(plane*vErrEffBP_Planes_runs.size())-1][plane]);
                    vErrEffPerPlaneNBP.push_back(vErrEffNBP_Planes_runs[i-(plane*vErrEffNBP_Planes_runs.size())-1][plane]);
                    //cout << "pushing back element " << vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())-1][plane] << endl;
                }
                vEffPerPlanePerRunBoth.push_back(vEffPerPlaneBoth);
                vEffPerPlanePerRunBP.push_back(vEffPerPlaneBP);
                vEffPerPlanePerRunNBP.push_back(vEffPerPlaneNBP);
                vErrEffPerPlanePerRunBoth.push_back(vErrEffPerPlaneBoth);
                vErrEffPerPlanePerRunBP.push_back(vErrEffPerPlaneBP);
                vErrEffPerPlanePerRunNBP.push_back(vErrEffPerPlaneNBP);
                
                vEffPerPlaneBoth.clear();
                vEffPerPlaneBP.clear();
                vEffPerPlaneNBP.clear();
                vErrEffPerPlaneBoth.clear();
                vErrEffPerPlaneBP.clear();
                vErrEffPerPlaneNBP.clear();
                
                plane++; //move to the next plane
            }
            
            else {
                //cout << "pushing back element " << vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())][plane] << "\t i \t" << i << "\t plane \t" << plane << endl;
                //cout << "pushing back element with index " << i-(plane*vEffBoth_Planes_runs.size()) << "\t i \t" << i << "\t plane \t" << plane << endl;
                if (plane == 0) {
                    //Eff per plane on both planes, BP and NBP
                    vEffPerPlaneBoth.push_back(vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())][plane]);
                    vEffPerPlaneBP.push_back(vEffBP_Planes_runs[i-(plane*vEffBP_Planes_runs.size())][plane]);
                    vEffPerPlaneNBP.push_back(vEffNBP_Planes_runs[i-(plane*vEffNBP_Planes_runs.size())][plane]);
                    //Err on eff per plane on both planes, BP and NBP
                    vErrEffPerPlaneBoth.push_back(vErrEffBoth_Planes_runs[i-(plane*vErrEffBoth_Planes_runs.size())][plane]);
                    vErrEffPerPlaneBP.push_back(vErrEffBP_Planes_runs[i-(plane*vErrEffBP_Planes_runs.size())][plane]);
                    vErrEffPerPlaneNBP.push_back(vErrEffNBP_Planes_runs[i-(plane*vErrEffNBP_Planes_runs.size())][plane]);
                }
                else {
                    //Eff per plane on both planes, BP and NBP
                    vEffPerPlaneBoth.push_back(vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())-1][plane]);
                    vEffPerPlaneBP.push_back(vEffBP_Planes_runs[i-(plane*vEffBP_Planes_runs.size())-1][plane]);
                    vEffPerPlaneNBP.push_back(vEffNBP_Planes_runs[i-(plane*vEffNBP_Planes_runs.size())-1][plane]);
                    //Err on eff per plane on both planes, BP and NBP
                    vErrEffPerPlaneBoth.push_back(vErrEffBoth_Planes_runs[i-(plane*vErrEffBoth_Planes_runs.size())-1][plane]);
                    vErrEffPerPlaneBP.push_back(vErrEffBP_Planes_runs[i-(plane*vErrEffBP_Planes_runs.size())-1][plane]);
                    vErrEffPerPlaneNBP.push_back(vErrEffNBP_Planes_runs[i-(plane*vErrEffNBP_Planes_runs.size())-1][plane]);
                }   
            }
        } //End of vector "conversion" for planes

        ////Convert in vectors with efficiency per RPC and per run
        vector<double> vEffPerRPCBoth, vEffPerRPCBP, vEffPerRPCNBP;
        vector<double> vErrEffPerRPCBoth, vErrEffPerRPCBP, vErrEffPerRPCNBP;

        vector<vector<double>> vEffPerRPCPerRunBoth, vEffPerRPCPerRunBP, vEffPerRPCPerRunNBP;
        vector<vector<double>> vErrEffPerRPCPerRunBoth, vErrEffPerRPCPerRunBP, vErrEffPerRPCPerRunNBP;
        
        int rpc = 0; //Variable to be used in the following loop to keep track of the RPC (0 -> 71)
        //----------------------------------------------------------//
        //eff of both planes per plane
        //vEffBoth_RPC_runs = [[0,.....,71]......[0,.....,71]]
        //                       run #1             run #n 
        //
        //We want to get efficiency per plane and per run so the following loop creates this
        //
        //vEffPerPlanePerRunBoth = [[0_run#1...0_run#n],.....[71run#1...71_run#n]]
        //
        //So that we can plot it later on and also the sae for BP and NBP alone

        for (int i = 0; i <= (vEffBoth_RPC_runs.size())*72; i++) {
            
            if ((i % vEffBoth_RPC_runs.size() == 0) && (i != 0)) {
                //cout << "Pushing back vector" << endl;
                //cout << i-(rpc*vEffBoth_RPC_runs.size()) << "\t" << rpc << endl;
                if (rpc > 0) {
                    //Eff per RPC on both planes, BP and NBP
                    vEffPerRPCBoth.push_back(vEffBoth_RPC_runs[i-(rpc*vEffBoth_RPC_runs.size())-1][rpc]);
                    vEffPerRPCBP.push_back(vEffBP_RPC_runs[i-(rpc*vEffBP_RPC_runs.size())-1][rpc]);
                    vEffPerRPCNBP.push_back(vEffNBP_RPC_runs[i-(rpc*vEffNBP_RPC_runs.size())-1][rpc]);
                    //Err on eff per RPC on both planes, BP and NBP
                    vErrEffPerRPCBoth.push_back(vErrEffBoth_RPC_runs[i-(rpc*vErrEffBoth_RPC_runs.size())-1][rpc]);
                    vErrEffPerRPCBP.push_back(vErrEffBP_RPC_runs[i-(rpc*vErrEffBP_RPC_runs.size())-1][rpc]);
                    vErrEffPerRPCNBP.push_back(vErrEffNBP_RPC_runs[i-(rpc*vErrEffNBP_RPC_runs.size())-1][rpc]);
                    //cout << "pushing back element " << vEffBoth_RPC_runs[i-(rpc*vEffBoth_RPC_runs.size())-1][rpc] << endl;
                }
                vEffPerRPCPerRunBoth.push_back(vEffPerRPCBoth);
                vEffPerRPCPerRunBP.push_back(vEffPerRPCBP);
                vEffPerRPCPerRunNBP.push_back(vEffPerRPCNBP);

                vErrEffPerRPCPerRunBoth.push_back(vErrEffPerRPCBoth);
                vErrEffPerRPCPerRunBP.push_back(vErrEffPerRPCBP);
                vErrEffPerRPCPerRunNBP.push_back(vErrEffPerRPCNBP);
                
                vEffPerRPCBoth.clear();
                vEffPerRPCBP.clear();
                vEffPerRPCNBP.clear();
                vErrEffPerRPCBoth.clear();
                vErrEffPerRPCBP.clear();
                vErrEffPerRPCNBP.clear();
                
                rpc++; //move to the next rpc
            }
            
            else {
                //cout << "pushing back element " << vEffBoth_RPC_runs[i-(rpc*vEffBoth_RPC_runs.size())][rpc] << "\t i \t" << i << "\t rpc \t" << rpc << endl;
                //cout << "pushing back element with index " << i-(rpc*vEffBoth_RPC_runs.size()) << "\t i \t" << i << "\t rpc \t" << rpc << endl;
                if (rpc == 0) {
                    //Eff per RPC on both planes, BP and NBP
                    vEffPerRPCBoth.push_back(vEffBoth_RPC_runs[i-(rpc*vEffBoth_RPC_runs.size())][rpc]);
                    vEffPerRPCBP.push_back(vEffBP_RPC_runs[i-(rpc*vEffBP_RPC_runs.size())][rpc]);
                    vEffPerRPCNBP.push_back(vEffNBP_RPC_runs[i-(rpc*vEffNBP_RPC_runs.size())][rpc]);
                    //Err on eff per RPC on both planes, BP and NBP
                    vErrEffPerRPCBoth.push_back(vErrEffBoth_RPC_runs[i-(rpc*vErrEffBoth_RPC_runs.size())][rpc]);
                    vErrEffPerRPCBP.push_back(vErrEffBP_RPC_runs[i-(rpc*vErrEffBP_RPC_runs.size())][rpc]);
                    vErrEffPerRPCNBP.push_back(vErrEffNBP_RPC_runs[i-(rpc*vErrEffNBP_RPC_runs.size())][rpc]);
                }
                else {
                    //Eff per RPC on both planes, BP and NBP
                    vEffPerRPCBoth.push_back(vEffBoth_RPC_runs[i-(rpc*vEffBoth_RPC_runs.size())-1][rpc]);
                    vEffPerRPCBP.push_back(vEffBP_RPC_runs[i-(rpc*vEffBP_RPC_runs.size())-1][rpc]);
                    vEffPerRPCNBP.push_back(vEffNBP_RPC_runs[i-(rpc*vEffNBP_RPC_runs.size())-1][rpc]);
                    //Err on eff per RPC on both planes, BP and NBP
                    vErrEffPerRPCBoth.push_back(vErrEffBoth_RPC_runs[i-(rpc*vErrEffBoth_RPC_runs.size())-1][rpc]);
                    vErrEffPerRPCBP.push_back(vErrEffBP_RPC_runs[i-(rpc*vErrEffBP_RPC_runs.size())-1][rpc]);
                    vErrEffPerRPCNBP.push_back(vErrEffNBP_RPC_runs[i-(rpc*vErrEffNBP_RPC_runs.size())-1][rpc]);
                }   
            }
        } //End of vector "conversion" for RPCs

        for (unsigned int kk = 0; kk < vRun.size(); kk++) {
            cout << "run: " << vRun.at(kk) << " eff both planes MT11: " << vEffPerPlanePerRunBoth[0][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff BP  MT11: " << vEffPerPlanePerRunBP[0][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff NBP MT11: " << vEffPerPlanePerRunNBP[0][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff both planes MT12: " << vEffPerPlanePerRunBoth[1][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff BP  MT12: " << vEffPerPlanePerRunBP[1][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff NBP MT12: " << vEffPerPlanePerRunNBP[1][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff both planes MT21: " << vEffPerPlanePerRunBoth[2][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff BP  MT21: " << vEffPerPlanePerRunBP[2][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff NBP MT21: " << vEffPerPlanePerRunNBP[2][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff both planes MT22: " << vEffPerPlanePerRunBoth[3][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff BP  MT22: " << vEffPerPlanePerRunBP[3][kk] << endl;
            cout << "run: " << vRun.at(kk) << " eff NBP MT22: " << vEffPerPlanePerRunNBP[3][kk] << endl;
        }

        //Create graphs of efficiency vs run number
        for (int i = 0; i < nBinsPlane; i++) {
            TGraphErrors *g1 = new TGraphErrors(vRun.size(),&vRun[0],&vEffPerPlanePerRunBoth[i][0],NULL,&vErrEffPerPlanePerRunBoth[i][0]);
            TGraphErrors *g2 = new TGraphErrors(vRun.size(),&vRun[0],&vEffPerPlanePerRunBP[i][0],NULL,&vErrEffPerPlanePerRunBP[i][0]);
            TGraphErrors *g3 = new TGraphErrors(vRun.size(),&vRun[0],&vEffPerPlanePerRunNBP[i][0],NULL,&vErrEffPerPlanePerRunNBP[i][0]);

            cout << endl << endl << vRun.size() << "\t" << vEffPerPlanePerRunNBP[i].size() << "\t" << vErrEffPerPlanePerRunNBP[i].size() << endl<<endl;

            if (isTime) {
                g1->GetXaxis()->SetTimeDisplay(1);
                g1->GetXaxis()->SetNdivisions(503);
                g1->GetXaxis()->SetTimeFormat("%Y-%m-%d");
                g1->GetXaxis()->SetTimeOffset(0,"gmt");
                g1->GetXaxis()->SetTitle("Time [UTC]");
            }
            else {
                g1->GetXaxis()->SetTitle("Run #");
                g1->GetXaxis()->SetNoExponent(1);
            }
            g1->SetMarkerSize(1.5);
            g1->SetMarkerColor(kBlack);
            g1->SetMarkerStyle(markerStyle.at(per));
            g1->GetXaxis()->SetTitleOffset(0.9);
            g1->GetXaxis()->SetTitleSize(0.04);
            g1->GetXaxis()->SetTitleFont(62);
            g1->GetXaxis()->SetLabelSize(0.04);
            g1->GetXaxis()->SetLabelFont(62);
            g1->GetYaxis()->SetTitle("Efficiency [%]");
            g1->GetYaxis()->SetTitleOffset(1.05);
            g1->GetYaxis()->SetTitleSize(0.04);
            g1->GetYaxis()->SetTitleFont(62);
            g1->GetYaxis()->SetLabelSize(0.04);
            g1->GetYaxis()->SetLabelFont(62);

            if (isTime) {
                g2->GetXaxis()->SetTimeDisplay(1);
                g2->GetXaxis()->SetNdivisions(503);
                g2->GetXaxis()->SetTimeFormat("%Y-%m-%d");
                g2->GetXaxis()->SetTimeOffset(0,"gmt");
                g2->GetXaxis()->SetTitle("Time [UTC]");
            }
            else {
                g2->GetXaxis()->SetTitle("Run #");
                g2->GetXaxis()->SetNoExponent(1);
            }
            
            g2->SetMarkerSize(1.5);
            g2->SetMarkerColor(kRed);
            g2->SetMarkerStyle(markerStyle.at(per));
            g2->GetXaxis()->SetTitleOffset(0.9);
            g2->GetXaxis()->SetTitleSize(0.04);
            g2->GetXaxis()->SetTitleFont(62);
            g2->GetXaxis()->SetLabelSize(0.04);
            g2->GetXaxis()->SetLabelFont(62);
            g2->GetYaxis()->SetTitle("Efficiency [%]");
            g2->GetYaxis()->SetTitleOffset(1.05);
            g2->GetYaxis()->SetTitleSize(0.04);
            g2->GetYaxis()->SetTitleFont(62);
            g2->GetYaxis()->SetLabelSize(0.04);
            g2->GetYaxis()->SetLabelFont(62);

            if (isTime) {
                g3->GetXaxis()->SetTimeDisplay(1);
                g3->GetXaxis()->SetNdivisions(503);
                g3->GetXaxis()->SetTimeFormat("%Y-%m-%d");
                g3->GetXaxis()->SetTimeOffset(0,"gmt");
                g3->GetXaxis()->SetTitle("Time [UTC]");
            }
            else {
                g3->GetXaxis()->SetTitle("Run #");
                g3->GetXaxis()->SetNoExponent(1);
            }
            g3->SetMarkerSize(1.5);
            g3->SetMarkerColor(kGreen+3);
            g3->SetMarkerStyle(markerStyle.at(per));
            g3->GetXaxis()->SetTitleOffset(0.9);
            g3->GetXaxis()->SetTitleSize(0.04);
            g3->GetXaxis()->SetTitleFont(62);
            g3->GetXaxis()->SetLabelSize(0.04);
            g3->GetXaxis()->SetLabelFont(62);
            g3->GetYaxis()->SetTitle("Efficiency [%]");
            g3->GetYaxis()->SetTitleOffset(1.05);
            g3->GetYaxis()->SetTitleSize(0.04);
            g3->GetYaxis()->SetTitleFont(62);
            g3->GetYaxis()->SetLabelSize(0.04);
            g3->GetYaxis()->SetLabelFont(62);

            if (i == 0) {
                lAllAnalyzedPeriodsBoth->AddEntry(g1,(analyzedPeriod.at(per)).c_str(),"p");
                lAllAnalyzedPeriodsBP->AddEntry(g2,(analyzedPeriod.at(per)).c_str(),"p");
                lAllAnalyzedPeriodsNBP->AddEntry(g3,(analyzedPeriod.at(per)).c_str(),"p");
            }
            
            mEffPlaneBoth.at(i)->Add(g1); 
            mEffPlaneBP.at(i)->Add(g2);
            mEffPlaneNBP.at(i)->Add(g3);
        } 

        vRun.clear();
        vIR.clear();
        vBField.clear();
        vStart.clear();
        vEnd.clear();

        vEffBoth_Planes_runs.clear();;
        vEffBP_Planes_runs.clear();        
        vEffNBP_Planes_runs.clear();
        
        vErrEffBoth_Planes_runs.clear();
        vErrEffBP_Planes_runs.clear(); 
        vErrEffNBP_Planes_runs.clear();

        for (unsigned int i = 0; i < nBinsPlane; i++) {
            vEffPerPlanePerRunBoth[i].clear();
            vEffPerPlanePerRunBP[i].clear();
            vEffPerPlanePerRunNBP[i].clear();
            vErrEffPerPlanePerRunBoth[i].clear();
            vErrEffPerPlanePerRunBP[i].clear();
            vErrEffPerPlanePerRunNBP[i].clear();
        }

    } //End of loop on all periods
    
    //cout << vEffPerPlanePerRunBoth_PbPb[0].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[1].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[2].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[3].size() << endl;
    //cout << vEffPerPlanePerRunBP_PbPb[0].size() << "\t" << vEffPerPlanePerRunBP_PbPb[1].size() << "\t" << vEffPerPlanePerRunBP_PbPb[2].size() << "\t" << vEffPerPlanePerRunBP_PbPb[3].size() << endl;
    //cout << vEffPerPlanePerRunNBP_PbPb[0].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[1].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[2].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[3].size() << endl;

    //Both
    //TCanvas *cEffPerPlaneBoth = new TCanvas();
    //cEffPerPlaneBoth->Divide(1,4);
    TCanvas *cEffPerPlaneBoth[4];
    TCanvas *cEffPerPlaneBP[4];
    TCanvas *cEffPerPlaneNBP[4];
    double x,y;
    TMarker *m;

    for (int i = 0; i < nBinsPlane; i++) {
        //cEffPerPlaneBoth->cd(i+1);
        cEffPerPlaneBoth[i] = new TCanvas();
        cEffPerPlaneBP[i] = new TCanvas();
        cEffPerPlaneNBP[i] = new TCanvas();
        
        cEffPerPlaneBoth[i]->cd();
        cEffPerPlaneBoth[i]->SetGridx();
        cEffPerPlaneBoth[i]->SetGridy();
        mEffPlaneBoth.at(i)->SetTitle((planeName[i] + " Both planes efficiency").c_str());
        mEffPlaneBoth.at(i)->GetXaxis()->SetTitle("Run #");
        mEffPlaneBoth.at(i)->GetXaxis()->SetTitleOffset(0.9);
        mEffPlaneBoth.at(i)->GetXaxis()->SetTitleSize(0.04);
        mEffPlaneBoth.at(i)->GetXaxis()->SetTitleFont(62);
        mEffPlaneBoth.at(i)->GetXaxis()->SetLabelSize(0.04);
        mEffPlaneBoth.at(i)->GetXaxis()->SetLabelFont(62);
        mEffPlaneBoth.at(i)->GetXaxis()->SetNoExponent(1);
        mEffPlaneBoth.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        mEffPlaneBoth.at(i)->GetYaxis()->SetTitleOffset(0.85);
        mEffPlaneBoth.at(i)->GetYaxis()->SetTitleSize(0.04);
        mEffPlaneBoth.at(i)->GetYaxis()->SetTitleFont(62);
        mEffPlaneBoth.at(i)->GetYaxis()->SetLabelSize(0.04);
        mEffPlaneBoth.at(i)->GetYaxis()->SetLabelFont(62); 
        mEffPlaneBoth.at(i)->GetYaxis()->SetRangeUser(80,100); 
        mEffPlaneBoth.at(i)->Draw("AP");
        TList *l1 = mEffPlaneBoth.at(i)->GetListOfGraphs();
        TIter next(l1);
        TObject *obj1;
        
        /*while ((obj1 = next())) {
            TGraphErrors *g = (TGraphErrors*)obj1;
            for (int j = 0; j < g->GetN(); j++) {
                g->GetPoint(j,x,y);
                m = new TMarker(x,y,g->GetMarkerStyle());
                if (vBFieldTot.at(j) == "minus") {
                    //m->SetMarkerColor(kOrange);
                    m->SetMarkerColor(kBlack);
                    m->SetMarkerSize(1);
                    m->Draw("SAME");
                }
            }
        }*/
        
        lAllAnalyzedPeriodsBoth->Draw("SAME");

        cEffPerPlaneBP[i]->cd();
        cEffPerPlaneBP[i]->SetGridx();
        cEffPerPlaneBP[i]->SetGridy();
        mEffPlaneBP.at(i)->SetTitle((planeName[i] + " BP efficiency").c_str());
        mEffPlaneBP.at(i)->GetXaxis()->SetTitle("Run #");
        mEffPlaneBP.at(i)->GetXaxis()->SetTitleOffset(0.9);
        mEffPlaneBP.at(i)->GetXaxis()->SetTitleSize(0.04);
        mEffPlaneBP.at(i)->GetXaxis()->SetTitleFont(62);
        mEffPlaneBP.at(i)->GetXaxis()->SetLabelSize(0.04);
        mEffPlaneBP.at(i)->GetXaxis()->SetLabelFont(62);
        mEffPlaneBP.at(i)->GetXaxis()->SetNoExponent(1);
        mEffPlaneBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        mEffPlaneBP.at(i)->GetYaxis()->SetTitleOffset(0.85);
        mEffPlaneBP.at(i)->GetYaxis()->SetTitleSize(0.04);
        mEffPlaneBP.at(i)->GetYaxis()->SetTitleFont(62);
        mEffPlaneBP.at(i)->GetYaxis()->SetLabelSize(0.04);
        mEffPlaneBP.at(i)->GetYaxis()->SetLabelFont(62); 
        mEffPlaneBP.at(i)->GetYaxis()->SetRangeUser(80,100);
        mEffPlaneBP.at(i)->Draw("AP");
        TList *l2 = mEffPlaneBP.at(i)->GetListOfGraphs();
        TIter next2(l2);
        TObject *obj2;
        
        /*while ((obj2 = next2())) {
            TGraphErrors *g = (TGraphErrors*)obj2;
            for (int j = 0; j < g->GetN(); j++) {
                g->GetPoint(j,x,y);
                m = new TMarker(x,y,g->GetMarkerStyle());
                if (vBFieldTot.at(j) == "minus") {
                    //m->SetMarkerColor(kOrange);
                    m->SetMarkerColor(kRed);
                    m->SetMarkerSize(1);
                    m->Draw("SAME");
                }
            }
        }*/

        lAllAnalyzedPeriodsBP->Draw("SAME");

        cEffPerPlaneNBP[i]->cd();
        cEffPerPlaneNBP[i]->SetGridx();
        cEffPerPlaneNBP[i]->SetGridy();
        mEffPlaneNBP.at(i)->SetTitle((planeName[i] + " NBP efficiency").c_str());
        mEffPlaneNBP.at(i)->GetXaxis()->SetTitle("Run #");
        mEffPlaneNBP.at(i)->GetXaxis()->SetTitleOffset(0.9);
        mEffPlaneNBP.at(i)->GetXaxis()->SetTitleSize(0.04);
        mEffPlaneNBP.at(i)->GetXaxis()->SetTitleFont(62);
        mEffPlaneNBP.at(i)->GetXaxis()->SetLabelSize(0.04);
        mEffPlaneNBP.at(i)->GetXaxis()->SetLabelFont(62);
        mEffPlaneNBP.at(i)->GetXaxis()->SetNoExponent(1);
        mEffPlaneNBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        mEffPlaneNBP.at(i)->GetYaxis()->SetTitleOffset(0.85);
        mEffPlaneNBP.at(i)->GetYaxis()->SetTitleSize(0.04);
        mEffPlaneNBP.at(i)->GetYaxis()->SetTitleFont(62);
        mEffPlaneNBP.at(i)->GetYaxis()->SetLabelSize(0.04);
        mEffPlaneNBP.at(i)->GetYaxis()->SetLabelFont(62);
        mEffPlaneNBP.at(i)->GetYaxis()->SetRangeUser(80,100); 
        mEffPlaneNBP.at(i)->Draw("AP"); 
        TList *l3 = mEffPlaneBP.at(i)->GetListOfGraphs();
        TIter next3(l3);
        TObject *obj3;
        
        /*while ((obj3 = next3())) {
            TGraphErrors *g = (TGraphErrors*)obj3;
            for (int j = 0; j < g->GetN(); j++) {
                g->GetPoint(j,x,y);
                m = new TMarker(x,y,g->GetMarkerStyle());
                if (vBFieldTot.at(j) == "minus") {
                    m->SetMarkerColor(kOrange);
                    cout << "here" << endl;
                    //m->SetMarkerColor(kGreen+3);
                    m->SetMarkerSize(1);
                    m->Draw("SAME");
                }
            }
        }*/

        lAllAnalyzedPeriodsNBP->Draw("SAME");
    }
    
    //BP
    //TCanvas *cEffPerPlaneBP = new TCanvas();
    //cEffPerPlaneBP->Divide(1,4);
    /*TCanvas *cEffPerPlaneBP[4];
    
    for (int i = 0; i < nBinsPlane; i++) {
        //cEffPerPlaneBP->cd(i+1);
        cEffPerPlaneBP[i] = new TCanvas();
        cEffPerPlaneBP[i]->cd();

        gEffPerRunPlaneBP.at(i)->SetMarkerStyle(8);
        gEffPerRunPlaneBP.at(i)->SetMarkerSize(1);
        gEffPerRunPlaneBP.at(i)->SetMarkerColor(kRed);
        gEffPerRunPlaneBP.at(i)->SetTitle((planeName[i] + " BP").c_str());
        if (isTime) {
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTimeDisplay(1);
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetNdivisions(503);
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTimeFormat("%Y-%m-%d");
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTimeOffset(0,"gmt");
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitle("Time [UTC]");
        }
        else {
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitle("Run #");
            gEffPerRunPlaneBP.at(i)->GetXaxis()->SetNoExponent(1);
        }
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitleOffset(0.9);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitleSize(0.04);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetLabelSize(0.04);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitleOffset(1.05);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitleSize(0.04);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetLabelSize(0.04);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunPlaneBP.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunPlaneBP.at(i)->Draw("AP");

        bool changedMarkerPlus = false, changedMarkerMinus = false; //this bool is used because the first run that is analyzed is the one from the data reconstructed on grid and then the
        //same one from AO2D (same run number) and in this way I can change the marker only for the run that I have reconstructed myself 

        for (int j = 0; j < gEffPerRunPlaneBP.at(i)->GetN(); j++) {
            gEffPerRunPlaneBP.at(i)->GetPoint(j,x,y);

            if (vBField_2023.at(j) == "plus") {
                if ((int)x == 544565 && changedMarkerPlus == false) {
                    //cout << "here!" << endl;
                    m = new TMarker(x,y,33);
                    m->SetMarkerColor(kRed);
                    m->SetMarkerSize(2);  
                    changedMarkerPlus = true;
                }
                else {
                    m = new TMarker(x,y,8);
                    m->SetMarkerColor(kRed);
                    m->SetMarkerSize(1);
                }
                m->Draw("SAME");
            }
            
            else if (vBField_2023.at(j) == "minus") {
                if ((int)x == 545185 && changedMarkerMinus == false) {
                    m = new TMarker(x,y,33);
                    m->SetMarkerColor(kOrange);
                    m->SetMarkerSize(2);
                    changedMarkerMinus = true;
                }
                else {
                    m = new TMarker(x,y,8);
                    m->SetMarkerColor(kOrange);
                    m->SetMarkerSize(1);
                }               
                m->Draw("SAME");
            }
        }

        TLegend *lEffPlaneBP = (TLegend*)lEffPlenBPIR->Clone();
        lEffPlaneBP->Draw("SAME");
        gPad->Update();
        if (i != 3) {
            lEffPlaneBP->SetX1NDC(0.799);
            lEffPlaneBP->SetY1NDC(0.605);
            lEffPlaneBP->SetX2NDC(0.876);
            lEffPlaneBP->SetY2NDC(0.859);
        }
        else {
            lEffPlaneBP->SetX1NDC(0.799);
            lEffPlaneBP->SetY1NDC(0.250);
            lEffPlaneBP->SetX2NDC(0.876);
            lEffPlaneBP->SetY2NDC(0.504);
        }
        gPad->Modified();
    }
    
    //NBP
    //TCanvas *cEffPerPlaneNBP = new TCanvas();
    //cEffPerPlaneNBP->Divide(1,4);
    TCanvas *cEffPerPlaneNBP[4];
    
    for (int i = 0; i < nBinsPlane; i++) {
        //cEffPerPlaneNBP->cd(i+1);
        cEffPerPlaneNBP[i] = new TCanvas();
        cEffPerPlaneNBP[i]->cd();

        gEffPerRunPlaneNBP.at(i)->SetMarkerStyle(8);
        gEffPerRunPlaneNBP.at(i)->SetMarkerSize(1);
        gEffPerRunPlaneNBP.at(i)->SetMarkerColor(kGreen+3);
        gEffPerRunPlaneNBP.at(i)->SetTitle((planeName[i] + " NBP").c_str());
        if (isTime) {
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTimeDisplay(1);
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetNdivisions(503);
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTimeFormat("%Y-%m-%d");
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTimeOffset(0,"gmt");
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitle("Time [UTC]");
        }
        else {
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitle("Run #");
            gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetNoExponent(1);
        }
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitleOffset(0.9);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitleSize(0.04);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetLabelSize(0.04);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitleOffset(1.05);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitleSize(0.04);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetLabelSize(0.04);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunPlaneNBP.at(i)->Draw("AP");

        bool changedMarkerPlus = false, changedMarkerMinus = false; //this bool is used because the first run that is analyzed is the one from the data reconstructed on grid and then the
        //same one from AO2D (same run number) and in this way I can change the marker only for the run that I have reconstructed myself 

        for (int j = 0; j < gEffPerRunPlaneNBP.at(i)->GetN(); j++) {
            gEffPerRunPlaneNBP.at(i)->GetPoint(j,x,y);

            if (vBField_2023.at(j) == "plus") {
                if ((int)x == 544565 && changedMarkerPlus == false) {
                    //cout << "here!" << endl;
                    m = new TMarker(x,y,33);
                    m->SetMarkerColor(kGreen+3);
                    m->SetMarkerSize(2);  
                    changedMarkerPlus = true;
                }
                else {
                    m = new TMarker(x,y,8);
                    m->SetMarkerColor(kGreen+3);
                    m->SetMarkerSize(1);
                }
                m->Draw("SAME");
            }
            
            else if (vBField_2023.at(j) == "minus") {
                if ((int)x == 545185 && changedMarkerMinus == false) {
                    m = new TMarker(x,y,33);
                    m->SetMarkerColor(kOrange);
                    m->SetMarkerSize(2);
                    changedMarkerMinus = true;
                }
                else {
                    m = new TMarker(x,y,8);
                    m->SetMarkerColor(kOrange);
                    m->SetMarkerSize(1);
                }               
                m->Draw("SAME");
            }
        }

        TLegend *lEffPlaneNBP = (TLegend*)lEffPlenNBPIR->Clone();
        lEffPlaneNBP->Draw("SAME");
        gPad->Update();
        if (i != 3) {
            lEffPlaneNBP->SetX1NDC(0.799);
            lEffPlaneNBP->SetY1NDC(0.605);
            lEffPlaneNBP->SetX2NDC(0.876);
            lEffPlaneNBP->SetY2NDC(0.859);
        }
        else {
            lEffPlaneNBP->SetX1NDC(0.799);
            lEffPlaneNBP->SetY1NDC(0.250);
            lEffPlaneNBP->SetX2NDC(0.876);
            lEffPlaneNBP->SetY2NDC(0.504);
        }
        gPad->Modified();
    }

    

    //Vector of TGraphErrors for all 72 RPCs on both, BP and NBP
    vector<TGraphErrors*> gEffPerRunRPCBoth, gEffPerRunRPCBP, gEffPerRunRPCNBP;
    
    for (int i = 0; i < nBinsRPC; i++) {
        //cout << i << "\t" << vEffPerRPCPerRunBoth[i].size() << endl;
        TGraphErrors *g1 = new TGraphErrors(vRun_2023.size(),&vRun_2023[0],&vEffPerRPCPerRunBoth[i][0],NULL,&vErrEffPerRPCPerRunBoth[i][0]);
        TGraphErrors *g2 = new TGraphErrors(vRun_2023.size(),&vRun_2023[0],&vEffPerRPCPerRunBP[i][0],NULL,&vErrEffPerRPCPerRunBP[i][0]);
        TGraphErrors *g3 = new TGraphErrors(vRun_2023.size(),&vRun_2023[0],&vEffPerRPCPerRunNBP[i][0],NULL,&vErrEffPerRPCPerRunNBP[i][0]);
        
        gEffPerRunRPCBoth.push_back(g1);
        gEffPerRunRPCBP.push_back(g2);
        gEffPerRunRPCNBP.push_back(g3);
    } 

    bool isTimeRPC = false; //if true -> plots have time on x axis otherwise run number

    //Both
    TCanvas *cEffPerRPCBoth_MT11 = new TCanvas();
    cEffPerRPCBoth_MT11->Divide(2,9);

    TCanvas *cEffPerRPCBoth_MT12 = new TCanvas();
    cEffPerRPCBoth_MT12->Divide(2,9);

    TCanvas *cEffPerRPCBoth_MT21 = new TCanvas();
    cEffPerRPCBoth_MT21->Divide(2,9);

    TCanvas *cEffPerRPCBoth_MT22 = new TCanvas();
    cEffPerRPCBoth_MT22->Divide(2,9);

    gStyle->SetTitleFontSize(0.07);
    
    for (int i = 0; i < nBinsRPC; i++) {        
        if (i >= 0 && i <= 8) { //MT11 ok
            cEffPerRPCBoth_MT11->cd(i+18-(3*i));
        }
        else if (i >= 36 && i <= 44) { //MT11 ok
            cEffPerRPCBoth_MT11->cd(i-19-(3*(i-36)));
        }
        //---//
        else if (i >= 9 && i <= 17){ //MT12 ok
            cEffPerRPCBoth_MT12->cd(i+9-(3*(i-9)));
        }
        else if (i >= 45 && i <= 53){ //MT12 ok
            cEffPerRPCBoth_MT12->cd(i-28-(3*(i-45)));
        }
        //---//
        else if (i >= 18 && i <= 26) { //MT21 ok
            cEffPerRPCBoth_MT21->cd(i-(3*(i-18)));
        }
        else if (i >= 54 && i <= 62) { //MT21 ok
            cEffPerRPCBoth_MT21->cd(i-37-(3*(i-54)));
        }
        //---//
        else if (i >= 27 && i <= 35) { //MT22 ok
            cEffPerRPCBoth_MT22->cd(i-9-(3*(i-27)));
        }
        else if (i >= 63 && i <= 71) { //MT22 ok
            cEffPerRPCBoth_MT22->cd(i-46-(3*(i-63)));
        }
        //---//
        gEffPerRunRPCBoth.at(i)->SetMarkerStyle(8);
        gEffPerRunRPCBoth.at(i)->SetMarkerSize(0.7);
        gEffPerRunRPCBoth.at(i)->SetMarkerColor(kBlack);
        string detName = o2::mid::detparams::getDEName(i);
        gEffPerRunRPCBoth.at(i)->SetTitle(("RPC" + to_string(i) + " Both planes -> " + detName).c_str());
        if (isTimeRPC) {
            gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTimeDisplay(1);
            gEffPerRunRPCBoth.at(i)->GetXaxis()->SetNdivisions(503);
            gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTimeFormat("%Y-%m-%d");
            gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTimeOffset(0,"gmt");
            gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTitle("Time [UTC]");
        }
        else {
            gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTitle("Run #");
        }
        gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPerRunRPCBoth.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunRPCBoth.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPerRunRPCBoth.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunRPCBoth.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunRPCBoth.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPerRunRPCBoth.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPerRunRPCBoth.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunRPCBoth.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPerRunRPCBoth.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunRPCBoth.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunRPCBoth.at(i)->Draw("AP");   
    }
    
    //BP
    TCanvas *cEffPerRPCBP_MT11 = new TCanvas();
    cEffPerRPCBP_MT11->Divide(2,9);

    TCanvas *cEffPerRPCBP_MT12 = new TCanvas();
    cEffPerRPCBP_MT12->Divide(2,9);

    TCanvas *cEffPerRPCBP_MT21 = new TCanvas();
    cEffPerRPCBP_MT21->Divide(2,9);

    TCanvas *cEffPerRPCBP_MT22 = new TCanvas();
    cEffPerRPCBP_MT22->Divide(2,9);

    for (int i = 0; i < nBinsRPC; i++) {
        if (i >= 0 && i <= 8) { //MT11 ok
            cEffPerRPCBP_MT11->cd(i+18-(3*i));
        }
        else if (i >= 36 && i <= 44) { //MT11 ok
            cEffPerRPCBP_MT11->cd(i-19-(3*(i-36)));
        }
        //---//
        else if (i >= 9 && i <= 17){ //MT12 ok
            cEffPerRPCBP_MT12->cd(i+9-(3*(i-9)));
        }
        else if (i >= 45 && i <= 53){ //MT12 ok
            cEffPerRPCBP_MT12->cd(i-28-(3*(i-45)));
        }
        //---//
        else if (i >= 18 && i <= 26) { //MT21 ok
            cEffPerRPCBP_MT21->cd(i-(3*(i-18)));
        }
        else if (i >= 54 && i <= 62) { //MT21 ok
            cEffPerRPCBP_MT21->cd(i-37-(3*(i-54)));
        }
        //---//
        else if (i >= 27 && i <= 35) { //MT22 ok
            cEffPerRPCBP_MT22->cd(i-9-(3*(i-27)));
        }
        else if (i >= 63 && i <= 71) { //MT22 ok
            cEffPerRPCBP_MT22->cd(i-46-(3*(i-63)));
        }
        
        gEffPerRunRPCBP.at(i)->SetMarkerStyle(8);
        gEffPerRunRPCBP.at(i)->SetMarkerSize(0.7);
        gEffPerRunRPCBP.at(i)->SetMarkerColor(kRed);
        string detName = o2::mid::detparams::getDEName(i);
        gEffPerRunRPCBP.at(i)->SetTitle(("RPC" + to_string(i) + " BP -> " + detName).c_str());
        if (isTimeRPC) {
            gEffPerRunRPCBP.at(i)->GetXaxis()->SetTimeDisplay(1);
            gEffPerRunRPCBP.at(i)->GetXaxis()->SetNdivisions(503);
            gEffPerRunRPCBP.at(i)->GetXaxis()->SetTimeFormat("%Y-%m-%d");
            gEffPerRunRPCBP.at(i)->GetXaxis()->SetTimeOffset(0,"gmt");
            gEffPerRunRPCBP.at(i)->GetXaxis()->SetTitle("Time [UTC]");
        }
        else {
            gEffPerRunRPCBP.at(i)->GetXaxis()->SetTitle("Run #");
        }
        gEffPerRunRPCBP.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPerRunRPCBP.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPerRunRPCBP.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunRPCBP.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPerRunRPCBP.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunRPCBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunRPCBP.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPerRunRPCBP.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPerRunRPCBP.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunRPCBP.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPerRunRPCBP.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunRPCBP.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunRPCBP.at(i)->Draw("AP");
    }
    
    //NBP
    TCanvas *cEffPerRPCNBP_MT11 = new TCanvas();
    cEffPerRPCNBP_MT11->Divide(2,9);

    TCanvas *cEffPerRPCNBP_MT12 = new TCanvas();
    cEffPerRPCNBP_MT12->Divide(2,9);

    TCanvas *cEffPerRPCNBP_MT21 = new TCanvas();
    cEffPerRPCNBP_MT21->Divide(2,9);

    TCanvas *cEffPerRPCNBP_MT22 = new TCanvas();
    cEffPerRPCNBP_MT22->Divide(2,9);

    for (int i = 0; i < nBinsRPC; i++) {
        if (i >= 0 && i <= 8) { //MT11 ok
            cEffPerRPCNBP_MT11->cd(i+18-(3*i));
        }
        else if (i >= 36 && i <= 44) { //MT11 ok
            cEffPerRPCNBP_MT11->cd(i-19-(3*(i-36)));
        }
        //---//
        else if (i >= 9 && i <= 17){ //MT12 ok
            cEffPerRPCNBP_MT12->cd(i+9-(3*(i-9)));
        }
        else if (i >= 45 && i <= 53){ //MT12 ok
            cEffPerRPCNBP_MT12->cd(i-28-(3*(i-45)));
        }
        //---//
        else if (i >= 18 && i <= 26) { //MT21 ok
            cEffPerRPCNBP_MT21->cd(i-(3*(i-18)));
        }
        else if (i >= 54 && i <= 62) { //MT21 ok
            cEffPerRPCNBP_MT21->cd(i-37-(3*(i-54)));
        }
        //---//
        else if (i >= 27 && i <= 35) { //MT22 ok
            cEffPerRPCNBP_MT22->cd(i-9-(3*(i-27)));
        }
        else if (i >= 63 && i <= 71) { //MT22 ok
            cEffPerRPCNBP_MT22->cd(i-46-(3*(i-63)));
        }
        
        gEffPerRunRPCNBP.at(i)->SetMarkerStyle(8);
        gEffPerRunRPCNBP.at(i)->SetMarkerSize(0.7);
        gEffPerRunRPCNBP.at(i)->SetMarkerColor(kGreen);
        string detName = o2::mid::detparams::getDEName(i);
        gEffPerRunRPCNBP.at(i)->SetTitle(("RPC" + to_string(i) + " NBP -> " + detName).c_str());
        if (isTimeRPC) {
            gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTimeDisplay(1);
            gEffPerRunRPCNBP.at(i)->GetXaxis()->SetNdivisions(503);
            gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTimeFormat("%Y-%m-%d");
            gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTimeOffset(0,"gmt");
            gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTitle("Time [UTC]");
        }
        else {
            gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTitle("Run #");
        }
        gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPerRunRPCNBP.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunRPCNBP.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPerRunRPCNBP.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunRPCNBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunRPCNBP.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPerRunRPCNBP.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPerRunRPCNBP.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunRPCNBP.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPerRunRPCNBP.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunRPCNBP.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunRPCNBP.at(i)->Draw("AP");
    }

    //hRun_pp2023.close();
    //hDate_pp2023.close();
    //hRun_PbPb2023.close();
    //hDate_PbPb2023.close();
    hIR_pp2023.close();
    hIR_PbPb2023.close();

    bool saveData = false;
    
    if (saveData) {
        fOutEffPlane->cd();
        cEffPlaneBothPlanesIR->Write("effPerPlaneBothPlanesVsIR");
        cEffPlaneBPIR->Write("effPerPlaneBPVsIR");
        cEffPlaneNBPIR->Write("effPerPlaneNBPVsIR");
        for (int i = 0; i < nBinsPlane; i++) {
            cEffPerPlaneBoth[i]->Write(("effPerPlaneBothPlanesVsRun_"+planeName[i]).c_str());
            cEffPerPlaneBP[i]->Write(("effPerPlaneBPVsRun_"+planeName[i]).c_str());
            cEffPerPlaneNBP[i]->Write(("effPerPlaneBPVsRun_"+planeName[i]).c_str());
        }
        fOutEffPlane->Close();
    }*/
}