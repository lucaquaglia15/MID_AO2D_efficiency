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
#include "MIDBase/DetectorParameters.h" //Detector parameter
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

bool debug = false;
const int trackGoal = 8e+6; //Normal value
//long int trackGoal = 8e+11;//To test, it goes faster since this number of tracks is never reached
long int tracks = 0;
long int cumulativeTracks = 0;
int nBinsPlane = 4; //Number of planes
int nBinsRPC = 72; //Number of RPCs
int nBinsBoard = 936; //Number of LBs

o2::ccdb::CcdbApi api; //CCDB API as global object
o2::mid::Mapping mapping; //MID mapping object to construct ccdb object

vector<int> markerStyle{50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,70,71,72,73,74,75};
vector<int> markerColor{1,2,3,4,5,6,7,8,9};

void effByTime() { //Main function

    //Per plane
    double effBothPlane = 0, effBPPlane = 0, effNBPPlane =0;
    double errEffBothPlane = 0, errEffBPPlane = 0, errEffNBPPlane = 0;

    //Per RPC
    double effBothRPC = 0, effBPRPC = 0, effNBPPRPC =0;
    double errEffBothRPC = 0, errEffBPRPC = 0, errEffNBPRPC = 0;

    //Run by run - Plane
    vector<double> vEffBoth_Planes, vEffBP_Planes, vEffNBP_Planes;
    vector<double> vErrEffBoth_Planes, vErrEffBP_Planes, vErrEffNBP_Planes;
    //Pb-Pb
    vector<double> vEffBoth_Planes_PbPb, vEffBP_Planes_PbPb, vEffNBP_Planes_PbPb;
    vector<double> vErrEffBoth_Planes_PbPb, vErrEffBP_Planes_PbPb, vErrEffNBP_Planes_PbPb;

    vector<vector<double>> vEffBoth_Planes_runs, vEffBP_Planes_runs, vEffNBP_Planes_runs;
    vector<vector<double>> vErrEffBoth_Planes_runs, vErrEffBP_Planes_runs, vErrEffNBP_Planes_runs;
    //Pb-Pb
    vector<vector<double>> vEffBoth_Planes_runs_PbPb, vEffBP_Planes_runs_PbPb, vEffNBP_Planes_runs_PbPb;
    vector<vector<double>> vErrEffBoth_Planes_runs_PbPb, vErrEffBP_Planes_runs_PbPb, vErrEffNBP_Planes_runs_PbPb;
    

    //Run by run - RPC
    vector<double> vEffBoth_RPC, vEffBP_RPC, vEffNBP_RPC;
    vector<double> vErrEffBoth_RPC, vErrEffBP_RPC, vErrEffNBP_RPC;
    
    vector<vector<double>> vEffBoth_RPC_runs, vEffBP_RPC_runs, vEffNBP_RPC_runs;
    vector<vector<double>> vErrEffBoth_RPC_runs, vErrEffBP_RPC_runs, vErrEffNBP_RPC_runs;
    
    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    string period_pp2023 = "LHC23_pass4_skimmed"; //pp skimmed QC data of 2023 pass 4
    string period_PbPb2023 = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    
    //pp
    string globalPath_pp2023 = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period_pp2023+"/";
    //PbPb
    string globalPath_PbPb2023 = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period_PbPb2023+"/";

    //Path of the merged file, run-by-run
    //pp
    string runPath_pp2023 = globalPath_pp2023+"runs/";
    //PbPb
    string runPath_PbPb2023 = globalPath_PbPb2023+"runs/";

    //Path for the .txt file of the run list of the period
    //pp
    string runNumbers_pp2023 = globalPath_pp2023+"run_list.txt"; 
    string runDates_pp2023 = globalPath_pp2023+"run_dates.txt";
    //PbPb
    string runNumbers_PbPb2023 = globalPath_PbPb2023+"run_list.txt"; 
    string runDates_PbPb2023 = globalPath_PbPb2023+"run_dates.txt";
    //IR PbPb
    string fIR_PbPb2023 = globalPath_PbPb2023+"run_IR.txt"; 

    //Open txt file of runs
    //pp
    ifstream hRun_pp2023;
    hRun_pp2023.open(runNumbers_pp2023.c_str());
    //PbPb
    ifstream hRun_PbPb2023;
    hRun_PbPb2023.open(runNumbers_PbPb2023.c_str());

    //Open txt file of start/end dates of the runs
    //pp
    ifstream hDate_pp2023;
    hDate_pp2023.open(runDates_pp2023.c_str());
    //PbPb
    ifstream hDate_PbPb2023;
    hDate_PbPb2023.open(runDates_PbPb2023.c_str());

    //Open txt file for IR PbPb 2023
    ifstream hIR_PbPb2023;
    hIR_PbPb2023.open(fIR_PbPb2023.c_str());
    
    //Get start and end of each run
    //pp
    long int runForDate_pp2023;
    double start_pp2023, end_pp2023;
    vector<long int> vRunForDate_pp2023;
    vector<double> vStart_pp2023, vEnd_pp2023;
    vector<double> vStart_2023;
    
    while (hDate_pp2023 >> runForDate_pp2023 >> start_pp2023 >> end_pp2023){
        vRunForDate_pp2023.push_back(runForDate_pp2023);
        vStart_pp2023.push_back(start_pp2023);
        vEnd_pp2023.push_back(end_pp2023);
        vStart_2023.push_back(start_pp2023);
        cout << start_pp2023 << "\t" << start_pp2023/1000 << endl;
        printf("start_pp2023: %f \t start_pp2023/1000: %f",start_pp2023,start_pp2023/1000);
    }
    //PbPb
    long int runForDate_PbPb2023;
    double start_PbPb2023, end_PbPb2023;
    vector<long int> vRunForDate_PbPb2023;
    vector<double> vStart_PbPb2023, vEnd_PbPb2023;
    
    while (hDate_PbPb2023 >> runForDate_PbPb2023 >> start_PbPb2023 >> end_PbPb2023){
        vRunForDate_PbPb2023.push_back(runForDate_PbPb2023);
        vStart_PbPb2023.push_back(start_PbPb2023);
        vEnd_PbPb2023.push_back(end_PbPb2023);
        vStart_2023.push_back(start_pp2023);
    }

    //IR only for PbPb
    double run_PbPb2023IR, IR_PbPb2023;
    vector<double> vRun_PbPb2023_IR, vIR_PbPb2023;
    while (hIR_PbPb2023 >> run_PbPb2023IR >> IR_PbPb2023) {
        vRun_PbPb2023_IR.push_back(run_PbPb2023IR);
        vIR_PbPb2023.push_back(IR_PbPb2023/1000);
    }

    //Push back to a vector of int (no need to care about size)
    //pp
    double run_pp2023;
    vector<double> vRun_pp2023;
    vector<double> vRun_2023;

    while(hRun_pp2023 >> run_pp2023) {
        vRun_pp2023.push_back(run_pp2023);
        vRun_2023.push_back(run_pp2023);
    }
    //sort in ascending order
    sort(vRun_pp2023.begin(), vRun_pp2023.end()); //pp 2023
    
    //PbPb
    double run_PbPb2023;
    vector<double> vRun_PbPb2023;

    while(hRun_PbPb2023 >> run_PbPb2023) {
        vRun_PbPb2023.push_back(run_PbPb2023);
        vRun_2023.push_back(run_PbPb2023);
    }
    //sort in ascending order
    sort(vRun_PbPb2023.begin(), vRun_PbPb2023.end()); //PbPb 2023
    sort(vRun_2023.begin(), vRun_2023.end()); //All 2023 (pp + PbPb)

    //General string name
    string fileName = "AnalysisResults.root";

    //Loop on all runs pp
    for (unsigned int iRun = 0; iRun < vRun_pp2023.size(); iRun++) {
        //Enter the folder
        string runFolder = runPath_pp2023+to_string((int)vRun_pp2023.at(iRun));

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

    } //End of loop on all runs pp

    //Loop on all runs PbPb
    for (unsigned int iRun = 0; iRun < vRun_PbPb2023.size(); iRun++) {
        //Enter the folder
        string runFolder = runPath_PbPb2023+to_string((int)vRun_PbPb2023.at(iRun));

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

        //Calculate eff per run at plane level
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
            //Only PbPb
            vEffBoth_Planes_PbPb.push_back(effBothPlane);
            vEffBP_Planes_PbPb.push_back(effBPPlane);
            vEffNBP_Planes_PbPb.push_back(effNBPPlane);
            
            //Fill vector for error on efficiency per plane in the run
            vErrEffBoth_Planes.push_back(errEffBothPlane);
            vErrEffBP_Planes.push_back(errEffBPPlane);
            vErrEffNBP_Planes.push_back(errEffNBPPlane);
            
            //Only PbPb
            vErrEffBoth_Planes_PbPb.push_back(errEffBothPlane);
            vErrEffBP_Planes_PbPb.push_back(errEffBPPlane);
            vErrEffNBP_Planes_PbPb.push_back(errEffNBPPlane);
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
        //Only Pb-Pb
        vEffBoth_Planes_runs_PbPb.push_back(vEffBoth_Planes_PbPb);
        vEffBP_Planes_runs_PbPb.push_back(vEffBP_Planes_PbPb);
        vEffNBP_Planes_runs_PbPb.push_back(vEffNBP_Planes_PbPb);
        
        //Push back the vector with the error on eff of LB to a larger vector of vectors (one element of this = one run) - per plane
        vErrEffBoth_Planes_runs.push_back(vErrEffBoth_Planes);
        vErrEffBP_Planes_runs.push_back(vErrEffBP_Planes); 
        vErrEffNBP_Planes_runs.push_back(vErrEffNBP_Planes);
        //Only Pb-Pb
        vErrEffBoth_Planes_runs_PbPb.push_back(vErrEffBoth_Planes_PbPb);
        vErrEffBP_Planes_runs_PbPb.push_back(vErrEffBP_Planes_PbPb);
        vErrEffNBP_Planes_runs_PbPb.push_back(vErrEffNBP_Planes_PbPb);

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
        vEffBoth_Planes_PbPb.clear();
        vEffBP_Planes_PbPb.clear();
        vEffNBP_Planes_PbPb.clear();
        //Clear vector of error on eff per Plane in a run
        vErrEffBoth_Planes.clear();
        vErrEffBP_Planes.clear();
        vErrEffNBP_Planes.clear();
        vErrEffBoth_Planes_PbPb.clear();
        vErrEffBP_Planes_PbPb.clear();
        vErrEffNBP_Planes_PbPb.clear();

        //Clear vector of eff per RPC in a run 
        vEffBoth_RPC.clear();
        vEffBP_RPC.clear();
        vEffNBP_RPC.clear();
        //Clear vector of error on eff per RPC in a run
        vErrEffBoth_RPC.clear();
        vErrEffBP_RPC.clear();
        vErrEffNBP_RPC.clear();

    } //End of loop on all runs PbPb

    cout << "Size of vStart_pp2023 " << vStart_pp2023.size() << endl;
    cout << "Size of vStart_PbPb2023 " << vStart_PbPb2023.size() << endl;
    cout << "Size of vStart_2023 " << vStart_2023.size() << endl;
    cout << "Size of vEffBoth_Planes_PbPb " << vEffBoth_Planes_runs_PbPb.size() << endl;

    //per plane
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
    } 

    //per plane PbPb
    vector<double> vEffPerPlaneBoth_PbPb, vEffPerPlaneBP_PbPb, vEffPerPlaneNBP_PbPb;
    vector<double> vErrEffPerPlaneBoth_PbPb, vErrEffPerPlaneBP_PbPb, vErrEffPerPlaneNBP_PbPb;

    vector<vector<double>> vEffPerPlanePerRunBoth_PbPb, vEffPerPlanePerRunBP_PbPb, vEffPerPlanePerRunNBP_PbPb;
    vector<vector<double>> vErrEffPerPlanePerRunBoth_PbPb, vErrEffPerPlanePerRunBP_PbPb, vErrEffPerPlanePerRunNBP_PbPb;
    
    plane = 0; //Variable to be used in the following loop to keep track of the plane (MT11, MT12, MT21, MT22)

    for (int i = 0; i <= (vEffBoth_Planes_runs_PbPb.size())*4; i++) {
        
        if ((i % vEffBoth_Planes_runs_PbPb.size() == 0) && (i != 0)) {
            //cout << "Pushing back vector" << endl;
            //cout << i-(plane*vEffBoth_Planes_runs.size()) << "\t" << plane << endl;
            if (plane > 0) {
                vEffPerPlaneBoth_PbPb.push_back(vEffBoth_Planes_runs_PbPb[i-(plane*vEffBoth_Planes_runs_PbPb.size())-1][plane]);
                vEffPerPlaneBP_PbPb.push_back(vEffBP_Planes_runs_PbPb[i-(plane*vEffBP_Planes_runs_PbPb.size())-1][plane]);
                vEffPerPlaneNBP_PbPb.push_back(vEffNBP_Planes_runs_PbPb[i-(plane*vEffNBP_Planes_runs_PbPb.size())-1][plane]);
                vErrEffPerPlaneBoth_PbPb.push_back(vErrEffBoth_Planes_runs_PbPb[i-(plane*vErrEffBoth_Planes_runs_PbPb.size())-1][plane]);
                vErrEffPerPlaneBP_PbPb.push_back(vErrEffBP_Planes_runs_PbPb[i-(plane*vErrEffBP_Planes_runs_PbPb.size())-1][plane]);
                vErrEffPerPlaneNBP_PbPb.push_back(vErrEffNBP_Planes_runs_PbPb[i-(plane*vErrEffNBP_Planes_runs_PbPb.size())-1][plane]);
                //cout << "pushing back element " << vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())-1][plane] << endl;
            }
            vEffPerPlanePerRunBoth_PbPb.push_back(vEffPerPlaneBoth_PbPb);
            vEffPerPlanePerRunBP_PbPb.push_back(vEffPerPlaneBP_PbPb);
            vEffPerPlanePerRunNBP_PbPb.push_back(vEffPerPlaneNBP_PbPb);
            vErrEffPerPlanePerRunBoth_PbPb.push_back(vErrEffPerPlaneBoth_PbPb);
            vErrEffPerPlanePerRunBP_PbPb.push_back(vErrEffPerPlaneBP_PbPb);
            vErrEffPerPlanePerRunNBP_PbPb.push_back(vErrEffPerPlaneNBP_PbPb);
            
            vEffPerPlaneBoth_PbPb.clear();
            vEffPerPlaneBP_PbPb.clear();
            vEffPerPlaneNBP_PbPb.clear();
            vErrEffPerPlaneBoth_PbPb.clear();
            vErrEffPerPlaneBP_PbPb.clear();
            vErrEffPerPlaneNBP_PbPb.clear();
            
            plane++; //move to the next plane
        }
        
        else {
            //cout << "pushing back element " << vEffBoth_Planes_runs[i-(plane*vEffBoth_Planes_runs.size())][plane] << "\t i \t" << i << "\t plane \t" << plane << endl;
            //cout << "pushing back element with index " << i-(plane*vEffBoth_Planes_runs.size()) << "\t i \t" << i << "\t plane \t" << plane << endl;
            if (plane == 0) {
                //Eff per plane on both planes, BP and NBP
                vEffPerPlaneBoth_PbPb.push_back(vEffBoth_Planes_runs_PbPb[i-(plane*vEffBoth_Planes_runs_PbPb.size())][plane]);
                vEffPerPlaneBP_PbPb.push_back(vEffBP_Planes_runs_PbPb[i-(plane*vEffBP_Planes_runs_PbPb.size())][plane]);
                vEffPerPlaneNBP_PbPb.push_back(vEffNBP_Planes_runs_PbPb[i-(plane*vEffNBP_Planes_runs_PbPb.size())][plane]);
                //Err on eff per plane on both planes, BP and NBP
                vErrEffPerPlaneBoth_PbPb.push_back(vErrEffBoth_Planes_runs_PbPb[i-(plane*vErrEffBoth_Planes_runs_PbPb.size())][plane]);
                vErrEffPerPlaneBP_PbPb.push_back(vErrEffBP_Planes_runs_PbPb[i-(plane*vErrEffBP_Planes_runs_PbPb.size())][plane]);
                vErrEffPerPlaneNBP_PbPb.push_back(vErrEffNBP_Planes_runs_PbPb[i-(plane*vErrEffNBP_Planes_runs_PbPb.size())][plane]);
            }
            else {
                //Eff per plane on both planes, BP and NBP
                vEffPerPlaneBoth_PbPb.push_back(vEffBoth_Planes_runs_PbPb[i-(plane*vEffBoth_Planes_runs_PbPb.size())-1][plane]);
                vEffPerPlaneBP_PbPb.push_back(vEffBP_Planes_runs_PbPb[i-(plane*vEffBP_Planes_runs_PbPb.size())-1][plane]);
                vEffPerPlaneNBP_PbPb.push_back(vEffNBP_Planes_runs_PbPb[i-(plane*vEffNBP_Planes_runs_PbPb.size())-1][plane]);
                //Err on eff per plane on both planes, BP and NBP
                vErrEffPerPlaneBoth_PbPb.push_back(vErrEffBoth_Planes_runs_PbPb[i-(plane*vErrEffBoth_Planes_runs_PbPb.size())-1][plane]);
                vErrEffPerPlaneBP_PbPb.push_back(vErrEffBP_Planes_runs_PbPb[i-(plane*vErrEffBP_Planes_runs_PbPb.size())-1][plane]);
                vErrEffPerPlaneNBP_PbPb.push_back(vErrEffNBP_Planes_runs_PbPb[i-(plane*vErrEffNBP_Planes_runs_PbPb.size())-1][plane]);
            }   
        }
    }

    cout << vEffPerPlanePerRunBoth_PbPb[0].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[1].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[2].size() << "\t" << vEffPerPlanePerRunBoth_PbPb[3].size() << endl;
    cout << vEffPerPlanePerRunBP_PbPb[0].size() << "\t" << vEffPerPlanePerRunBP_PbPb[1].size() << "\t" << vEffPerPlanePerRunBP_PbPb[2].size() << "\t" << vEffPerPlanePerRunBP_PbPb[3].size() << endl;
    cout << vEffPerPlanePerRunNBP_PbPb[0].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[1].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[2].size() << "\t" << vEffPerPlanePerRunNBP_PbPb[3].size() << endl;

    //Plane eff vs IR
    vector<TGraphErrors*> gEffPlaneBothPlanesIR, gEffPlaneBPIR, gEffPlaneNBPIR;
    
    for (int i = 0; i < nBinsPlane; i++) {
        TGraphErrors *g1 = new TGraphErrors(vIR_PbPb2023.size(),&vIR_PbPb2023[0],&vEffPerPlanePerRunBoth_PbPb[i][0],NULL,&vErrEffPerPlanePerRunBoth_PbPb[i][0]);
        TGraphErrors *g2 = new TGraphErrors(vIR_PbPb2023.size(),&vIR_PbPb2023[0],&vEffPerPlanePerRunBP_PbPb[i][0],NULL,&vErrEffPerPlanePerRunBP_PbPb[i][0]);
        TGraphErrors *g3 = new TGraphErrors(vIR_PbPb2023.size(),&vIR_PbPb2023[0],&vEffPerPlanePerRunNBP_PbPb[i][0],NULL,&vErrEffPerPlanePerRunNBP_PbPb[i][0]);
        
        gEffPlaneBothPlanesIR.push_back(g1);
        gEffPlaneBPIR.push_back(g2);
        gEffPlaneNBPIR.push_back(g3);
    } 

    //From run 544868 onwards the magnet polarity has been changed from +/+ to -/-
    //Let's try to color the markers differently in the Eff v IR plot

    //Both planes
    TCanvas *cEffPlaneBothPlanesIR = new TCanvas();
    cEffPlaneBothPlanesIR->Divide(1,4);
    for (int i = 0; i < nBinsPlane; i++) {
        cEffPlaneBothPlanesIR->cd(i+1);
        gEffPlaneBothPlanesIR.at(i)->SetMarkerStyle(8);
        gEffPlaneBothPlanesIR.at(i)->SetMarkerSize(1);
        gEffPlaneBothPlanesIR.at(i)->SetMarkerColor(kBlack);
        gEffPlaneBothPlanesIR.at(i)->SetTitle((planeName[i] + " Both planes vs IR").c_str());
        gEffPlaneBothPlanesIR.at(i)->GetXaxis()->SetTitle("IR [kHz]]");
        gEffPlaneBothPlanesIR.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPlaneBothPlanesIR.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPlaneBothPlanesIR.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPlaneBothPlanesIR.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPlaneBothPlanesIR.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPlaneBothPlanesIR.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPlaneBothPlanesIR.at(i)->Draw("AP");
    }
    //BP
    TCanvas *cEffPlaneBPIR = new TCanvas();
    cEffPlaneBPIR->Divide(1,4);
    for (int i = 0; i < nBinsPlane; i++) {
        cEffPlaneBPIR->cd(i+1);
        gEffPlaneBPIR.at(i)->SetMarkerStyle(8);
        gEffPlaneBPIR.at(i)->SetMarkerSize(1);
        gEffPlaneBPIR.at(i)->SetMarkerColor(kRed);
        gEffPlaneBPIR.at(i)->SetTitle((planeName[i] + " BP vs IR").c_str());
        gEffPlaneBPIR.at(i)->GetXaxis()->SetTitle("IR [kHz]]");
        gEffPlaneBPIR.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPlaneBPIR.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPlaneBPIR.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPlaneBPIR.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPlaneBPIR.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPlaneBPIR.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPlaneBPIR.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPlaneBPIR.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPlaneBPIR.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPlaneBPIR.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPlaneBPIR.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPlaneBPIR.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPlaneBPIR.at(i)->Draw("AP");
    }
    //NBP
    TCanvas *cEffPlaneNBPIR = new TCanvas();
    cEffPlaneNBPIR->Divide(1,4);
    for (int i = 0; i < nBinsPlane; i++) {
        cEffPlaneNBPIR->cd(i+1);
        gEffPlaneNBPIR.at(i)->SetMarkerStyle(8);
        gEffPlaneNBPIR.at(i)->SetMarkerSize(1);
        gEffPlaneNBPIR.at(i)->SetMarkerColor(kBlack);
        gEffPlaneNBPIR.at(i)->SetTitle((planeName[i] + " NBP vs IR").c_str());
        gEffPlaneNBPIR.at(i)->GetXaxis()->SetTitle("IR [kHz]]");
        gEffPlaneNBPIR.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPlaneNBPIR.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPlaneNBPIR.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPlaneNBPIR.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPlaneNBPIR.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPlaneNBPIR.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPlaneNBPIR.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPlaneNBPIR.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPlaneNBPIR.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPlaneNBPIR.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPlaneNBPIR.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPlaneNBPIR.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPlaneNBPIR.at(i)->Draw("AP");
    }

    
    //Vector of TGraphErrors for all 4 planes on both, BP and NBP - plane eff vs run #
    vector<TGraphErrors*> gEffPerRunPlaneBoth, gEffPerRunPlaneBP, gEffPerRunPlaneNBP;
    
    for (int i = 0; i < nBinsPlane; i++) {
        TGraphErrors *g1 = new TGraphErrors(vRun_2023.size(),&vRun_2023[0],&vEffPerPlanePerRunBoth[i][0],NULL,&vErrEffPerPlanePerRunBoth[i][0]);
        TGraphErrors *g2 = new TGraphErrors(vRun_2023.size(),&vRun_2023[0],&vEffPerPlanePerRunBP[i][0],NULL,&vErrEffPerPlanePerRunBP[i][0]);
        TGraphErrors *g3 = new TGraphErrors(vRun_2023.size(),&vRun_2023[0],&vEffPerPlanePerRunNBP[i][0],NULL,&vErrEffPerPlanePerRunNBP[i][0]);
        
        gEffPerRunPlaneBoth.push_back(g1);
        gEffPerRunPlaneBP.push_back(g2);
        gEffPerRunPlaneNBP.push_back(g3);
    } 

    bool isTime = false;

    //Both
    TCanvas *cEffPerPlaneBoth = new TCanvas();
    cEffPerPlaneBoth->Divide(1,4);
    for (int i = 0; i < nBinsPlane; i++) {
        cEffPerPlaneBoth->cd(i+1);
        gEffPerRunPlaneBoth.at(i)->SetMarkerStyle(8);
        gEffPerRunPlaneBoth.at(i)->SetMarkerSize(1);
        gEffPerRunPlaneBoth.at(i)->SetMarkerColor(kBlack);
        gEffPerRunPlaneBoth.at(i)->SetTitle((planeName[i] + " Both planes").c_str());
        if (isTime) {
            gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTimeDisplay(1);
            gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetNdivisions(503);
            gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTimeFormat("%Y-%m-%d");
            gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTimeOffset(0,"gmt");
            gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTitle("Time [UTC]");
        }
        else {
            gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTitle("Run #");
        }
        gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPerRunPlaneBoth.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunPlaneBoth.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunPlaneBoth.at(i)->Draw("AP");
    }
    
    //BP
    TCanvas *cEffPerPlaneBP = new TCanvas();
    cEffPerPlaneBP->Divide(1,4);
    for (int i = 0; i < nBinsPlane; i++) {
        cEffPerPlaneBP->cd(i+1);
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
        }
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPerRunPlaneBP.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPerRunPlaneBP.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunPlaneBP.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunPlaneBP.at(i)->Draw("AP");
    }
    
    //NBP
    TCanvas *cEffPerPlaneNBP = new TCanvas();
    cEffPerPlaneNBP->Divide(1,4);
    for (int i = 0; i < nBinsPlane; i++) {
        cEffPerPlaneNBP->cd(i+1);
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
        }
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitleOffset(0.5);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitleSize(0.07);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetTitleFont(62);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetLabelSize(0.07);
        gEffPerRunPlaneNBP.at(i)->GetXaxis()->SetLabelFont(62);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitle("Efficiency [%]");
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitleOffset(0.35);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitleSize(0.07);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetTitleFont(62);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetLabelSize(0.07);
        gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetLabelFont(62);
        //gEffPerRunPlaneNBP.at(i)->GetYaxis()->SetRangeUser(85,100);
        gEffPerRunPlaneNBP.at(i)->Draw("AP");
    }

    //per RPC
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
    
    hRun_pp2023.close();
    hDate_pp2023.close();
    hRun_PbPb2023.close();
    hDate_PbPb2023.close();
}