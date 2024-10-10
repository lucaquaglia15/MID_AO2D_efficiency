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

    double effBothPlane = 0, effBPPlane = 0, effNBPPlane =0;
    double errEffBothPlane = 0, errEffBPPlane = 0, errEffNBPPlane = 0;

    //Run by run
    vector<double> vEffBoth_Planes, vEffBP_Planes, vEffNBP_Planes;
    vector<double> vErrEffBoth_Planes, vErrEffBP_Planes, vErrEffNBP_Planes;

    //Run by run
    vector<vector<double>> vEffBoth_Planes_runs, vEffBP_Planes_runs, vEffNBP_Planes_runs;
    vector<vector<double>> vErrEffBoth_Planes_runs, vErrEffBP_Planes_runs, vErrEffNBP_Planes_runs;
    
    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    
    string globalPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/";

    /*ifstream hLowEff;
    hLowEff.open(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+".txt").c_str());

    string lowEffType, planeType, plusMinus;
    int LBnumber;
    double lowEfficiency, errLowEfficiency;
    vector<int> LBfiftyBP, LBfiftyNBP; //LB with efficiency < 50% on BP and NBP
    vector<int> LBfifty_eightyBP, LBfifty_eightyNBP; //LB with 50% < efficiency < 80% on BP and NBP
    vector<int> LBeighty_ninetyBP, LBeighty_ninetyNBP; //LB with 80% < efficiency < 90% on BP and NBP
    vector<int> LBninety_ninetyfiveBP, LBninety_ninetyfiveNBP; //LB with 90% < efficiency < 95% on BP and NBP
    vector<int> LBninetyfiveBP, LBninetyfiveNBP; //LB with efficiency > 95% on BP and NBP

    while (hLowEff >> lowEffType >> LBnumber >> planeType >> lowEfficiency >> errLowEfficiency) {
        //cout << lowEffType << "\t" << LBnumber << "\t" << planeType << "\t" << lowEfficiency << "\t" << errLowEfficiency << endl;
        if (lowEffType == "50" && planeType == "BP") {
            LBfiftyBP.push_back(LBnumber);
        }
        else if (lowEffType == "50" && planeType == "NBP") {
            LBfiftyNBP.push_back(LBnumber);
        }
        else if(lowEffType == "50_80" && planeType == "BP") {
            LBfifty_eightyBP.push_back(LBnumber);
        }
        else if(lowEffType == "50_80" && planeType == "NBP") {
            LBfifty_eightyNBP.push_back(LBnumber);
        }
        else if(lowEffType == "80_90" && planeType == "BP") {
            LBeighty_ninetyBP.push_back(LBnumber);
        }
        else if(lowEffType == "80_90" && planeType == "NBP") {
            LBeighty_ninetyNBP.push_back(LBnumber);
        }
        else if(lowEffType == "90_95" && planeType == "BP") {
            LBninety_ninetyfiveBP.push_back(LBnumber);
        }
        else if(lowEffType == "90_95" && planeType == "NBP") {
            LBninety_ninetyfiveNBP.push_back(LBnumber);
        }
        else if(lowEffType == "95" && planeType == "BP") {
            LBninetyfiveBP.push_back(LBnumber);
        }
        else if(lowEffType == "95" && planeType == "NBP") {
            LBninetyfiveNBP.push_back(LBnumber);
        }
    }
    //cout << LBfiftyBP.size() << "\t" << LBfiftyNBP.size() << endl;
    //Placeholder labels for histogram
    vector<int> binsX;
    for (int i = 0; i < 5; i++) {
        binsX.push_back(i);
    }

    vector<string> labelsEff {"50","50_80","80_90","90_95","95"};

    //Eff ranges per LB on BP
    TH1F *BPboardsEff = new TH1F("BPboardsEff","BPboardsEff",binsX.size(),binsX.front()+0.5,binsX.back()+0.5);
    BPboardsEff->SetBinContent(1,LBfiftyBP.size());
    BPboardsEff->SetBinContent(2,LBfifty_eightyBP.size());
    BPboardsEff->SetBinContent(3,LBeighty_ninetyBP.size());
    BPboardsEff->SetBinContent(4,LBninety_ninetyfiveBP.size());
    BPboardsEff->SetBinContent(5,LBninetyfiveBP.size());

    BPboardsEff->GetXaxis()->SetTitle("BP efficiency [%]");
    BPboardsEff->GetYaxis()->SetTitle("Normalized counts");
    BPboardsEff->Scale(1/BPboardsEff->Integral());
    BPboardsEff->GetYaxis()->SetRangeUser(0,1);

    for (int i = 0; i < binsX.size(); i++) {
        BPboardsEff->GetXaxis()->SetBinLabel(i+1,(labelsEff.at(i)).c_str());
    }

    TCanvas *cEffRangeBP = new TCanvas();
    cEffRangeBP->cd();
    cEffRangeBP->SetTitle(("BPboardsEff_"+ period).c_str());
    cEffRangeBP->SetCanvasSize(1200,1200);
    BPboardsEff->SetStats(0);
    BPboardsEff->SetLineWidth(3);    
    BPboardsEff->Draw("HISTO");
    cEffRangeBP->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/BPboardsEff_"+ period+".png").c_str());

    //Eff ranges per LB on NBP
    TH1F *NBPboardsEff = new TH1F("NBPboardsEff","NBPboardsEff",binsX.size(),binsX.front()+0.5,binsX.back()+0.5);
    NBPboardsEff->SetBinContent(1,LBfiftyNBP.size());
    NBPboardsEff->SetBinContent(2,LBfifty_eightyNBP.size());
    NBPboardsEff->SetBinContent(3,LBeighty_ninetyNBP.size());
    NBPboardsEff->SetBinContent(4,LBninety_ninetyfiveNBP.size());
    NBPboardsEff->SetBinContent(5,LBninetyfiveNBP.size());

    NBPboardsEff->GetXaxis()->SetTitle("NBP efficiency [%]");
    NBPboardsEff->GetYaxis()->SetTitle("Normalized counts");
    NBPboardsEff->Scale(1/NBPboardsEff->Integral());
    NBPboardsEff->GetYaxis()->SetRangeUser(0,1);

    for (int i = 0; i < binsX.size(); i++) {
        NBPboardsEff->GetXaxis()->SetBinLabel(i+1,(labelsEff.at(i)).c_str());
    }

    TCanvas *cEffRangeNBP = new TCanvas();
    cEffRangeNBP->cd();
    cEffRangeNBP->SetTitle(("NBPboardsEff_"+ period).c_str());
    cEffRangeNBP->SetCanvasSize(1200,1200);
    NBPboardsEff->SetStats(0);
    NBPboardsEff->SetLineWidth(3);
    NBPboardsEff->Draw("HISTO");
    cEffRangeNBP->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/NBPboardsEff_"+ period+".png").c_str());*/

    //Path of the merged file, run-by-run
    string runPath = globalPath+"runs/";

    //Path for the .txt file of the run list of the period
    string runNumbers = globalPath+"run_list.txt"; 
    string runDates = globalPath+"run_dates.txt";

    //Open txt file of runs
    ifstream hRun;
    hRun.open(runNumbers.c_str());

    //Open txt file of start/end dates of the runs
    ifstream hDate;
    hDate.open(runDates.c_str());
    //hDate.open("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/LHC23_PbPb_pass3_I-A11/run_dates.txt");
    
    //Get start and end of each run
    long int runForDate;
    double start, end;
    vector<long int> vRunForDate;
    vector<double> vStart, vEnd;
    
    while (hDate >> runForDate >> start >> end){
        vRunForDate.push_back(runForDate);
        vStart.push_back(start);
        vEnd.push_back(end);
    }

    //Push back to a vector of int (no need to care about size)
    double run;
    vector<double> vRun;

    while(hRun >> run) {
        vRun.push_back(run);
    }
    //sort in ascending order
    sort(vRun.begin(), vRun.end()); 

    //General string name
    string fileName = "AnalysisResults.root";

    //Loop on all runs
    for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
        //Enter the folder
        string runFolder = runPath+to_string((int)vRun.at(iRun));

        //run file name = path of the folder + run number (runFolder) + fileName
        string runFileName = runFolder+"/"+fileName;
        TFile *fRun = new TFile(runFileName.c_str(),"READ");

        TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

        TH1F *hFiredBoth_Planes = (TH1F*)d->Get("nFiredBothperPlane");
        TH1F *hFiredBP_Planes = (TH1F*)d->Get("nFiredBPperPlane");
        TH1F *hFiredNBP_Planes = (TH1F*)d->Get("nFiredNBPperPlane");
        TH1F *hTotPlanes = (TH1F*)d->Get("nTotperPlane");    

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
            
            //Fill vector for error on efficiency per LB in the run
            vErrEffBoth_Planes.push_back(errEffBothPlane);
            vErrEffBP_Planes.push_back(errEffBPPlane);
            vErrEffNBP_Planes.push_back(errEffNBPPlane);
        
        }

        //Push back the vector with the eff of LB to a larger vector of vectors (one element of this = one run)
        vEffBoth_Planes_runs.push_back(vEffBoth_Planes);
        vEffBP_Planes_runs.push_back(vEffBP_Planes);
        vEffNBP_Planes_runs.push_back(vEffNBP_Planes);
        
        //Push back the vector with the error on eff of LB to a larger vector of vectors (one element of this = one run)
        vErrEffBoth_Planes_runs.push_back(vErrEffBoth_Planes);
        vErrEffBP_Planes_runs.push_back(vErrEffBP_Planes); 
        vErrEffNBP_Planes_runs.push_back(vErrEffNBP_Planes);

        //Clear vector of eff for LB in a run
        vEffBoth_Planes.clear();
        vEffBP_Planes.clear();
        vEffNBP_Planes.clear();
        //Clear vector of error on eff for LB in a run
        vErrEffBoth_Planes.clear();
        vErrEffBP_Planes.clear();
        vErrEffNBP_Planes.clear();

    } //End of loop on all runs
    cout << "Size of vStart " << vStart.size() << endl;
    cout << "Size of vEffBoth_Planes " << vEffBoth_Planes_runs.size() << endl;

    for (unsigned int i = 0; i < vEffBoth_Planes_runs.size(); i++) {
        cout << vEffBoth_Planes_runs[i][0] << "+-" << vErrEffBoth_Planes_runs[i][0] << endl;
    }

    vector<double> vEffPerPlaneBoth, vEffPerPlaneBP, vEffPerPlaneNBP;
    vector<double> vErrEffPerPlaneBoth, vErrEffPerPlaneBP, vErrEffPerPlaneNBP;

    vector<vector<double>> vEffPerPlanePerRunBoth, vEffPerPlanePerRunBP, vEffPerPlanePerRunNBP;
    vector<vector<double>> vErrEffPerPlanePerRunBoth, vErrEffPerPlanePerRunBP, vErrEffPerPlanePerRunNBP;
    
    int j = 0; //Variable to be used in the following loop to keep track of the plane (MT11, MT12, MT21, MT22)
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
            //cout << i-(j*vEffBoth_Planes_runs.size()) << "\t" << j << endl;
            if (j > 0) {
                vEffPerPlaneBoth.push_back(vEffBoth_Planes_runs[i-(j*vEffBoth_Planes_runs.size())-1][j]);
                //cout << "pushing back element " << vEffBoth_Planes_runs[i-(j*vEffBoth_Planes_runs.size())-1][j] << endl;
            }
            vEffPerPlanePerRunBoth.push_back(vEffPerPlaneBoth);
            vEffPerPlaneBoth.clear();
            j++;
        }
        
        else {
            //cout << "pushing back element " << vEffBoth_Planes_runs[i-(j*vEffBoth_Planes_runs.size())][j] << "\t i \t" << i << "\t j \t" << j << endl;
            //cout << "pushing back element with index " << i-(j*vEffBoth_Planes_runs.size()) << "\t i \t" << i << "\t j \t" << j << endl;
            if (j == 0) {
                vEffPerPlaneBoth.push_back(vEffBoth_Planes_runs[i-(j*vEffBoth_Planes_runs.size())][j]);
            }
            else {
                vEffPerPlaneBoth.push_back(vEffBoth_Planes_runs[i-(j*vEffBoth_Planes_runs.size())-1][j]);
            }   
        }
    } 

    cout << vEffPerPlanePerRunBoth[0].size() << "\t" << vEffPerPlanePerRunBoth[1].size() << "\t" << vEffPerPlanePerRunBoth[2].size() << "\t" << vEffPerPlanePerRunBoth[3].size() << endl;

    TGraphErrors *gEffPerRunPlane11Both = new TGraphErrors(vStart.size(),&vStart[0],&vEffPerPlanePerRunBoth[0][0],NULL,NULL);
    TGraphErrors *gEffPerRunPlane12Both = new TGraphErrors(vStart.size(),&vStart[0],&vEffPerPlanePerRunBoth[1][0],NULL,NULL);
    TGraphErrors *gEffPerRunPlane21Both = new TGraphErrors(vStart.size(),&vStart[0],&vEffPerPlanePerRunBoth[2][0],NULL,NULL);
    TGraphErrors *gEffPerRunPlane22Both = new TGraphErrors(vStart.size(),&vStart[0],&vEffPerPlanePerRunBoth[3][0],NULL,NULL);

    /*TGraphErrors *gEffPerRunPlane11BP = new TGraphErrors();
    TGraphErrors *gEffPerRunPlane12BP = new TGraphErrors();
    TGraphErrors *gEffPerRunPlane21BP = new TGraphErrors();
    TGraphErrors *gEffPerRunPlane22BP = new TGraphErrors();

    TGraphErrors *gEffPerRunPlane11NBP = new TGraphErrors();
    TGraphErrors *gEffPerRunPlane12NBP = new TGraphErrors();
    TGraphErrors *gEffPerRunPlane21NBP = new TGraphErrors();
    TGraphErrors *gEffPerRunPlane22NBP = new TGraphErrors();*/

    TCanvas *cEffPerPlaneBoth = new TCanvas();
    //cEffPerPlaneBoth->SetCanvasSize(600,600);
    cEffPerPlaneBoth->Divide(1,4);
    cEffPerPlaneBoth->cd(1); //MT11
    gEffPerRunPlane11Both->SetMarkerStyle(8);
    gEffPerRunPlane11Both->SetMarkerSize(2);
    gEffPerRunPlane11Both->Draw("AP");
    
    cEffPerPlaneBoth->cd(2); //MT12
    gEffPerRunPlane12Both->SetMarkerStyle(8);
    gEffPerRunPlane12Both->SetMarkerSize(2);
    gEffPerRunPlane12Both->Draw("AP");

    cEffPerPlaneBoth->cd(3); //MT21
    gEffPerRunPlane21Both->SetMarkerStyle(8);
    gEffPerRunPlane21Both->SetMarkerSize(2);
    gEffPerRunPlane21Both->Draw("AP");

    cEffPerPlaneBoth->cd(4); //MT22
    gEffPerRunPlane22Both->SetMarkerStyle(8);
    gEffPerRunPlane22Both->SetMarkerSize(2);
    gEffPerRunPlane22Both->Draw("AP");

    
    hRun.close();
    hDate.close();
}