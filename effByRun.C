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

//#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

bool debug = false;
const int trackGoal = 8e+6;
int tracks = 0;
int cumulativeTracks = 0;
int nBinsPlane = 4; //Number of planes
int nBinsRPC = 72; //Number of RPCs
int nBinsBoard = 936; //Number of LBs

//Function to delete folders with data from multiple runs merge -> called everytime the code is launched to clean up in case we change
//target number of tracks, hence number of merges
void clearFolders(string folder) {

    string extension = "mergedRuns";

    for (auto const& dir_entry :std::filesystem::directory_iterator{folder}) { //Loop through all folders in the merged_files folder
        if ((string(dir_entry.path())).find(extension) != string::npos) { //If the given folder contains the extension in its name
            filesystem::remove_all((string(dir_entry.path())).c_str()); //Delete folder to be ready for new iteration 
        }
    }
} //End of clearFolders

//Function to calculate efficiency -> defined as a function not to write over and over the for loop
void calculateEfficiencyByRun(string runFolder) {

    bool enablePrint = false; //To enable eff printout
    
    TFile *fIn = new TFile(runFolder.c_str(),"READ"); //Open the merged files from several runs

    TDirectoryFile *d = (TDirectoryFile*)fIn->Get("mid-efficiency");

    TH1F *hFiredBothPlanesLB = (TH1F*)d->Get("nFiredBothperBoard");
    TH1F *hFiredBPLB = (TH1F*)d->Get("nFiredBPperBoard");
    TH1F *hFiredNBPLB = (TH1F*)d->Get("nFiredNBPperBoard");
    TH1F *hTotLB = (TH1F*)d->Get("nTotperBoard");

    //To calculate efficiency and error on efficiency
    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    //Histograms to visualize easily the efficiency distribution
    TH1F *hEffLB_both = new TH1F("hEffLB_both","hEffLB_both",220,0.,110.);
    TH1F *hEffLB_BP = new TH1F("hEffLB_BP","hEffLB_BP",220,0.,110.);
    TH1F *hEffLB_NBP = new TH1F("hEffLB_NBP","hEffLB_NBP",220,0.,110.);
    //Histograms to visualize easily the error on efficiency distribution
    TH1F *hErrEffLB_both = new TH1F("hErrEffLB_both","hErrEffLB_both",200,0.,10.);
    TH1F *hErrEffLB_BP = new TH1F("hErrEffLB_BP","hErrEffLB_BP",200,0.,10.);
    TH1F *hErrEffLB_NBP = new TH1F("hErrEffLB_NBP","hErrEffLB_NBP",200,0.,10.);

    //Loop through Local Boards
    for (int i = 1; i <= nBinsBoard; i++) {
        //cout << "LB Both planes " <<  i << "\t" << hFiredBothPlanesLB->GetBinContent(i) << "\t" << hTotLB->GetBinContent(i) << endl;

        if (hTotLB->GetBinContent(i) != 0) {

            effBothLB = (hFiredBothPlanesLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
            effBPLB = (hFiredBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
            effNBPLB = (hFiredNBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;

            errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLB->GetBinContent(i));
            errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLB->GetBinContent(i));
            errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLB->GetBinContent(i));

            if (enablePrint) {
                cout << "LB " << i+1 << endl;
                cout << "Eff both " << effBothLB << " +- " << errEffBothLB << endl; 
                cout << "Eff BP " << effBPLB << " +- " << errEffBPLB << endl; 
                cout << "Eff NBP " << effNBPLB << " +- " << errEffNBPLB << endl; 
            }

            hEffLB_both->Fill(effBothLB); 
            hEffLB_BP->Fill(effBPLB); 
            hEffLB_NBP->Fill(effNBPLB); 

            hErrEffLB_both->Fill(errEffBothLB); 
            hErrEffLB_BP->Fill(errEffBPLB); 
            hErrEffLB_NBP->Fill(errEffNBPLB); 
        }
    }

    //Eff on both planes
    TCanvas *cEffLB_both = new TCanvas();
    cEffLB_both->cd();
    hEffLB_both->Draw("HISTO");
    //Eff on BP
    TCanvas *cEffLB_BP = new TCanvas();
    cEffLB_BP->cd();
    hEffLB_BP->Draw("HISTO");
    //Eff on NBP
    TCanvas *cEffLB_NBP = new TCanvas();
    cEffLB_NBP->cd();
    hEffLB_NBP->Draw("HISTO");
    //Err on Eff on both planes
    TCanvas *cErrEffLB_both = new TCanvas();
    cErrEffLB_both->cd();
    hErrEffLB_both->Draw("HISTO");
    //Error on Eff on BP
    TCanvas *cErrEffLB_BP = new TCanvas();
    cErrEffLB_BP->cd();
    hErrEffLB_BP->Draw("HISTO");
    //Error on Eff on NBP
    TCanvas *cErrEffLB_NBP = new TCanvas();
    cErrEffLB_NBP->cd();
    hErrEffLB_NBP->Draw("HISTO");
    
    //We need the path to create an output .root file on which to save the eff and error on eff hisotgrams
    //We start from thr merged file path and we remove the "AnalysisResults.root" part at the end (20 characters)
    //Using the erase function from the string library
    runFolder.erase(runFolder.end()-20,runFolder.end());
    cout << "Test of removing string portion: " << runFolder << endl;

    //Create output file in the merged files directory and save the hisotgrams on it
    TFile *fOut = new TFile((runFolder+"mergedRuns.root").c_str(),"RECREATE");
    cEffLB_both->Write("Eff_LB_both");
    cEffLB_BP->Write("Eff_LB_BP");
    cEffLB_NBP->Write("Eff_LB_NBP");
    cErrEffLB_both->Write("Err_Eff_LB_both");
    cErrEffLB_BP->Write("Err_Eff_LB_BP");
    cErrEffLB_NBP->Write("Err_Eff_LB_NBP");

} //End of calculateEfficieny function

void effByRun() { //Main function

    bool open = false; //used to keep track if the .dat file (updated to contain the list of .root files to be merged) has been opened once

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    vector<float> vEffBothLB, vEffBPLB, vEffNBPLB;
    vector<float> vErrEffBothLB, vErrEffBPLB, vErrEffNBPLB;

    vector<vector<float>> vEffBothLB_runs, vEffBPLB_runs, vEffNBPLB_runs;
    vector<vector<float>> vErrEffBothLB_runs, vErrEffBPLB_runs, vErrEffNBPLB_runs;

    //Plane name
    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //General path to add flexibility to the code + period name
    string period = "LHC23_pass4_skimmed_QC1";
    string globalPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/";

    //Path of the merged file, run-by-run
    string runPath = globalPath+"merged_files/"; 

    //Path for the .txt file of the run list of the period
    string runNumbers = globalPath+"run_list.txt"; 

    //Call function to clear all the folders with the data merged from multiple runs
    clearFolders(runPath);

    //Open txt file of runs
    ifstream hRun;
    hRun.open(runNumbers.c_str());

    //Push back to a vector of int (no need to care about size)
    float run;
    vector<float> vRun;

    while(hRun >> run) {
        vRun.push_back(run);
    }

    //Output file for the merge of root files if the number of tracks reaches the desired goal
    ofstream hMergeRuns;

    //General string name
    string fileName = "AnalysisResults.root";

    //Load hadd.C macro to merge the root files from different runs
    gROOT->ProcessLine(".L /home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/hadd.C");

    //Test - to keep track of the number of merges
    int mergeCounter = 0;

    //Loop on all runs
    for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {
        //Enter the folder
        string runFolder = runPath+to_string((int)vRun.at(iRun));

        //run file name = path of the folder + run number (runFolder) + fileName
        string runFileName = runFolder+"/"+fileName;
        TFile *fRun = new TFile(runFileName.c_str(),"READ");

        TDirectoryFile *d = (TDirectoryFile*)fRun->Get("mid-efficiency");

        TH1F *hFiredBothPlanesLB = (TH1F*)d->Get("nFiredBothperBoard");
        TH1F *hFiredBPLB = (TH1F*)d->Get("nFiredBPperBoard");
        TH1F *hFiredNBPLB = (TH1F*)d->Get("nFiredNBPperBoard");
        TH1F *hTotLB = (TH1F*)d->Get("nTotperBoard");

        tracks = hTotLB->GetEntries();
        cumulativeTracks += tracks;
        cout << "Run number " << vRun.at(iRun) << " tot tracks in all LB " << tracks << " cumulative " << cumulativeTracks << endl;

        if (cumulativeTracks < trackGoal) { //If total track number is below the target -> Fill the file with the path of each AnalysisResults.root from each run
            //Open the output file only if it has not been opened (i.e. open == false)
            //meaning that it's the first run to be analyzed
            if (!open) {
                hMergeRuns.open("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/merged_files/runs.dat");
                open = true;
            }
            
            hMergeRuns << runFileName << "\n";
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
            hMergeRuns << runFileName << "\n";
            mergeCounter++;
            cumulativeTracks = 0;
            hMergeRuns.close();
            open = false;
            //Test - create a sub-folder inside the merged_files directory
            gSystem->mkdir((runPath+"mergedRuns"+to_string(mergeCounter)).c_str());
            //Execute the hadd code (already loaded before the loop)
            string mergeFilesForHadd = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root"+'"';
            string mergeFiles = runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root";
            cout << "mergeFilesForHadd: " << mergeFilesForHadd << endl;
            gROOT->ProcessLine(Form("hadd(%s)",mergeFilesForHadd.c_str()));
            calculateEfficiencyByRun(mergeFiles);
        }

        //Edge case, it's the last run of the list and the number of cumulative tracks has not yet reached the goal
            // we still merge those files and that's it so we copied the same code as the merges earlier on
        if (vRun.at(iRun) == vRun.back() && cumulativeTracks < trackGoal) {
            cout << "Last run of the list and track goal not yet reached, merging anyway" << endl;
            hMergeRuns << runFileName << "\n";
            mergeCounter++;
            cumulativeTracks = 0;
            hMergeRuns.close();
            open = false; //no need to set it to false but for redudancy we do it
            //Test - create a sub-folder inside the merged_files directory
            gSystem->mkdir((runPath+"mergedRuns"+to_string(mergeCounter)).c_str());
            //Execute the hadd code (already loaded before the loop)
            string mergeFiles = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root"+'"';
            //cout << mergeFiles << endl;
            gROOT->ProcessLine(Form("hadd(%s)",mergeFiles.c_str()));
        }

        for (int i = 1; i <= nBinsBoard; i++) {
            //cout << "LB Both planes " <<  i << "\t" << hFiredBothPlanesLB->GetBinContent(i) << "\t" << hTotLB->GetBinContent(i) << endl;

            if (hTotLB->GetBinContent(i) != 0) {

                effBothLB = (hFiredBothPlanesLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effBPLB = (hFiredBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effNBPLB = (hFiredNBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;

                errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLB->GetBinContent(i));
                errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLB->GetBinContent(i));
                errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLB->GetBinContent(i));

                //Fill vector for efficiency per LB in the run
                vEffBothLB.push_back(effBothLB);
                vEffBPLB.push_back(effBPLB);
                vEffNBPLB.push_back(effNBPLB);
                
                //Fill vector for error on efficiency per LB in the run
                vErrEffBothLB.push_back(errEffBothLB);
                vErrEffBPLB.push_back(errEffBPLB);
                vErrEffNBPLB.push_back(errEffNBPLB);
            }
        }

        //Push back the vector with the eff of LB to a larger vector of vectors (one element of this = one run)
        vEffBothLB_runs.push_back(vEffBothLB);
        vEffBPLB_runs.push_back(vEffBPLB);
        vEffNBPLB_runs.push_back(vEffNBPLB);
        
        //Push back the vector with the error on eff of LB to a larger vector of vectors (one element of this = one run)
        vErrEffBothLB_runs.push_back(vErrEffBothLB);
        vErrEffBPLB_runs.push_back(vErrEffBPLB); 
        vErrEffNBPLB_runs.push_back(vErrEffNBPLB);

        //Clear vector of eff for LB in a run
        vEffBothLB.clear();
        vEffBPLB.clear();
        vEffNBPLB.clear();
        //Clear vector of error on eff for LB in a run
        vErrEffBothLB.clear();
        vErrEffBPLB.clear();
        vErrEffNBPLB.clear();

    } //End of loop on all runs

    //Close .dat file of runs to be merged
    hMergeRuns.close();

    //cout << vEffBothLB_runs[0].size() << "\t" << vEffBPLB_runs[0].size() << "\t" << vEffNBPLB_runs[0].size() << endl;

    //First [xx] is the run number in the period and the second [xx] is the LB number
    TGraphErrors *gExample = new TGraphErrors(vRun.size(),&vRun[0],&vEffBPLB_runs[0][30],NULL,&vErrEffBPLB_runs[0][30]);
    TGraphErrors *gExample2 = new TGraphErrors(vRun.size(),&vRun[0],&vEffBPLB_runs[22][30],NULL,&vErrEffBPLB_runs[22][30]);
    TGraphErrors *gExample3 = new TGraphErrors(vRun.size(),&vRun[0],&vEffBPLB_runs[30][30],NULL,&vErrEffBPLB_runs[30][30]);

    gExample->SetMarkerStyle(8);
    gExample2->SetMarkerColor(kGreen);
    gExample2->SetMarkerStyle(8);
    gExample2->SetMarkerColor(kRed);
    gExample3->SetMarkerStyle(8);
    gExample3->SetMarkerColor(kGreen);

    TMultiGraph *m = new TMultiGraph();
    m->Add(gExample);
    m->Add(gExample2);
    m->Add(gExample3);

    TCanvas *cExample = new TCanvas();
    cExample->cd();
    m->Draw("AP");
}