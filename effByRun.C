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
const int trackGoal = 8e+6; //Normal value
//const int trackGoal = 8e+7;//To test, it goes faster since this number of tracks is never reached
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
//vectors passed by reference in order to modify them when the function is called
void calculateEfficiencyByRun(string runFolder, vector<float> &effBP, vector<float> &effNBP, vector<float> &effboth,
vector<vector<float>> &effBPmerged, vector<vector<float>> &effNBPmerged, vector<vector<float>> &effbothmerged,
vector<float> &errEffBP, vector<float> &errEffNBP, vector<float> &errEffboth,
vector<vector<float>> &errEffBPmerged, vector<vector<float>> &errEffNBPmerged, vector<vector<float>> &errEffbothmerged) {

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

            effboth.push_back(effBothLB);
            effBP.push_back(effBPLB);
            effNBP.push_back(effNBPLB);

            hErrEffLB_both->Fill(errEffBothLB); 
            hErrEffLB_BP->Fill(errEffBPLB); 
            hErrEffLB_NBP->Fill(errEffNBPLB); 

            errEffboth.push_back(errEffBothLB);
            errEffBP.push_back(errEffBPLB);
            errEffNBP.push_back(errEffNBPLB);
        }
    }

    effbothmerged.push_back(effboth);
    effBPmerged.push_back(effBP);
    effNBPmerged.push_back(effNBP);

    effboth.clear();
    effBP.clear();
    effNBP.clear();

    errEffbothmerged.push_back(errEffboth);
    errEffBPmerged.push_back(errEffBP);
    errEffNBPmerged.push_back(errEffNBP);

    errEffboth.clear();
    errEffBP.clear();
    errEffNBP.clear();

    //We need the path to create an output .root file on which to save the eff and error on eff hisotgrams
    //We start from thr merged file path and we remove the "AnalysisResults.root" part at the end (20 characters)
    //Using the erase function from the string library
    runFolder.erase(runFolder.end()-20,runFolder.end());
    cout << "Test of removing string portion: " << runFolder << endl;

    //Eff on both planes
    bool makePlots = false;
    if (makePlots) {
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

        TFile *fOut = new TFile((runFolder+"mergedRuns.root").c_str(),"RECREATE");
        cEffLB_both->Write("Eff_LB_both");
        cEffLB_BP->Write("Eff_LB_BP");
        cEffLB_NBP->Write("Eff_LB_NBP");
        cErrEffLB_both->Write("Err_Eff_LB_both");
        cErrEffLB_BP->Write("Err_Eff_LB_BP");
        cErrEffLB_NBP->Write("Err_Eff_LB_NBP");
    }
    
} //End of calculateEfficieny function

void effByRun() { //Main function

    bool open = false; //used to keep track if the .dat file (updated to contain the list of .root files to be merged) has been opened once

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    //Run by run
    vector<float> vEffBothLB, vEffBPLB, vEffNBPLB;
    vector<float> vErrEffBothLB, vErrEffBPLB, vErrEffNBPLB;

    //Merged runs
    vector<float> vEffBothLBmerged, vEffBPLBmerged, vEffNBPLBmerged;
    vector<float> vErrEffBothLBmerged, vErrEffBPLBmerged, vErrEffNBPLBmerged;

    //Run by run
    vector<vector<float>> vEffBothLB_runs, vEffBPLB_runs, vEffNBPLB_runs;
    vector<vector<float>> vErrEffBothLB_runs, vErrEffBPLB_runs, vErrEffNBPLB_runs;
    
    //Merged runs
    vector<vector<float>> vEffBothLB_merged, vEffBPLB_merged, vEffNBPLB_merged;
    vector<vector<float>> vErrEffBothLB_merged, vErrEffBPLB_merged, vErrEffNBPLB_merged;

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
    //sort in ascending order
    sort(vRun.begin(), vRun.end()); 

    //Output file for the merge of root files if the number of tracks reaches the desired goal
    ofstream hMergeRuns;

    //General string name
    string fileName = "AnalysisResults.root";

    //Load hadd.C macro to merge the root files from different runs
    gROOT->ProcessLine(".L /home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/hadd.C");

    //To keep track of the number of merges and the number of times the track number was below the goal (to calculate average run number)
    int mergeCounter = 0, avgCalculation = 0;
    //To sum runs to calculate "average" run number
    float mergeRun = 0;

    //Each element of this vector is the average run number for which the merge is done
    //e.g. the merge is done between runs 340, 342 and 345 -> (340+342+345)/3=342
    vector<float> averageRun;

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

        mergeRun+=vRun.at(iRun); //sum run number to calculate "average" run number
        avgCalculation++; //increase by one for average calculations

        if (cumulativeTracks < trackGoal) { //If total track number is below the target -> Fill the file with the path of each AnalysisResults.root from each run
            //Open the output file only if it has not been opened (i.e. open == false)
            //meaning that it's the first run to be analyzed
            if (!open) {
                hMergeRuns.open("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/merged_files/runs.dat");
                open = true;
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
            hMergeRuns << runFileName << "\n";
            mergeCounter++;
            //mergeRun+=vRun.at(iRun);
            //avgCalculation++;
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
            //Call calculateEfficiency function on merged files
            calculateEfficiencyByRun(mergeFiles,vEffBPLBmerged,vEffNBPLBmerged,vEffBothLBmerged,vEffBPLB_merged,vEffNBPLB_merged,vEffBothLB_merged,
            vErrEffBPLBmerged, vErrEffNBPLBmerged, vErrEffBothLBmerged,vErrEffBPLB_merged, vErrEffNBPLB_merged, vErrEffBothLB_merged);
            averageRun.push_back(mergeRun/avgCalculation);
            cout << "avg run " << mergeRun/avgCalculation << "\t mergeRun " << mergeRun << "\t avgCalculation " << avgCalculation << endl;
            mergeRun = 0;
            avgCalculation = 0;
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
            //Test - create a sub-folder inside the merged_files directory
            gSystem->mkdir((runPath+"mergedRuns"+to_string(mergeCounter)).c_str());
            //Execute the hadd code (already loaded before the loop)
            string mergeFilesForHadd = '"'+runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root"+'"';
            string mergeFiles = runPath+"mergedRuns"+to_string(mergeCounter)+"/AnalysisResults.root";
            cout << "mergeFilesForHadd: " << mergeFilesForHadd << endl;
            gROOT->ProcessLine(Form("hadd(%s)",mergeFilesForHadd.c_str()));
            //Call calculateEfficiency function on merged files
            calculateEfficiencyByRun(mergeFiles,vEffBPLBmerged,vEffNBPLBmerged,vEffBothLBmerged,vEffBPLB_merged,vEffNBPLB_merged,vEffBothLB_merged,
            vErrEffBPLBmerged, vErrEffNBPLBmerged, vErrEffBothLBmerged,vErrEffBPLB_merged, vErrEffNBPLB_merged, vErrEffBothLB_merged);
            averageRun.push_back(mergeRun/avgCalculation);
            cout << "avg run " << mergeRun/avgCalculation << "\t mergeRun " << mergeRun << "\t avgCalculation " << avgCalculation << endl;
            mergeRun = 0;
            avgCalculation = 0;
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

    cout << "Size of the vector of vectors run by run " << vEffBothLB_runs.size() << "\t" << vEffBPLB_runs.size() << "\t" << vEffNBPLB_runs.size() << endl;
    cout << "Size of the vector of vectors merged " << vEffBothLB_merged[0].size() << "\t" << vEffBPLB_merged.size() << "\t" << vEffNBPLB_merged.size() << endl;

    cout << "size of average run: " << averageRun.size() << endl;
    for (unsigned int i = 0; i < averageRun.size(); i++) {
        cout << averageRun.at(i) << endl;
    }
    
    //Structure of the vector is the following
    // (LB1...............LB936) -> first run
    // ..
    // ..
    // ..
    // ..
    // ..
    // ..
    // ..
    // (LB1...............LB936) -> Last run   
    //In the plot we want to show columns, in order to do that we have to declare new vectors and save the data of a column inside them
    //Get efficiency for a given LB as a function of the run number
    //run by run
    vector<float> effLB1,effLB2,effLB3,effLB4;
    vector<float> errEffLB1,errEffLB2,errEffLB3,errEffLB4;

    for (unsigned int i = 0; i < vRun.size(); i++) {
        effLB1.push_back(vEffBPLB_runs[i][10]);
        effLB2.push_back(vEffBPLB_runs[i][35]);
        effLB3.push_back(vEffBPLB_runs[i][55]);
        effLB4.push_back(vEffBPLB_runs[i][145]);

        errEffLB1.push_back(vErrEffBPLB_runs[i][10]);
        errEffLB2.push_back(vErrEffBPLB_runs[i][35]);
        errEffLB3.push_back(vErrEffBPLB_runs[i][55]);
        errEffLB4.push_back(vErrEffBPLB_runs[i][145]);
    }

    //merged runs
    vector<float> effLB1merged,effLB2merged,effLB3merged,effLB4merged;
    vector<float> errEffLB1merged,errEffLB2merged,errEffLB3merged,errEffLB4merged;

    for (unsigned int i = 0; i < averageRun.size(); i++) {
        effLB1merged.push_back(vEffBPLB_merged[i][10]);
        effLB2merged.push_back(vEffBPLB_merged[i][35]);
        effLB3merged.push_back(vEffBPLB_merged[i][55]);
        effLB4merged.push_back(vEffBPLB_merged[i][145]);

        errEffLB1merged.push_back(vErrEffBPLB_merged[i][10]);
        errEffLB2merged.push_back(vErrEffBPLB_merged[i][35]);
        errEffLB3merged.push_back(vErrEffBPLB_merged[i][55]);
        errEffLB4merged.push_back(vErrEffBPLB_merged[i][145]);
    }

    TGraphErrors *gEffLBrun = new TGraphErrors(vRun.size(),&vRun[0],&effLB1[0],NULL,&errEffLB1[0]);
    TGraphErrors *gEffLBrun2 = new TGraphErrors(vRun.size(),&vRun[0],&effLB2[0],NULL,&errEffLB2[0]);
    TGraphErrors *gEffLBrun3 = new TGraphErrors(vRun.size(),&vRun[0],&effLB3[0],NULL,&errEffLB3[0]);
    TGraphErrors *gEffLBrun4 = new TGraphErrors(vRun.size(),&vRun[0],&effLB4[0],NULL,&errEffLB4[0]);

    TGraphErrors *gEffLBmerged = new TGraphErrors(averageRun.size(),&averageRun[0],&effLB1merged[0],NULL,&errEffLB1merged[0]);
    TGraphErrors *gEffLBmerged2 = new TGraphErrors(averageRun.size(),&averageRun[0],&effLB2merged[0],NULL,&errEffLB2merged[0]);
    TGraphErrors *gEffLBmerged3 = new TGraphErrors(averageRun.size(),&averageRun[0],&effLB3merged[0],NULL,&errEffLB3merged[0]);
    TGraphErrors *gEffLBmerged4 = new TGraphErrors(averageRun.size(),&averageRun[0],&effLB4merged[0],NULL,&errEffLB4merged[0]);

    gEffLBrun->SetMarkerStyle(8);
    gEffLBrun2->SetMarkerColor(kGreen);
    gEffLBrun2->SetMarkerStyle(8);
    gEffLBrun2->SetMarkerColor(kRed);
    gEffLBrun3->SetMarkerStyle(8);
    gEffLBrun3->SetMarkerColor(kGreen);
    gEffLBrun4->SetMarkerStyle(8);
    gEffLBrun4->SetMarkerColor(kMagenta);

    gEffLBmerged->SetMarkerStyle(24);
    gEffLBmerged2->SetMarkerColor(kGreen);
    gEffLBmerged2->SetMarkerStyle(24);
    gEffLBmerged2->SetMarkerColor(kRed);
    gEffLBmerged3->SetMarkerStyle(24);
    gEffLBmerged3->SetMarkerColor(kGreen);
    gEffLBmerged4->SetMarkerStyle(24);
    gEffLBmerged4->SetMarkerColor(kMagenta);

    TMultiGraph *m = new TMultiGraph();
    m->Add(gEffLBrun);
    //m->Add(gEffLBrun2);
    //m->Add(gEffLBrun3);
    //m->Add(gEffLBrun4);
    m->Add(gEffLBmerged);
    //m->Add(gEffLBmerged2);
    //m->Add(gEffLBmerged3);
    //m->Add(gEffLBmerged4);

    TCanvas *cExample = new TCanvas();
    cExample->cd();
    m->GetXaxis()->SetNoExponent(1);
    m->GetYaxis()->SetRangeUser(0,105);
    m->Draw("AP");
}