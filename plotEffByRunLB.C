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

void plotEffByRunLB() {

    float effBothLB = 0, effBPLB = 0, effNBPLB =0;
    float errEffBothLB = 0, errEffBPLB = 0, errEffNBPLB = 0;

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    //string period = "LHC22_pass7_skimmed";
    string period = "LHC23_pass4_skimmed";
    //string period = "LHC23_PbPb_pass4";
    //string period = "LHC24_pass1_skimmed";
    string globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/";

    string fileName = "";
    
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
    vector<float> vRunFloat;

    //LB225
    vector<float> vRunFloat_225_Both_before;
    vector<float> vRunFloat_225_Both_after;
    vector<float> vRunFloat_225_BP_before;
    vector<float> vRunFloat_225_BP_after;
    vector<float> vRunFloat_225_NBP_before;
    vector<float> vRunFloat_225_NBP_after;

    //LB30
    vector<float> vRunFloat_30_Both_before;
    vector<float> vRunFloat_30_Both_after;
    vector<float> vRunFloat_30_BP_before;
    vector<float> vRunFloat_30_BP_after;
    vector<float> vRunFloat_30_NBP_before;
    vector<float> vRunFloat_30_NBP_after;

    //LB64
    vector<float> vRunFloat_64_Both_before;
    vector<float> vRunFloat_64_Both_after;
    vector<float> vRunFloat_64_BP_before;
    vector<float> vRunFloat_64_BP_after;
    vector<float> vRunFloat_64_NBP_before;
    vector<float> vRunFloat_64_NBP_after;

    //LB459
    vector<float> vRunFloat_459_Both_before;
    vector<float> vRunFloat_459_Both_after;
    vector<float> vRunFloat_459_BP_before;
    vector<float> vRunFloat_459_BP_after;
    vector<float> vRunFloat_459_NBP_before;
    vector<float> vRunFloat_459_NBP_after;

    //LB264
    vector<float> vRunFloat_264_Both_before;
    vector<float> vRunFloat_264_Both_after;
    vector<float> vRunFloat_264_BP_before;
    vector<float> vRunFloat_264_BP_after;
    vector<float> vRunFloat_264_NBP_before;
    vector<float> vRunFloat_264_NBP_after;

    //LB298
    vector<float> vRunFloat_298_Both_before;
    vector<float> vRunFloat_298_Both_after;
    vector<float> vRunFloat_298_BP_before;
    vector<float> vRunFloat_298_BP_after;
    vector<float> vRunFloat_298_NBP_before;
    vector<float> vRunFloat_298_NBP_after;

    //LB693
    vector<float> vRunFloat_693_Both_before;
    vector<float> vRunFloat_693_Both_after;
    vector<float> vRunFloat_693_BP_before;
    vector<float> vRunFloat_693_BP_after;
    vector<float> vRunFloat_693_NBP_before;
    vector<float> vRunFloat_693_NBP_after;

    //LB498
    vector<float> vRunFloat_498_Both_before;
    vector<float> vRunFloat_498_Both_after;
    vector<float> vRunFloat_498_BP_before;
    vector<float> vRunFloat_498_BP_after;
    vector<float> vRunFloat_498_NBP_before;
    vector<float> vRunFloat_498_NBP_after;

    //LB532
    vector<float> vRunFloat_532_Both_before;
    vector<float> vRunFloat_532_Both_after;
    vector<float> vRunFloat_532_BP_before;
    vector<float> vRunFloat_532_BP_after;
    vector<float> vRunFloat_532_NBP_before;
    vector<float> vRunFloat_532_NBP_after;

    //LB927
    vector<float> vRunFloat_927_Both_before;
    vector<float> vRunFloat_927_Both_after;
    vector<float> vRunFloat_927_BP_before;
    vector<float> vRunFloat_927_BP_after;
    vector<float> vRunFloat_927_NBP_before;
    vector<float> vRunFloat_927_NBP_after;

    //LB732
    vector<float> vRunFloat_732_Both_before;
    vector<float> vRunFloat_732_Both_after;
    vector<float> vRunFloat_732_BP_before;
    vector<float> vRunFloat_732_BP_after;
    vector<float> vRunFloat_732_NBP_before;
    vector<float> vRunFloat_732_NBP_after;

    //LB766
    vector<float> vRunFloat_766_Both_before;
    vector<float> vRunFloat_766_Both_after;
    vector<float> vRunFloat_766_BP_before;
    vector<float> vRunFloat_766_BP_after;
    vector<float> vRunFloat_766_NBP_before;
    vector<float> vRunFloat_766_NBP_after;


    while (hRun >> isIn >> run >> IR >> bField >> start >> end){
        if (isIn) {
            vRun.push_back(run);
            vRunFloat.push_back(run);
            vStart.push_back(start);
            vEnd.push_back(end);
        }
    }

    cout << vRun.size() << "\t" << vStart.size() << "\t" << vEnd.size() << endl;

    vector<float> vEffBothLB_225_before, vEffBothLB_459_before, vEffBothLB_693_before, vEffBothLB_927_before;
    vector<float> vEffBPLB_225_before, vEffBPLB_459_before, vEffBPLB_693_before, vEffBPLB_927_before; 
    vector<float> vEffNBPLB_225_before, vEffNBPLB_459_before, vEffNBPLB_693_before, vEffNBPLB_927_before;

    vector<float> vEffBothLB_225_after, vEffBothLB_459_after, vEffBothLB_693_after, vEffBothLB_927_after;
    vector<float> vEffBPLB_225_after, vEffBPLB_459_after, vEffBPLB_693_after, vEffBPLB_927_after; 
    vector<float> vEffNBPLB_225_after, vEffNBPLB_459_after, vEffNBPLB_693_after, vEffNBPLB_927_after;

    vector<float> vEffBothLB_30_before, vEffBothLB_264_before, vEffBothLB_498_before, vEffBothLB_732_before; 
    vector<float> vEffBPLB_30_before, vEffBPLB_264_before, vEffBPLB_498_before, vEffBPLB_732_before; 
    vector<float> vEffNBPLB_30_before, vEffNBPLB_264_before, vEffNBPLB_498_before, vEffNBPLB_732_before;

    vector<float> vEffBothLB_30_after, vEffBothLB_264_after, vEffBothLB_498_after, vEffBothLB_732_after; 
    vector<float> vEffBPLB_30_after, vEffBPLB_264_after, vEffBPLB_498_after, vEffBPLB_732_after; 
    vector<float> vEffNBPLB_30_after, vEffNBPLB_264_after, vEffNBPLB_498_after, vEffNBPLB_732_after;

    vector<float> vEffBothLB_64_before, vEffBothLB_298_before, vEffBothLB_532_before, vEffBothLB_766_before; 
    vector<float> vEffBPLB_64_before, vEffBPLB_298_before, vEffBPLB_532_before, vEffBPLB_766_before; 
    vector<float> vEffNBPLB_64_before, vEffNBPLB_298_before, vEffNBPLB_532_before, vEffNBPLB_766_before;

    vector<float> vEffBothLB_64_after, vEffBothLB_298_after, vEffBothLB_532_after, vEffBothLB_766_after; 
    vector<float> vEffBPLB_64_after, vEffBPLB_298_after, vEffBPLB_532_after, vEffBPLB_766_after; 
    vector<float> vEffNBPLB_64_after, vEffNBPLB_298_after, vEffNBPLB_532_after, vEffNBPLB_766_after;

    //Tot counts per LB before
    vector<float> vTotCountsLB_225_before_both, vTotCountsLB_459_before_both, vTotCountsLB_693_before_both, vTotCountsLB_927_before_both;
    vector<float> vTotCountsLB_30_before_both, vTotCountsLB_264_before_both, vTotCountsLB_498_before_both, vTotCountsLB_732_before_both;
    vector<float> vTotCountsLB_64_before_both, vTotCountsLB_298_before_both, vTotCountsLB_532_before_both, vTotCountsLB_766_before_both;

    vector<float> vTotCountsLB_225_before_BP, vTotCountsLB_459_before_BP, vTotCountsLB_693_before_BP, vTotCountsLB_927_before_BP;
    vector<float> vTotCountsLB_30_before_BP, vTotCountsLB_264_before_BP, vTotCountsLB_498_before_BP, vTotCountsLB_732_before_BP;
    vector<float> vTotCountsLB_64_before_BP, vTotCountsLB_298_before_BP, vTotCountsLB_532_before_BP, vTotCountsLB_766_before_BP;

    vector<float> vTotCountsLB_225_before_NBP, vTotCountsLB_459_before_NBP, vTotCountsLB_693_before_NBP, vTotCountsLB_927_before_NBP;
    vector<float> vTotCountsLB_30_before_NBP, vTotCountsLB_264_before_NBP, vTotCountsLB_498_before_NBP, vTotCountsLB_732_before_NBP;
    vector<float> vTotCountsLB_64_before_NBP, vTotCountsLB_298_before_NBP, vTotCountsLB_532_before_NBP, vTotCountsLB_766_before_NBP;

    //Tot counts per LB after
    vector<float> vTotCountsLB_225_after_both, vTotCountsLB_459_after_both, vTotCountsLB_693_after_both, vTotCountsLB_927_after_both;
    vector<float> vTotCountsLB_30_after_both, vTotCountsLB_264_after_both, vTotCountsLB_498_after_both, vTotCountsLB_732_after_both;
    vector<float> vTotCountsLB_64_after_both, vTotCountsLB_298_after_both, vTotCountsLB_532_after_both, vTotCountsLB_766_after_both;

    vector<float> vTotCountsLB_225_after_BP, vTotCountsLB_459_after_BP, vTotCountsLB_693_after_BP, vTotCountsLB_927_after_BP;
    vector<float> vTotCountsLB_30_after_BP, vTotCountsLB_264_after_BP, vTotCountsLB_498_after_BP, vTotCountsLB_732_after_BP;
    vector<float> vTotCountsLB_64_after_BP, vTotCountsLB_298_after_BP, vTotCountsLB_532_after_BP, vTotCountsLB_766_after_BP;

    vector<float> vTotCountsLB_225_after_NBP, vTotCountsLB_459_after_NBP, vTotCountsLB_693_after_NBP, vTotCountsLB_927_after_NBP;
    vector<float> vTotCountsLB_30_after_NBP, vTotCountsLB_264_after_NBP, vTotCountsLB_498_after_NBP, vTotCountsLB_732_after_NBP;
    vector<float> vTotCountsLB_64_after_NBP, vTotCountsLB_298_after_NBP, vTotCountsLB_532_after_NBP, vTotCountsLB_766_after_NBP;

    //Error on eff
    vector<float> vErrEffBothLB_225_before, vErrEffBothLB_459_before, vErrEffBothLB_693_before, vErrEffBothLB_927_before;
    vector<float> vErrEffBPLB_225_before, vErrEffBPLB_459_before, vErrEffBPLB_693_before, vErrEffBPLB_927_before; 
    vector<float> vErrEffNBPLB_225_before, vErrEffNBPLB_459_before, vErrEffNBPLB_693_before, vErrEffNBPLB_927_before;

    vector<float> vErrEffBothLB_225_after, vErrEffBothLB_459_after, vErrEffBothLB_693_after, vErrEffBothLB_927_after;
    vector<float> vErrEffBPLB_225_after, vErrEffBPLB_459_after, vErrEffBPLB_693_after, vErrEffBPLB_927_after; 
    vector<float> vErrEffNBPLB_225_after, vErrEffNBPLB_459_after, vErrEffNBPLB_693_after, vErrEffNBPLB_927_after;

    vector<float> vErrEffBothLB_30_before, vErrEffBothLB_264_before, vErrEffBothLB_498_before, vErrEffBothLB_732_before; 
    vector<float> vErrEffBPLB_30_before, vErrEffBPLB_264_before, vErrEffBPLB_498_before, vErrEffBPLB_732_before; 
    vector<float> vErrEffNBPLB_30_before, vErrEffNBPLB_264_before, vErrEffNBPLB_498_before, vErrEffNBPLB_732_before;

    vector<float> vErrEffBothLB_30_after, vErrEffBothLB_264_after, vErrEffBothLB_498_after, vErrEffBothLB_732_after; 
    vector<float> vErrEffBPLB_30_after, vErrEffBPLB_264_after, vErrEffBPLB_498_after, vErrEffBPLB_732_after; 
    vector<float> vErrEffNBPLB_30_after, vErrEffNBPLB_264_after, vErrEffNBPLB_498_after, vErrEffNBPLB_732_after;

    vector<float> vErrEffBothLB_64_before, vErrEffBothLB_298_before, vErrEffBothLB_532_before, vErrEffBothLB_766_before; 
    vector<float> vErrEffBPLB_64_before, vErrEffBPLB_298_before, vErrEffBPLB_532_before, vErrEffBPLB_766_before; 
    vector<float> vErrEffNBPLB_64_before, vErrEffNBPLB_298_before, vErrEffNBPLB_532_before, vErrEffNBPLB_766_before;

    vector<float> vErrEffBothLB_64_after, vErrEffBothLB_298_after, vErrEffBothLB_532_after, vErrEffBothLB_766_after; 
    vector<float> vErrEffBPLB_64_after, vErrEffBPLB_298_after, vErrEffBPLB_532_after, vErrEffBPLB_766_after; 
    vector<float> vErrEffNBPLB_64_after, vErrEffNBPLB_298_after, vErrEffNBPLB_532_after, vErrEffNBPLB_766_after;

    //vector<int> LBs = {225,459,693,927,30,264,498,732,64,298,532,766};

    //Loop on all runs
    for (unsigned int iRun = 0; iRun < vRun.size(); iRun++) {

        fileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/" + period + "/detailedOutput/details_run_" + to_string(vRun.at(iRun)) + ".root";

        TFile *fIn = new TFile(fileName.c_str(),"READ");

        TH1F *hTotCountsBefore = (TH1F*)fIn->Get("tot");
        TH1F *hTotCountsAfter = (TH1F*)fIn->Get("totAfter");
        TH1F *hBPBefore = (TH1F*)fIn->Get("BP");
        TH1F *hBPAfter = (TH1F*)fIn->Get("BPAfter");
        TH1F *hNBPBefore = (TH1F*)fIn->Get("NBP");
        TH1F *hNBPAfter = (TH1F*)fIn->Get("NBPAfter");
        TH1F *hBothBefore = (TH1F*)fIn->Get("Both");
        TH1F *hBothAfter = (TH1F*)fIn->Get("BothAfter");

        //LB 225
        effBothLB = (hBothBefore->GetBinContent(225)/hTotCountsBefore->GetBinContent(225))*100;
        effBPLB = (hBPBefore->GetBinContent(225)/hTotCountsBefore->GetBinContent(225))*100;
        effNBPLB = (hNBPBefore->GetBinContent(225)/hTotCountsBefore->GetBinContent(225))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(225));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(225));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(225));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_225_before.push_back(effBothLB);
            vErrEffBothLB_225_before.push_back(errEffBothLB);
            vRunFloat_225_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_225_before_both.push_back(hTotCountsBefore->GetBinContent(225));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_225_before.push_back(effBPLB);
            vErrEffBPLB_225_before.push_back(errEffBPLB);
            vRunFloat_225_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_225_before_BP.push_back(hTotCountsBefore->GetBinContent(225));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_225_before.push_back(effNBPLB);
            vErrEffNBPLB_225_before.push_back(errEffNBPLB);
            vRunFloat_225_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_225_before_NBP.push_back(hTotCountsBefore->GetBinContent(225));
        }
        
        effBothLB = (hBothAfter->GetBinContent(225)/hTotCountsAfter->GetBinContent(225))*100;
        effBPLB = (hBPAfter->GetBinContent(225)/hTotCountsAfter->GetBinContent(225))*100;
        effNBPLB = (hNBPAfter->GetBinContent(225)/hTotCountsAfter->GetBinContent(225))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(225));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(225));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(225));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_225_after.push_back(effBothLB);
            vErrEffBothLB_225_after.push_back(errEffBothLB);
            vRunFloat_225_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_225_after_both.push_back(hTotCountsAfter->GetBinContent(225));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_225_after.push_back(effBPLB);
            vErrEffBPLB_225_after.push_back(errEffBPLB);
            vRunFloat_225_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_225_after_BP.push_back(hTotCountsAfter->GetBinContent(225));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_225_after.push_back(effNBPLB);
            vErrEffNBPLB_225_after.push_back(errEffNBPLB);
            vRunFloat_225_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_225_after_NBP.push_back(hTotCountsAfter->GetBinContent(225));
        }

        //LB 30
        effBothLB = (hBothBefore->GetBinContent(30)/hTotCountsBefore->GetBinContent(30))*100;
        effBPLB = (hBPBefore->GetBinContent(30)/hTotCountsBefore->GetBinContent(30))*100;
        effNBPLB = (hNBPBefore->GetBinContent(30)/hTotCountsBefore->GetBinContent(30))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(30));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(30));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(30));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_30_before.push_back(effBothLB);
            vErrEffBothLB_30_before.push_back(errEffBothLB);
            vRunFloat_30_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_30_before_both.push_back(hTotCountsBefore->GetBinContent(30));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_30_before.push_back(effBPLB);
            vErrEffBPLB_30_before.push_back(errEffBPLB);
            vRunFloat_30_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_30_before_BP.push_back(hTotCountsBefore->GetBinContent(30));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_30_before.push_back(effNBPLB);
            vErrEffNBPLB_30_before.push_back(errEffNBPLB);
            vRunFloat_30_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_30_before_NBP.push_back(hTotCountsBefore->GetBinContent(30));
        }
        
        effBothLB = (hBothAfter->GetBinContent(30)/hTotCountsAfter->GetBinContent(30))*100;
        effBPLB = (hBPAfter->GetBinContent(30)/hTotCountsAfter->GetBinContent(30))*100;
        effNBPLB = (hNBPAfter->GetBinContent(30)/hTotCountsAfter->GetBinContent(30))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(30));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(30));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(30));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_30_after.push_back(effBothLB);
            vErrEffBothLB_30_after.push_back(errEffBothLB);
            vRunFloat_30_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_30_after_both.push_back(hTotCountsAfter->GetBinContent(30));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_30_after.push_back(effBPLB);
            vErrEffBPLB_30_after.push_back(errEffBPLB);
            vRunFloat_30_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_30_after_BP.push_back(hTotCountsAfter->GetBinContent(30));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_30_after.push_back(effNBPLB);
            vErrEffNBPLB_30_after.push_back(errEffNBPLB);
            vRunFloat_30_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_30_after_NBP.push_back(hTotCountsAfter->GetBinContent(30));
        }

        //LB 64
        effBothLB = (hBothBefore->GetBinContent(64)/hTotCountsBefore->GetBinContent(64))*100;
        effBPLB = (hBPBefore->GetBinContent(64)/hTotCountsBefore->GetBinContent(64))*100;
        effNBPLB = (hNBPBefore->GetBinContent(64)/hTotCountsBefore->GetBinContent(64))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(64));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(64));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(64));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_64_before.push_back(effBothLB);
            vErrEffBothLB_64_before.push_back(errEffBothLB);
            vRunFloat_64_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_64_before_both.push_back(hTotCountsBefore->GetBinContent(64));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_64_before.push_back(effBPLB);
            vErrEffBPLB_64_before.push_back(errEffBPLB);
            vRunFloat_64_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_64_before_BP.push_back(hTotCountsBefore->GetBinContent(64));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_64_before.push_back(effNBPLB);
            vErrEffNBPLB_64_before.push_back(errEffNBPLB);
            vRunFloat_64_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_64_before_NBP.push_back(hTotCountsBefore->GetBinContent(64));
        }
        
        effBothLB = (hBothAfter->GetBinContent(64)/hTotCountsAfter->GetBinContent(64))*100;
        effBPLB = (hBPAfter->GetBinContent(64)/hTotCountsAfter->GetBinContent(64))*100;
        effNBPLB = (hNBPAfter->GetBinContent(64)/hTotCountsAfter->GetBinContent(64))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(64));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(64));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(64));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_64_after.push_back(effBothLB);
            vErrEffBothLB_64_after.push_back(errEffBothLB);
            vRunFloat_64_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_64_after_both.push_back(hTotCountsAfter->GetBinContent(64));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_64_after.push_back(effBPLB);
            vErrEffBPLB_64_after.push_back(errEffBPLB);
            vRunFloat_64_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_64_after_BP.push_back(hTotCountsAfter->GetBinContent(64));

        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_64_after.push_back(effNBPLB);
            vErrEffNBPLB_64_after.push_back(errEffNBPLB);
            vRunFloat_64_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_64_after_NBP.push_back(hTotCountsAfter->GetBinContent(64));

        }

        //LB 459
        effBothLB = (hBothBefore->GetBinContent(459)/hTotCountsBefore->GetBinContent(459))*100;
        effBPLB = (hBPBefore->GetBinContent(459)/hTotCountsBefore->GetBinContent(459))*100;
        effNBPLB = (hNBPBefore->GetBinContent(459)/hTotCountsBefore->GetBinContent(459))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(459));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(459));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(459));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_459_before.push_back(effBothLB);
            vErrEffBothLB_459_before.push_back(errEffBothLB);
            vRunFloat_459_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_459_before_both.push_back(hTotCountsBefore->GetBinContent(459));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_459_before.push_back(effBPLB);
            vErrEffBPLB_459_before.push_back(errEffBPLB);
            vRunFloat_459_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_459_before_BP.push_back(hTotCountsBefore->GetBinContent(459));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_459_before.push_back(effNBPLB);
            vErrEffNBPLB_459_before.push_back(errEffNBPLB);
            vRunFloat_459_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_459_before_NBP.push_back(hTotCountsBefore->GetBinContent(459));
        }
        
        effBothLB = (hBothAfter->GetBinContent(459)/hTotCountsAfter->GetBinContent(459))*100;
        effBPLB = (hBPAfter->GetBinContent(459)/hTotCountsAfter->GetBinContent(459))*100;
        effNBPLB = (hNBPAfter->GetBinContent(459)/hTotCountsAfter->GetBinContent(459))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(459));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(459));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(459));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_459_after.push_back(effBothLB);
            vErrEffBothLB_459_after.push_back(errEffBothLB);
            vRunFloat_459_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_459_after_both.push_back(hTotCountsAfter->GetBinContent(459));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_459_after.push_back(effBPLB);
            vErrEffBPLB_459_after.push_back(errEffBPLB);
            vRunFloat_459_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_459_after_BP.push_back(hTotCountsAfter->GetBinContent(459));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_459_after.push_back(effNBPLB);
            vErrEffNBPLB_459_after.push_back(errEffNBPLB);
            vRunFloat_459_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_459_after_NBP.push_back(hTotCountsAfter->GetBinContent(459));
        }

        //LB 264
        effBothLB = (hBothBefore->GetBinContent(264)/hTotCountsBefore->GetBinContent(264))*100;
        effBPLB = (hBPBefore->GetBinContent(264)/hTotCountsBefore->GetBinContent(264))*100;
        effNBPLB = (hNBPBefore->GetBinContent(264)/hTotCountsBefore->GetBinContent(264))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(264));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(264));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(264));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_264_before.push_back(effBothLB);
            vErrEffBothLB_264_before.push_back(errEffBothLB);
            vRunFloat_264_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_264_before_both.push_back(hTotCountsBefore->GetBinContent(264));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_264_before.push_back(effBPLB);
            vErrEffBPLB_264_before.push_back(errEffBPLB);
            vRunFloat_264_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_264_before_BP.push_back(hTotCountsBefore->GetBinContent(264));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_264_before.push_back(effNBPLB);
            vErrEffNBPLB_264_before.push_back(errEffNBPLB);
            vRunFloat_264_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_264_before_NBP.push_back(hTotCountsBefore->GetBinContent(264));
        }
        
        effBothLB = (hBothAfter->GetBinContent(264)/hTotCountsAfter->GetBinContent(264))*100;
        effBPLB = (hBPAfter->GetBinContent(264)/hTotCountsAfter->GetBinContent(264))*100;
        effNBPLB = (hNBPAfter->GetBinContent(264)/hTotCountsAfter->GetBinContent(264))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(264));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(264));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(264));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_264_after.push_back(effBothLB);
            vErrEffBothLB_264_after.push_back(errEffBothLB);
            vRunFloat_264_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_264_after_both.push_back(hTotCountsAfter->GetBinContent(264));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_264_after.push_back(effBPLB);
            vErrEffBPLB_264_after.push_back(errEffBPLB);
            vRunFloat_264_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_264_after_BP.push_back(hTotCountsAfter->GetBinContent(264));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_264_after.push_back(effNBPLB);
            vErrEffNBPLB_264_after.push_back(errEffNBPLB);
            vRunFloat_264_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_264_after_NBP.push_back(hTotCountsAfter->GetBinContent(264));
        }

        //LB 298
        effBothLB = (hBothBefore->GetBinContent(298)/hTotCountsBefore->GetBinContent(298))*100;
        effBPLB = (hBPBefore->GetBinContent(298)/hTotCountsBefore->GetBinContent(298))*100;
        effNBPLB = (hNBPBefore->GetBinContent(298)/hTotCountsBefore->GetBinContent(298))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(298));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(298));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(298));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_298_before.push_back(effBothLB);
            vErrEffBothLB_298_before.push_back(errEffBothLB);
            vRunFloat_298_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_298_before_both.push_back(hTotCountsBefore->GetBinContent(298));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_298_before.push_back(effBPLB);
            vErrEffBPLB_298_before.push_back(errEffBPLB);
            vRunFloat_298_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_298_before_BP.push_back(hTotCountsBefore->GetBinContent(298));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_298_before.push_back(effNBPLB);
            vErrEffNBPLB_298_before.push_back(errEffNBPLB);
            vRunFloat_298_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_298_before_NBP.push_back(hTotCountsBefore->GetBinContent(298));
        }
        
        effBothLB = (hBothAfter->GetBinContent(298)/hTotCountsAfter->GetBinContent(298))*100;
        effBPLB = (hBPAfter->GetBinContent(298)/hTotCountsAfter->GetBinContent(298))*100;
        effNBPLB = (hNBPAfter->GetBinContent(298)/hTotCountsAfter->GetBinContent(298))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(298));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(298));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(298));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_298_after.push_back(effBothLB);
            vErrEffBothLB_298_after.push_back(errEffBothLB);
            vRunFloat_298_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_298_after_both.push_back(hTotCountsAfter->GetBinContent(298));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_298_after.push_back(effBPLB);
            vErrEffBPLB_298_after.push_back(errEffBPLB);
            vRunFloat_298_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_298_after_BP.push_back(hTotCountsAfter->GetBinContent(298));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_298_after.push_back(effNBPLB);
            vErrEffNBPLB_298_after.push_back(errEffNBPLB);
            vRunFloat_298_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_298_after_NBP.push_back(hTotCountsAfter->GetBinContent(298));
        }

        //LB 693
        effBothLB = (hBothBefore->GetBinContent(693)/hTotCountsBefore->GetBinContent(693))*100;
        effBPLB = (hBPBefore->GetBinContent(693)/hTotCountsBefore->GetBinContent(693))*100;
        effNBPLB = (hNBPBefore->GetBinContent(693)/hTotCountsBefore->GetBinContent(693))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(693));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(693));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(693));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_693_before.push_back(effBothLB);
            vErrEffBothLB_693_before.push_back(errEffBothLB);
            vRunFloat_693_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_693_before_both.push_back(hTotCountsBefore->GetBinContent(693));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_693_before.push_back(effBPLB);
            vErrEffBPLB_693_before.push_back(errEffBPLB);
            vRunFloat_693_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_693_before_BP.push_back(hTotCountsBefore->GetBinContent(693));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_693_before.push_back(effNBPLB);
            vErrEffNBPLB_693_before.push_back(errEffNBPLB);
            vRunFloat_693_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_693_before_NBP.push_back(hTotCountsBefore->GetBinContent(693));
        }
        
        effBothLB = (hBothAfter->GetBinContent(693)/hTotCountsAfter->GetBinContent(693))*100;
        effBPLB = (hBPAfter->GetBinContent(693)/hTotCountsAfter->GetBinContent(693))*100;
        effNBPLB = (hNBPAfter->GetBinContent(693)/hTotCountsAfter->GetBinContent(693))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(693));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(693));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(693));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_693_after.push_back(effBothLB);
            vErrEffBothLB_693_after.push_back(errEffBothLB);
            vRunFloat_693_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_693_after_both.push_back(hTotCountsAfter->GetBinContent(693));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_693_after.push_back(effBPLB);
            vErrEffBPLB_693_after.push_back(errEffBPLB);
            vRunFloat_693_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_693_after_BP.push_back(hTotCountsAfter->GetBinContent(693));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_693_after.push_back(effNBPLB);
            vErrEffNBPLB_693_after.push_back(errEffNBPLB);
            vRunFloat_693_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_693_after_NBP.push_back(hTotCountsAfter->GetBinContent(693));
        }

        //LB 498
        effBothLB = (hBothBefore->GetBinContent(498)/hTotCountsBefore->GetBinContent(498))*100;
        effBPLB = (hBPBefore->GetBinContent(498)/hTotCountsBefore->GetBinContent(498))*100;
        effNBPLB = (hNBPBefore->GetBinContent(498)/hTotCountsBefore->GetBinContent(498))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(498));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(498));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(498));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_498_before.push_back(effBothLB);
            vErrEffBothLB_498_before.push_back(errEffBothLB);
            vRunFloat_498_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_498_before_both.push_back(hTotCountsBefore->GetBinContent(498));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_498_before.push_back(effBPLB);
            vErrEffBPLB_498_before.push_back(errEffBPLB);
            vRunFloat_498_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_498_before_BP.push_back(hTotCountsBefore->GetBinContent(498));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_498_before.push_back(effNBPLB);
            vErrEffNBPLB_498_before.push_back(errEffNBPLB);
            vRunFloat_498_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_498_before_NBP.push_back(hTotCountsBefore->GetBinContent(498));
        }
        
        effBothLB = (hBothAfter->GetBinContent(498)/hTotCountsAfter->GetBinContent(498))*100;
        effBPLB = (hBPAfter->GetBinContent(498)/hTotCountsAfter->GetBinContent(498))*100;
        effNBPLB = (hNBPAfter->GetBinContent(498)/hTotCountsAfter->GetBinContent(498))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(498));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(498));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(498));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_498_after.push_back(effBothLB);
            vErrEffBothLB_498_after.push_back(errEffBothLB);
            vRunFloat_498_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_498_after_both.push_back(hTotCountsAfter->GetBinContent(498));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_498_after.push_back(effBPLB);
            vErrEffBPLB_498_after.push_back(errEffBPLB);
            vRunFloat_498_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_498_after_BP.push_back(hTotCountsAfter->GetBinContent(498));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_498_after.push_back(effNBPLB);
            vErrEffNBPLB_498_after.push_back(errEffNBPLB);
            vRunFloat_498_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_498_after_NBP.push_back(hTotCountsAfter->GetBinContent(498));
        }

        //LB 532
        effBothLB = (hBothBefore->GetBinContent(532)/hTotCountsBefore->GetBinContent(532))*100;
        effBPLB = (hBPBefore->GetBinContent(532)/hTotCountsBefore->GetBinContent(532))*100;
        effNBPLB = (hNBPBefore->GetBinContent(532)/hTotCountsBefore->GetBinContent(532))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(532));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(532));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(532));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_532_before.push_back(effBothLB);
            vErrEffBothLB_532_before.push_back(errEffBothLB);
            vRunFloat_532_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_532_before_both.push_back(hTotCountsBefore->GetBinContent(532));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_532_before.push_back(effBPLB);
            vErrEffBPLB_532_before.push_back(errEffBPLB);
            vRunFloat_532_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_532_before_BP.push_back(hTotCountsBefore->GetBinContent(532));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_532_before.push_back(effNBPLB);
            vErrEffNBPLB_532_before.push_back(errEffNBPLB);
            vRunFloat_532_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_532_before_NBP.push_back(hTotCountsBefore->GetBinContent(532));
        }
        
        effBothLB = (hBothAfter->GetBinContent(532)/hTotCountsAfter->GetBinContent(532))*100;
        effBPLB = (hBPAfter->GetBinContent(532)/hTotCountsAfter->GetBinContent(532))*100;
        effNBPLB = (hNBPAfter->GetBinContent(532)/hTotCountsAfter->GetBinContent(532))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(532));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(532));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(532));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_532_after.push_back(effBothLB);
            vErrEffBothLB_532_after.push_back(errEffBothLB);
            vRunFloat_532_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_532_after_both.push_back(hTotCountsAfter->GetBinContent(532));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_532_after.push_back(effBPLB);
            vErrEffBPLB_532_after.push_back(errEffBPLB);
            vRunFloat_532_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_532_after_BP.push_back(hTotCountsAfter->GetBinContent(532));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_532_after.push_back(effNBPLB);
            vErrEffNBPLB_532_after.push_back(errEffNBPLB);
            vRunFloat_532_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_532_after_NBP.push_back(hTotCountsAfter->GetBinContent(532));
        }

        //LB 927
        effBothLB = (hBothBefore->GetBinContent(927)/hTotCountsBefore->GetBinContent(927))*100;
        effBPLB = (hBPBefore->GetBinContent(927)/hTotCountsBefore->GetBinContent(927))*100;
        effNBPLB = (hNBPBefore->GetBinContent(927)/hTotCountsBefore->GetBinContent(927))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(927));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(927));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(927));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_927_before.push_back(effBothLB);
            vErrEffBothLB_927_before.push_back(errEffBothLB);
            vRunFloat_927_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_927_before_both.push_back(hTotCountsBefore->GetBinContent(927));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_927_before.push_back(effBPLB);
            vErrEffBPLB_927_before.push_back(errEffBPLB);
            vRunFloat_927_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_927_before_BP.push_back(hTotCountsBefore->GetBinContent(927));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_927_before.push_back(effNBPLB);
            vErrEffNBPLB_927_before.push_back(errEffNBPLB);
            vRunFloat_927_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_927_before_NBP.push_back(hTotCountsBefore->GetBinContent(927));
        }
        
        effBothLB = (hBothAfter->GetBinContent(927)/hTotCountsAfter->GetBinContent(927))*100;
        effBPLB = (hBPAfter->GetBinContent(927)/hTotCountsAfter->GetBinContent(927))*100;
        effNBPLB = (hNBPAfter->GetBinContent(927)/hTotCountsAfter->GetBinContent(927))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(927));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(927));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(927));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_927_after.push_back(effBothLB);
            vErrEffBothLB_927_after.push_back(errEffBothLB);
            vRunFloat_927_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_927_after_both.push_back(hTotCountsAfter->GetBinContent(927));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_927_after.push_back(effBPLB);
            vErrEffBPLB_927_after.push_back(errEffBPLB);
            vRunFloat_927_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_927_after_BP.push_back(hTotCountsAfter->GetBinContent(927));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_927_after.push_back(effNBPLB);
            vErrEffNBPLB_927_after.push_back(errEffNBPLB);
            vRunFloat_927_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_927_after_NBP.push_back(hTotCountsAfter->GetBinContent(927));
        }

        //LB 732
        effBothLB = (hBothBefore->GetBinContent(732)/hTotCountsBefore->GetBinContent(732))*100;
        effBPLB = (hBPBefore->GetBinContent(732)/hTotCountsBefore->GetBinContent(732))*100;
        effNBPLB = (hNBPBefore->GetBinContent(732)/hTotCountsBefore->GetBinContent(732))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(732));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(732));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(732));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_732_before.push_back(effBothLB);
            vErrEffBothLB_732_before.push_back(errEffBothLB);
            vRunFloat_732_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_732_before_both.push_back(hTotCountsBefore->GetBinContent(732));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_732_before.push_back(effBPLB);
            vErrEffBPLB_732_before.push_back(errEffBPLB);
            vRunFloat_732_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_732_before_BP.push_back(hTotCountsBefore->GetBinContent(732));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_732_before.push_back(effNBPLB);
            vErrEffNBPLB_732_before.push_back(errEffNBPLB);
            vRunFloat_732_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_732_before_NBP.push_back(hTotCountsBefore->GetBinContent(732));
        }
        
        effBothLB = (hBothAfter->GetBinContent(732)/hTotCountsAfter->GetBinContent(732))*100;
        effBPLB = (hBPAfter->GetBinContent(732)/hTotCountsAfter->GetBinContent(732))*100;
        effNBPLB = (hNBPAfter->GetBinContent(732)/hTotCountsAfter->GetBinContent(732))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(732));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(732));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(732));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_732_after.push_back(effBothLB);
            vErrEffBothLB_732_after.push_back(errEffBothLB);
            vRunFloat_732_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_732_after_both.push_back(hTotCountsAfter->GetBinContent(732));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_732_after.push_back(effBPLB);
            vErrEffBPLB_732_after.push_back(errEffBPLB);
            vRunFloat_732_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_732_after_BP.push_back(hTotCountsAfter->GetBinContent(732));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_732_after.push_back(effNBPLB);
            vErrEffNBPLB_732_after.push_back(errEffNBPLB);
            vRunFloat_732_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_732_after_NBP.push_back(hTotCountsAfter->GetBinContent(732));
        }

        //LB 766
        effBothLB = (hBothBefore->GetBinContent(766)/hTotCountsBefore->GetBinContent(766))*100;
        effBPLB = (hBPBefore->GetBinContent(766)/hTotCountsBefore->GetBinContent(766))*100;
        effNBPLB = (hNBPBefore->GetBinContent(766)/hTotCountsBefore->GetBinContent(766))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsBefore->GetBinContent(766));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsBefore->GetBinContent(766));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsBefore->GetBinContent(766));

        cout << "before " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "before " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "before " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_766_before.push_back(effBothLB);
            vErrEffBothLB_766_before.push_back(errEffBothLB);
            vRunFloat_766_Both_before.push_back(vRun.at(iRun));
            vTotCountsLB_766_before_both.push_back(hTotCountsBefore->GetBinContent(766));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_766_before.push_back(effBPLB);
            vErrEffBPLB_766_before.push_back(errEffBPLB);
            vRunFloat_766_BP_before.push_back(vRun.at(iRun));
            vTotCountsLB_766_before_BP.push_back(hTotCountsBefore->GetBinContent(766));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_766_before.push_back(effNBPLB);
            vErrEffNBPLB_766_before.push_back(errEffNBPLB);
            vRunFloat_766_NBP_before.push_back(vRun.at(iRun));
            vTotCountsLB_766_before_NBP.push_back(hTotCountsBefore->GetBinContent(766));
        }
        
        effBothLB = (hBothAfter->GetBinContent(766)/hTotCountsAfter->GetBinContent(766))*100;
        effBPLB = (hBPAfter->GetBinContent(766)/hTotCountsAfter->GetBinContent(766))*100;
        effNBPLB = (hNBPAfter->GetBinContent(766)/hTotCountsAfter->GetBinContent(766))*100;
        errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotCountsAfter->GetBinContent(766));
        errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotCountsAfter->GetBinContent(766));
        errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotCountsAfter->GetBinContent(766));

        cout << "after " << effBothLB << "+-" << errEffBothLB << endl;
        cout << "after " << effBPLB << "+-" << errEffBPLB << endl;
        cout << "after " << effNBPLB << "+-" << errEffNBPLB << endl;

        if (!std::isnan(effBothLB)) {
            vEffBothLB_766_after.push_back(effBothLB);
            vErrEffBothLB_766_after.push_back(errEffBothLB);
            vRunFloat_766_Both_after.push_back(vRun.at(iRun));
            vTotCountsLB_766_after_both.push_back(hTotCountsAfter->GetBinContent(766));
        }

        if (!std::isnan(effBPLB)) {
            vEffBPLB_766_after.push_back(effBPLB);
            vErrEffBPLB_766_after.push_back(errEffBPLB);
            vRunFloat_766_BP_after.push_back(vRun.at(iRun));
            vTotCountsLB_766_after_BP.push_back(hTotCountsAfter->GetBinContent(766));
        }

        if (!std::isnan(effNBPLB)) {
            vEffNBPLB_766_after.push_back(effNBPLB);
            vErrEffNBPLB_766_after.push_back(errEffNBPLB);
            vRunFloat_766_NBP_after.push_back(vRun.at(iRun));
            vTotCountsLB_766_after_NBP.push_back(hTotCountsAfter->GetBinContent(766));
        }
        

        fIn->Close();
        
    } //End of loop on runs

    TLegend *lLB225 = new TLegend(0.65,0.12,0.85,0.37,"","rNDC");
    lLB225->SetBorderSize(0); //No borders in legend
	lLB225->SetFillStyle(0); //Transparent background
	lLB225->SetTextFont(62);

    TLegend *lLB30 = new TLegend(0.65,0.12,0.85,0.37,"","rNDC");
    lLB30->SetBorderSize(0); //No borders in legend
	lLB30->SetFillStyle(0); //Transparent background
	lLB30->SetTextFont(62);

    TLegend *lLB64 = new TLegend(0.65,0.12,0.85,0.37,"","rNDC");
    lLB64->SetBorderSize(0); //No borders in legend
	lLB64->SetFillStyle(0); //Transparent background
	lLB64->SetTextFont(62);


    //Both planes MT11: 225, 30, 64
    TGraphErrors *gEffLB_225_Both_before = new TGraphErrors(vRunFloat_225_Both_before.size(),&vRunFloat_225_Both_before[0],&vEffBothLB_225_before[0],NULL,&vErrEffBothLB_225_before[0]); //MT11
    gEffLB_225_Both_before->SetName("gEffLB_225_Both_before");
    gEffLB_225_Both_before->SetTitle("gEffLB_225_Both_before");
    gEffLB_225_Both_before->SetMarkerSize(1.4);
    gEffLB_225_Both_before->SetMarkerColor(kBlack);
    gEffLB_225_Both_before->SetMarkerStyle(8);
    lLB225->AddEntry(gEffLB_225_Both_before,"Before p_{t} cut","p");
    TGraphErrors *gEffLB_225_Both_after = new TGraphErrors(vRunFloat_225_Both_after.size(),&vRunFloat_225_Both_after[0],&vEffBothLB_225_after[0],NULL,&vErrEffBothLB_225_after[0]); //MT11
    gEffLB_225_Both_after->SetName("gEffLB_225_Both_before");
    gEffLB_225_Both_after->SetTitle("gEffLB_225_Both_before");
    gEffLB_225_Both_after->SetMarkerSize(1.4);
    gEffLB_225_Both_after->SetMarkerColor(kRed);
    gEffLB_225_Both_after->SetMarkerStyle(24);
    lLB225->AddEntry(gEffLB_225_Both_after,"After p_{t} cut","p");
    //
    TGraphErrors *gEffLB_30_Both_before = new TGraphErrors(vRunFloat_30_Both_before.size(),&vRunFloat_30_Both_before[0],&vEffBothLB_30_before[0],NULL,&vErrEffBothLB_30_before[0]); //MT11
    gEffLB_30_Both_before->SetName("gEffLB_30_Both_before");
    gEffLB_30_Both_before->SetTitle("gEffLB_30_Both_before");
    gEffLB_30_Both_before->SetMarkerSize(1.4);
    gEffLB_30_Both_before->SetMarkerColor(kGreen+3);
    gEffLB_30_Both_before->SetMarkerStyle(8);
    lLB30->AddEntry(gEffLB_30_Both_before,"Before p_{t} cut","p");
    TGraphErrors *gEffLB_30_Both_after = new TGraphErrors(vRunFloat_30_Both_after.size(),&vRunFloat_30_Both_after[0],&vEffBothLB_30_after[0],NULL,&vErrEffBothLB_30_after[0]); //MT11
    gEffLB_30_Both_after->SetName("gEffLB_30_Both_before");
    gEffLB_30_Both_after->SetTitle("gEffLB_30_Both_before");
    gEffLB_30_Both_after->SetMarkerSize(1.4);
    gEffLB_30_Both_after->SetMarkerColor(kBlue);
    gEffLB_30_Both_after->SetMarkerStyle(24);
    lLB30->AddEntry(gEffLB_30_Both_after,"After p_{t} cut","p");
    //
    TGraphErrors *gEffLB_64_Both_before = new TGraphErrors(vRunFloat_64_Both_before.size(),&vRunFloat_64_Both_before[0],&vEffBothLB_64_before[0],NULL,&vErrEffBothLB_64_before[0]); //MT11
    gEffLB_64_Both_before->SetName("gEffLB_64_Both_before");
    gEffLB_64_Both_before->SetTitle("gEffLB_64_Both_before");
    gEffLB_64_Both_before->SetMarkerSize(1.4);
    gEffLB_64_Both_before->SetMarkerColor(kMagenta);
    gEffLB_64_Both_before->SetMarkerStyle(8);
    lLB64->AddEntry(gEffLB_64_Both_before,"Before p_{t} cut","p");
    TGraphErrors *gEffLB_64_Both_after = new TGraphErrors(vRunFloat_64_Both_after.size(),&vRunFloat_64_Both_after[0],&vEffBothLB_64_after[0],NULL,&vErrEffBothLB_64_after[0]); //MT11
    gEffLB_64_Both_after->SetName("gEffLB_64_Both_before");
    gEffLB_64_Both_after->SetTitle("gEffLB_64_Both_before");
    gEffLB_64_Both_after->SetMarkerSize(1.4);
    gEffLB_64_Both_after->SetMarkerColor(kOrange);
    gEffLB_64_Both_after->SetMarkerStyle(24);
    lLB64->AddEntry(gEffLB_64_Both_after,"After p_{t} cut","p");

    //Both planes MT12 459, 264, 298
    TGraphErrors *gEffLB_459_Both_before = new TGraphErrors(vRunFloat_459_Both_before.size(),&vRunFloat_459_Both_before[0],&vEffBothLB_459_before[0],NULL,&vErrEffBothLB_459_before[0]); //MT11
    gEffLB_459_Both_before->SetName("gEffLB_459_Both_before");
    gEffLB_459_Both_before->SetTitle("gEffLB_459_Both_before");
    gEffLB_459_Both_before->SetMarkerSize(1.4);
    gEffLB_459_Both_before->SetMarkerColor(kBlack);
    gEffLB_459_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_459_Both_after = new TGraphErrors(vRunFloat_459_Both_after.size(),&vRunFloat_459_Both_after[0],&vEffBothLB_459_after[0],NULL,&vErrEffBothLB_459_after[0]); //MT11
    gEffLB_459_Both_after->SetName("gEffLB_459_Both_after");
    gEffLB_459_Both_after->SetTitle("gEffLB_459_Both_after");
    gEffLB_459_Both_after->SetMarkerSize(1.4);
    gEffLB_459_Both_after->SetMarkerColor(kRed);
    gEffLB_459_Both_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_264_Both_before = new TGraphErrors(vRunFloat_264_Both_before.size(),&vRunFloat_264_Both_before[0],&vEffBothLB_264_before[0],NULL,&vErrEffBothLB_264_before[0]); //MT11
    gEffLB_264_Both_before->SetName("gEffLB_264_Both_before");
    gEffLB_264_Both_before->SetTitle("gEffLB_264_Both_before");
    gEffLB_264_Both_before->SetMarkerSize(1.4);
    gEffLB_264_Both_before->SetMarkerColor(kGreen+3);
    gEffLB_264_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_264_Both_after = new TGraphErrors(vRunFloat_264_Both_after.size(),&vRunFloat_264_Both_after[0],&vEffBothLB_264_after[0],NULL,&vErrEffBothLB_264_after[0]); //MT11
    gEffLB_264_Both_after->SetName("gEffLB_264_Both_before");
    gEffLB_264_Both_after->SetTitle("gEffLB_264_Both_before");
    gEffLB_264_Both_after->SetMarkerSize(1.4);
    gEffLB_264_Both_after->SetMarkerColor(kBlue);
    gEffLB_264_Both_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_298_Both_before = new TGraphErrors(vRunFloat_298_Both_before.size(),&vRunFloat_298_Both_before[0],&vEffBothLB_298_before[0],NULL,&vErrEffBothLB_298_before[0]); //MT11
    gEffLB_298_Both_before->SetName("gEffLB_298_Both_before");
    gEffLB_298_Both_before->SetTitle("gEffLB_298_Both_before");
    gEffLB_298_Both_before->SetMarkerSize(1.4);
    gEffLB_298_Both_before->SetMarkerColor(kMagenta);
    gEffLB_298_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_298_Both_after = new TGraphErrors(vRunFloat_298_Both_after.size(),&vRunFloat_298_Both_after[0],&vEffBothLB_298_after[0],NULL,&vErrEffBothLB_298_after[0]); //MT11
    gEffLB_298_Both_after->SetName("gEffLB_298_Both_before");
    gEffLB_298_Both_after->SetTitle("gEffLB_298_Both_before");
    gEffLB_298_Both_after->SetMarkerSize(1.4);
    gEffLB_298_Both_after->SetMarkerColor(kOrange);
    gEffLB_298_Both_after->SetMarkerStyle(24);

    //Both planes MT21 693, 498, 532
    TGraphErrors *gEffLB_693_Both_before = new TGraphErrors(vRunFloat_693_Both_before.size(),&vRunFloat_693_Both_before[0],&vEffBothLB_693_before[0],NULL,&vErrEffBothLB_693_before[0]); //MT11
    gEffLB_693_Both_before->SetName("gEffLB_693_Both_before");
    gEffLB_693_Both_before->SetTitle("gEffLB_693_Both_before");
    gEffLB_693_Both_before->SetMarkerSize(1.4);
    gEffLB_693_Both_before->SetMarkerColor(kBlack);
    gEffLB_693_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_693_Both_after = new TGraphErrors(vRunFloat_693_Both_after.size(),&vRunFloat_693_Both_after[0],&vEffBothLB_693_after[0],NULL,&vErrEffBothLB_693_after[0]); //MT11
    gEffLB_693_Both_after->SetName("gEffLB_693_Both_after");
    gEffLB_693_Both_after->SetTitle("gEffLB_693_Both_after");
    gEffLB_693_Both_after->SetMarkerSize(1.4);
    gEffLB_693_Both_after->SetMarkerColor(kRed);
    gEffLB_693_Both_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_498_Both_before = new TGraphErrors(vRunFloat_498_Both_before.size(),&vRunFloat_498_Both_before[0],&vEffBothLB_498_before[0],NULL,&vErrEffBothLB_498_before[0]); //MT11
    gEffLB_498_Both_before->SetName("gEffLB_498_Both_before");
    gEffLB_498_Both_before->SetTitle("gEffLB_498_Both_before");
    gEffLB_498_Both_before->SetMarkerSize(1.4);
    gEffLB_498_Both_before->SetMarkerColor(kGreen+3);
    gEffLB_498_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_498_Both_after = new TGraphErrors(vRunFloat_498_Both_after.size(),&vRunFloat_498_Both_after[0],&vEffBothLB_498_after[0],NULL,&vErrEffBothLB_498_after[0]); //MT11
    gEffLB_498_Both_after->SetName("gEffLB_498_Both_before");
    gEffLB_498_Both_after->SetTitle("gEffLB_498_Both_before");
    gEffLB_498_Both_after->SetMarkerSize(1.4);
    gEffLB_498_Both_after->SetMarkerColor(kBlue);
    gEffLB_498_Both_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_532_Both_before = new TGraphErrors(vRunFloat_532_Both_before.size(),&vRunFloat_532_Both_before[0],&vEffBothLB_532_before[0],NULL,&vErrEffBothLB_532_before[0]); //MT11
    gEffLB_532_Both_before->SetName("gEffLB_532_Both_before");
    gEffLB_532_Both_before->SetTitle("gEffLB_532_Both_before");
    gEffLB_532_Both_before->SetMarkerSize(1.4);
    gEffLB_532_Both_before->SetMarkerColor(kMagenta);
    gEffLB_532_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_532_Both_after = new TGraphErrors(vRunFloat_532_Both_after.size(),&vRunFloat_532_Both_after[0],&vEffBothLB_532_after[0],NULL,&vErrEffBothLB_532_after[0]); //MT11
    gEffLB_532_Both_after->SetName("gEffLB_532_Both_before");
    gEffLB_532_Both_after->SetTitle("gEffLB_532_Both_before");
    gEffLB_532_Both_after->SetMarkerSize(1.4);
    gEffLB_532_Both_after->SetMarkerColor(kOrange);
    gEffLB_532_Both_after->SetMarkerStyle(24);

    //Both planes MT22 927, 732, 766
    TGraphErrors *gEffLB_927_Both_before = new TGraphErrors(vRunFloat_927_Both_before.size(),&vRunFloat_927_Both_before[0],&vEffBothLB_927_before[0],NULL,&vErrEffBothLB_927_before[0]); //MT11
    gEffLB_927_Both_before->SetName("gEffLB_927_Both_before");
    gEffLB_927_Both_before->SetTitle("gEffLB_927_Both_before");
    gEffLB_927_Both_before->SetMarkerSize(1.4);
    gEffLB_927_Both_before->SetMarkerColor(kBlack);
    gEffLB_927_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_927_Both_after = new TGraphErrors(vRunFloat_927_Both_after.size(),&vRunFloat_927_Both_after[0],&vEffBothLB_927_after[0],NULL,&vErrEffBothLB_927_after[0]); //MT11
    gEffLB_927_Both_after->SetName("gEffLB_927_Both_after");
    gEffLB_927_Both_after->SetTitle("gEffLB_927_Both_after");
    gEffLB_927_Both_after->SetMarkerSize(1.4);
    gEffLB_927_Both_after->SetMarkerColor(kRed);
    gEffLB_927_Both_after->SetMarkerStyle(24);

    //
    TGraphErrors *gEffLB_732_Both_before = new TGraphErrors(vRunFloat_732_Both_before.size(),&vRunFloat_732_Both_before[0],&vEffBothLB_732_before[0],NULL,&vErrEffBothLB_732_before[0]); //MT11
    gEffLB_732_Both_before->SetName("gEffLB_732_Both_before");
    gEffLB_732_Both_before->SetTitle("gEffLB_732_Both_before");
    gEffLB_732_Both_before->SetMarkerSize(1.4);
    gEffLB_732_Both_before->SetMarkerColor(kGreen+3);
    gEffLB_732_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_732_Both_after = new TGraphErrors(vRunFloat_732_Both_after.size(),&vRunFloat_732_Both_after[0],&vEffBothLB_732_after[0],NULL,&vErrEffBothLB_732_after[0]); //MT11
    gEffLB_732_Both_after->SetName("gEffLB_732_Both_before");
    gEffLB_732_Both_after->SetTitle("gEffLB_732_Both_before");
    gEffLB_732_Both_after->SetMarkerSize(1.4);
    gEffLB_732_Both_after->SetMarkerColor(kBlue);
    gEffLB_732_Both_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_766_Both_before = new TGraphErrors(vRunFloat_766_Both_before.size(),&vRunFloat_766_Both_before[0],&vEffBothLB_766_before[0],NULL,&vErrEffBothLB_766_before[0]); //MT11
    gEffLB_766_Both_before->SetName("gEffLB_766_Both_before");
    gEffLB_766_Both_before->SetTitle("gEffLB_766_Both_before");
    gEffLB_766_Both_before->SetMarkerSize(1.4);
    gEffLB_766_Both_before->SetMarkerColor(kMagenta);
    gEffLB_766_Both_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_766_Both_after = new TGraphErrors(vRunFloat_766_Both_after.size(),&vRunFloat_766_Both_after[0],&vEffBothLB_766_after[0],NULL,&vErrEffBothLB_766_after[0]); //MT11
    gEffLB_766_Both_after->SetName("gEffLB_766_Both_before");
    gEffLB_766_Both_after->SetTitle("gEffLB_766_Both_before");
    gEffLB_766_Both_after->SetMarkerSize(1.4);
    gEffLB_766_Both_after->SetMarkerColor(kOrange);
    gEffLB_766_Both_after->SetMarkerStyle(24);

    //BP MT11: 225, 30, 64
    TGraphErrors *gEffLB_225_BP_before = new TGraphErrors(vRunFloat_225_BP_before.size(),&vRunFloat_225_BP_before[0],&vEffBPLB_225_before[0],NULL,&vErrEffBPLB_225_before[0]); //MT11
    gEffLB_225_BP_before->SetName("gEffLB_225_BP_before");
    gEffLB_225_BP_before->SetTitle("gEffLB_225_BP_before");
    gEffLB_225_BP_before->SetMarkerSize(1.4);
    gEffLB_225_BP_before->SetMarkerColor(kBlack);
    gEffLB_225_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_225_BP_after = new TGraphErrors(vRunFloat_225_BP_after.size(),&vRunFloat_225_BP_after[0],&vEffBPLB_225_after[0],NULL,&vErrEffBPLB_225_after[0]); //MT11
    gEffLB_225_BP_after->SetName("gEffLB_225_BP_before");
    gEffLB_225_BP_after->SetTitle("gEffLB_225_BP_before");
    gEffLB_225_BP_after->SetMarkerSize(1.4);
    gEffLB_225_BP_after->SetMarkerColor(kRed);
    gEffLB_225_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_30_BP_before = new TGraphErrors(vRunFloat_30_BP_before.size(),&vRunFloat_30_BP_before[0],&vEffBPLB_30_before[0],NULL,&vErrEffBPLB_30_before[0]); //MT11
    gEffLB_30_BP_before->SetName("gEffLB_30_BP_before");
    gEffLB_30_BP_before->SetTitle("gEffLB_30_BP_before");
    gEffLB_30_BP_before->SetMarkerSize(1.4);
    gEffLB_30_BP_before->SetMarkerColor(kGreen+3);
    gEffLB_30_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_30_BP_after = new TGraphErrors(vRunFloat_30_BP_after.size(),&vRunFloat_30_BP_after[0],&vEffBPLB_30_after[0],NULL,&vErrEffBPLB_30_after[0]); //MT11
    gEffLB_30_BP_after->SetName("gEffLB_30_BP_before");
    gEffLB_30_BP_after->SetTitle("gEffLB_30_BP_before");
    gEffLB_30_BP_after->SetMarkerSize(1.4);
    gEffLB_30_BP_after->SetMarkerColor(kBlue);
    gEffLB_30_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_64_BP_before = new TGraphErrors(vRunFloat_64_BP_before.size(),&vRunFloat_64_BP_before[0],&vEffBPLB_64_before[0],NULL,&vErrEffBPLB_64_before[0]); //MT11
    gEffLB_64_BP_before->SetName("gEffLB_64_BP_before");
    gEffLB_64_BP_before->SetTitle("gEffLB_64_BP_before");
    gEffLB_64_BP_before->SetMarkerSize(1.4);
    gEffLB_64_BP_before->SetMarkerColor(kMagenta);
    gEffLB_64_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_64_BP_after = new TGraphErrors(vRunFloat_64_BP_after.size(),&vRunFloat_64_BP_after[0],&vEffBPLB_64_after[0],NULL,&vErrEffBPLB_64_after[0]); //MT11
    gEffLB_64_BP_after->SetName("gEffLB_64_BP_before");
    gEffLB_64_BP_after->SetTitle("gEffLB_64_BP_before");
    gEffLB_64_BP_after->SetMarkerSize(1.4);
    gEffLB_64_BP_after->SetMarkerColor(kOrange);
    gEffLB_64_BP_after->SetMarkerStyle(24);

    //BP planes MT12 459, 264, 298
    TGraphErrors *gEffLB_459_BP_before = new TGraphErrors(vRunFloat_459_BP_before.size(),&vRunFloat_459_BP_before[0],&vEffBPLB_459_before[0],NULL,&vErrEffBPLB_459_before[0]); //MT11
    gEffLB_459_BP_before->SetName("gEffLB_459_BP_before");
    gEffLB_459_BP_before->SetTitle("gEffLB_459_BP_before");
    gEffLB_459_BP_before->SetMarkerSize(1.4);
    gEffLB_459_BP_before->SetMarkerColor(kBlack);
    gEffLB_459_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_459_BP_after = new TGraphErrors(vRunFloat_459_BP_after.size(),&vRunFloat_459_BP_after[0],&vEffBPLB_459_after[0],NULL,&vErrEffBPLB_459_after[0]); //MT11
    gEffLB_459_BP_after->SetName("gEffLB_459_BP_after");
    gEffLB_459_BP_after->SetTitle("gEffLB_459_BP_after");
    gEffLB_459_BP_after->SetMarkerSize(1.4);
    gEffLB_459_BP_after->SetMarkerColor(kRed);
    gEffLB_459_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_264_BP_before = new TGraphErrors(vRunFloat_264_BP_before.size(),&vRunFloat_264_BP_before[0],&vEffBPLB_264_before[0],NULL,&vErrEffBPLB_264_before[0]); //MT11
    gEffLB_264_BP_before->SetName("gEffLB_264_BP_before");
    gEffLB_264_BP_before->SetTitle("gEffLB_264_BP_before");
    gEffLB_264_BP_before->SetMarkerSize(1.4);
    gEffLB_264_BP_before->SetMarkerColor(kGreen+3);
    gEffLB_264_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_264_BP_after = new TGraphErrors(vRunFloat_264_BP_after.size(),&vRunFloat_264_BP_after[0],&vEffBPLB_264_after[0],NULL,&vErrEffBPLB_264_after[0]); //MT11
    gEffLB_264_BP_after->SetName("gEffLB_264_BP_before");
    gEffLB_264_BP_after->SetTitle("gEffLB_264_BP_before");
    gEffLB_264_BP_after->SetMarkerSize(1.4);
    gEffLB_264_BP_after->SetMarkerColor(kBlue);
    gEffLB_264_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_298_BP_before = new TGraphErrors(vRunFloat_298_BP_before.size(),&vRunFloat_298_BP_before[0],&vEffBPLB_298_before[0],NULL,&vErrEffBPLB_298_before[0]); //MT11
    gEffLB_298_BP_before->SetName("gEffLB_298_BP_before");
    gEffLB_298_BP_before->SetTitle("gEffLB_298_BP_before");
    gEffLB_298_BP_before->SetMarkerSize(1.4);
    gEffLB_298_BP_before->SetMarkerColor(kMagenta);
    gEffLB_298_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_298_BP_after = new TGraphErrors(vRunFloat_298_BP_after.size(),&vRunFloat_298_BP_after[0],&vEffBPLB_298_after[0],NULL,&vErrEffBPLB_298_after[0]); //MT11
    gEffLB_298_BP_after->SetName("gEffLB_298_BP_before");
    gEffLB_298_BP_after->SetTitle("gEffLB_298_BP_before");
    gEffLB_298_BP_after->SetMarkerSize(1.4);
    gEffLB_298_BP_after->SetMarkerColor(kOrange);
    gEffLB_298_BP_after->SetMarkerStyle(24);

    //BP planes MT21 693, 498, 532
    TGraphErrors *gEffLB_693_BP_before = new TGraphErrors(vRunFloat_693_BP_before.size(),&vRunFloat_693_BP_before[0],&vEffBPLB_693_before[0],NULL,&vErrEffBPLB_693_before[0]); //MT11
    gEffLB_693_BP_before->SetName("gEffLB_693_BP_before");
    gEffLB_693_BP_before->SetTitle("gEffLB_693_BP_before");
    gEffLB_693_BP_before->SetMarkerSize(1.4);
    gEffLB_693_BP_before->SetMarkerColor(kBlack);
    gEffLB_693_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_693_BP_after = new TGraphErrors(vRunFloat_693_BP_after.size(),&vRunFloat_693_BP_after[0],&vEffBPLB_693_after[0],NULL,&vErrEffBPLB_693_after[0]); //MT11
    gEffLB_693_BP_after->SetName("gEffLB_693_BP_after");
    gEffLB_693_BP_after->SetTitle("gEffLB_693_BP_after");
    gEffLB_693_BP_after->SetMarkerSize(1.4);
    gEffLB_693_BP_after->SetMarkerColor(kRed);
    gEffLB_693_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_498_BP_before = new TGraphErrors(vRunFloat_498_BP_before.size(),&vRunFloat_498_BP_before[0],&vEffBPLB_498_before[0],NULL,&vErrEffBPLB_498_before[0]); //MT11
    gEffLB_498_BP_before->SetName("gEffLB_498_BP_before");
    gEffLB_498_BP_before->SetTitle("gEffLB_498_BP_before");
    gEffLB_498_BP_before->SetMarkerSize(1.4);
    gEffLB_498_BP_before->SetMarkerColor(kGreen+3);
    gEffLB_498_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_498_BP_after = new TGraphErrors(vRunFloat_498_BP_after.size(),&vRunFloat_498_BP_after[0],&vEffBPLB_498_after[0],NULL,&vErrEffBPLB_498_after[0]); //MT11
    gEffLB_498_BP_after->SetName("gEffLB_498_BP_before");
    gEffLB_498_BP_after->SetTitle("gEffLB_498_BP_before");
    gEffLB_498_BP_after->SetMarkerSize(1.4);
    gEffLB_498_BP_after->SetMarkerColor(kBlue);
    gEffLB_498_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_532_BP_before = new TGraphErrors(vRunFloat_532_BP_before.size(),&vRunFloat_532_BP_before[0],&vEffBPLB_532_before[0],NULL,&vErrEffBPLB_532_before[0]); //MT11
    gEffLB_532_BP_before->SetName("gEffLB_532_BP_before");
    gEffLB_532_BP_before->SetTitle("gEffLB_532_BP_before");
    gEffLB_532_BP_before->SetMarkerSize(1.4);
    gEffLB_532_BP_before->SetMarkerColor(kMagenta);
    gEffLB_532_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_532_BP_after = new TGraphErrors(vRunFloat_532_BP_after.size(),&vRunFloat_532_BP_after[0],&vEffBPLB_532_after[0],NULL,&vErrEffBPLB_532_after[0]); //MT11
    gEffLB_532_BP_after->SetName("gEffLB_532_BP_before");
    gEffLB_532_BP_after->SetTitle("gEffLB_532_BP_before");
    gEffLB_532_BP_after->SetMarkerSize(1.4);
    gEffLB_532_BP_after->SetMarkerColor(kOrange);
    gEffLB_532_BP_after->SetMarkerStyle(24);

    //BP planes MT22 927, 732, 766
    TGraphErrors *gEffLB_927_BP_before = new TGraphErrors(vRunFloat_927_BP_before.size(),&vRunFloat_927_BP_before[0],&vEffBPLB_927_before[0],NULL,&vErrEffBPLB_927_before[0]); //MT11
    gEffLB_927_BP_before->SetName("gEffLB_927_BP_before");
    gEffLB_927_BP_before->SetTitle("gEffLB_927_BP_before");
    gEffLB_927_BP_before->SetMarkerSize(1.4);
    gEffLB_927_BP_before->SetMarkerColor(kBlack);
    gEffLB_927_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_927_BP_after = new TGraphErrors(vRunFloat_927_BP_after.size(),&vRunFloat_927_BP_after[0],&vEffBPLB_927_after[0],NULL,&vErrEffBPLB_927_after[0]); //MT11
    gEffLB_927_BP_after->SetName("gEffLB_927_BP_after");
    gEffLB_927_BP_after->SetTitle("gEffLB_927_BP_after");
    gEffLB_927_BP_after->SetMarkerSize(1.4);
    gEffLB_927_BP_after->SetMarkerColor(kRed);
    gEffLB_927_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_732_BP_before = new TGraphErrors(vRunFloat_732_BP_before.size(),&vRunFloat_732_BP_before[0],&vEffBPLB_732_before[0],NULL,&vErrEffBPLB_732_before[0]); //MT11
    gEffLB_732_BP_before->SetName("gEffLB_732_BP_before");
    gEffLB_732_BP_before->SetTitle("gEffLB_732_BP_before");
    gEffLB_732_BP_before->SetMarkerSize(1.4);
    gEffLB_732_BP_before->SetMarkerColor(kGreen+3);
    gEffLB_732_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_732_BP_after = new TGraphErrors(vRunFloat_732_BP_after.size(),&vRunFloat_732_BP_after[0],&vEffBPLB_732_after[0],NULL,&vErrEffBPLB_732_after[0]); //MT11
    gEffLB_732_BP_after->SetName("gEffLB_732_BP_before");
    gEffLB_732_BP_after->SetTitle("gEffLB_732_BP_before");
    gEffLB_732_BP_after->SetMarkerSize(1.4);
    gEffLB_732_BP_after->SetMarkerColor(kBlue);
    gEffLB_732_BP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_766_BP_before = new TGraphErrors(vRunFloat_766_BP_before.size(),&vRunFloat_766_BP_before[0],&vEffBPLB_766_before[0],NULL,&vErrEffBPLB_766_before[0]); //MT11
    gEffLB_766_BP_before->SetName("gEffLB_766_BP_before");
    gEffLB_766_BP_before->SetTitle("gEffLB_766_BP_before");
    gEffLB_766_BP_before->SetMarkerSize(1.4);
    gEffLB_766_BP_before->SetMarkerColor(kMagenta);
    gEffLB_766_BP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_766_BP_after = new TGraphErrors(vRunFloat_766_BP_after.size(),&vRunFloat_766_BP_after[0],&vEffBPLB_766_after[0],NULL,&vErrEffBPLB_766_after[0]); //MT11
    gEffLB_766_BP_after->SetName("gEffLB_766_BP_before");
    gEffLB_766_BP_after->SetTitle("gEffLB_766_BP_before");
    gEffLB_766_BP_after->SetMarkerSize(1.4);
    gEffLB_766_BP_after->SetMarkerColor(kOrange);
    gEffLB_766_BP_after->SetMarkerStyle(24);

    //NBP MT11: 225, 30, 64
    TGraphErrors *gEffLB_225_NBP_before = new TGraphErrors(vRunFloat_225_NBP_before.size(),&vRunFloat_225_NBP_before[0],&vEffNBPLB_225_before[0],NULL,&vErrEffNBPLB_225_before[0]); //MT11
    gEffLB_225_NBP_before->SetName("gEffLB_225_NBP_before");
    gEffLB_225_NBP_before->SetTitle("gEffLB_225_NBP_before");
    gEffLB_225_NBP_before->SetMarkerSize(1.4);
    gEffLB_225_NBP_before->SetMarkerColor(kBlack);
    gEffLB_225_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_225_NBP_after = new TGraphErrors(vRunFloat_225_NBP_after.size(),&vRunFloat_225_NBP_after[0],&vEffNBPLB_225_after[0],NULL,&vErrEffNBPLB_225_after[0]); //MT11
    gEffLB_225_NBP_after->SetName("gEffLB_225_NBP_before");
    gEffLB_225_NBP_after->SetTitle("gEffLB_225_NBP_before");
    gEffLB_225_NBP_after->SetMarkerSize(1.4);
    gEffLB_225_NBP_after->SetMarkerColor(kRed);
    gEffLB_225_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_30_NBP_before = new TGraphErrors(vRunFloat_30_NBP_before.size(),&vRunFloat_30_NBP_before[0],&vEffNBPLB_30_before[0],NULL,&vErrEffNBPLB_30_before[0]); //MT11
    gEffLB_30_NBP_before->SetName("gEffLB_30_NBP_before");
    gEffLB_30_NBP_before->SetTitle("gEffLB_30_NBP_before");
    gEffLB_30_NBP_before->SetMarkerSize(1.4);
    gEffLB_30_NBP_before->SetMarkerColor(kGreen+3);
    gEffLB_30_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_30_NBP_after = new TGraphErrors(vRunFloat_30_NBP_after.size(),&vRunFloat_30_NBP_after[0],&vEffNBPLB_30_after[0],NULL,&vErrEffNBPLB_30_after[0]); //MT11
    gEffLB_30_NBP_after->SetName("gEffLB_30_NBP_before");
    gEffLB_30_NBP_after->SetTitle("gEffLB_30_NBP_before");
    gEffLB_30_NBP_after->SetMarkerSize(1.4);
    gEffLB_30_NBP_after->SetMarkerColor(kBlue);
    gEffLB_30_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_64_NBP_before = new TGraphErrors(vRunFloat_64_NBP_before.size(),&vRunFloat_64_NBP_before[0],&vEffNBPLB_64_before[0],NULL,&vErrEffNBPLB_64_before[0]); //MT11
    gEffLB_64_NBP_before->SetName("gEffLB_64_NBP_before");
    gEffLB_64_NBP_before->SetTitle("gEffLB_64_NBP_before");
    gEffLB_64_NBP_before->SetMarkerSize(1.4);
    gEffLB_64_NBP_before->SetMarkerColor(kMagenta);
    gEffLB_64_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_64_NBP_after = new TGraphErrors(vRunFloat_64_NBP_after.size(),&vRunFloat_64_NBP_after[0],&vEffNBPLB_64_after[0],NULL,&vErrEffNBPLB_64_after[0]); //MT11
    gEffLB_64_NBP_after->SetName("gEffLB_64_NBP_before");
    gEffLB_64_NBP_after->SetTitle("gEffLB_64_NBP_before");
    gEffLB_64_NBP_after->SetMarkerSize(1.4);
    gEffLB_64_NBP_after->SetMarkerColor(kOrange);
    gEffLB_64_NBP_after->SetMarkerStyle(24);

    //NBP planes MT12 459, 264, 298
    TGraphErrors *gEffLB_459_NBP_before = new TGraphErrors(vRunFloat_459_NBP_before.size(),&vRunFloat_459_NBP_before[0],&vEffNBPLB_459_before[0],NULL,&vErrEffNBPLB_459_before[0]); //MT11
    gEffLB_459_NBP_before->SetName("gEffLB_459_NBP_before");
    gEffLB_459_NBP_before->SetTitle("gEffLB_459_NBP_before");
    gEffLB_459_NBP_before->SetMarkerSize(1.4);
    gEffLB_459_NBP_before->SetMarkerColor(kBlack);
    gEffLB_459_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_459_NBP_after = new TGraphErrors(vRunFloat_459_NBP_after.size(),&vRunFloat_459_NBP_after[0],&vEffNBPLB_459_after[0],NULL,&vErrEffNBPLB_459_after[0]); //MT11
    gEffLB_459_NBP_after->SetName("gEffLB_459_NBP_after");
    gEffLB_459_NBP_after->SetTitle("gEffLB_459_NBP_after");
    gEffLB_459_NBP_after->SetMarkerSize(1.4);
    gEffLB_459_NBP_after->SetMarkerColor(kRed);
    gEffLB_459_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_264_NBP_before = new TGraphErrors(vRunFloat_264_NBP_before.size(),&vRunFloat_264_NBP_before[0],&vEffNBPLB_264_before[0],NULL,&vErrEffNBPLB_264_before[0]); //MT11
    gEffLB_264_NBP_before->SetName("gEffLB_264_NBP_before");
    gEffLB_264_NBP_before->SetTitle("gEffLB_264_NBP_before");
    gEffLB_264_NBP_before->SetMarkerSize(1.4);
    gEffLB_264_NBP_before->SetMarkerColor(kGreen+3);
    gEffLB_264_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_264_NBP_after = new TGraphErrors(vRunFloat_264_NBP_after.size(),&vRunFloat_264_NBP_after[0],&vEffNBPLB_264_after[0],NULL,&vErrEffNBPLB_264_after[0]); //MT11
    gEffLB_264_NBP_after->SetName("gEffLB_264_NBP_before");
    gEffLB_264_NBP_after->SetTitle("gEffLB_264_NBP_before");
    gEffLB_264_NBP_after->SetMarkerSize(1.4);
    gEffLB_264_NBP_after->SetMarkerColor(kBlue);
    gEffLB_264_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_298_NBP_before = new TGraphErrors(vRunFloat_298_NBP_before.size(),&vRunFloat_298_NBP_before[0],&vEffNBPLB_298_before[0],NULL,&vErrEffNBPLB_298_before[0]); //MT11
    gEffLB_298_NBP_before->SetName("gEffLB_298_NBP_before");
    gEffLB_298_NBP_before->SetTitle("gEffLB_298_NBP_before");
    gEffLB_298_NBP_before->SetMarkerSize(1.4);
    gEffLB_298_NBP_before->SetMarkerColor(kMagenta);
    gEffLB_298_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_298_NBP_after = new TGraphErrors(vRunFloat_298_NBP_after.size(),&vRunFloat_298_NBP_after[0],&vEffNBPLB_298_after[0],NULL,&vErrEffNBPLB_298_after[0]); //MT11
    gEffLB_298_NBP_after->SetName("gEffLB_298_NBP_before");
    gEffLB_298_NBP_after->SetTitle("gEffLB_298_NBP_before");
    gEffLB_298_NBP_after->SetMarkerSize(1.4);
    gEffLB_298_NBP_after->SetMarkerColor(kOrange);
    gEffLB_298_NBP_after->SetMarkerStyle(24);

    //NBP planes MT21 693, 498, 532
    TGraphErrors *gEffLB_693_NBP_before = new TGraphErrors(vRunFloat_693_NBP_before.size(),&vRunFloat_693_NBP_before[0],&vEffNBPLB_693_before[0],NULL,&vErrEffNBPLB_693_before[0]); //MT11
    gEffLB_693_NBP_before->SetName("gEffLB_693_NBP_before");
    gEffLB_693_NBP_before->SetTitle("gEffLB_693_NBP_before");
    gEffLB_693_NBP_before->SetMarkerSize(1.4);
    gEffLB_693_NBP_before->SetMarkerColor(kBlack);
    gEffLB_693_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_693_NBP_after = new TGraphErrors(vRunFloat_693_NBP_after.size(),&vRunFloat_693_NBP_after[0],&vEffNBPLB_693_after[0],NULL,&vErrEffNBPLB_693_after[0]); //MT11
    gEffLB_693_NBP_after->SetName("gEffLB_693_NBP_after");
    gEffLB_693_NBP_after->SetTitle("gEffLB_693_NBP_after");
    gEffLB_693_NBP_after->SetMarkerSize(1.4);
    gEffLB_693_NBP_after->SetMarkerColor(kRed);
    gEffLB_693_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_498_NBP_before = new TGraphErrors(vRunFloat_498_NBP_before.size(),&vRunFloat_498_NBP_before[0],&vEffNBPLB_498_before[0],NULL,&vErrEffNBPLB_498_before[0]); //MT11
    gEffLB_498_NBP_before->SetName("gEffLB_498_NBP_before");
    gEffLB_498_NBP_before->SetTitle("gEffLB_498_NBP_before");
    gEffLB_498_NBP_before->SetMarkerSize(1.4);
    gEffLB_498_NBP_before->SetMarkerColor(kGreen+3);
    gEffLB_498_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_498_NBP_after = new TGraphErrors(vRunFloat_498_NBP_after.size(),&vRunFloat_498_NBP_after[0],&vEffNBPLB_498_after[0],NULL,&vErrEffNBPLB_498_after[0]); //MT11
    gEffLB_498_NBP_after->SetName("gEffLB_498_NBP_before");
    gEffLB_498_NBP_after->SetTitle("gEffLB_498_NBP_before");
    gEffLB_498_NBP_after->SetMarkerSize(1.4);
    gEffLB_498_NBP_after->SetMarkerColor(kBlue);
    gEffLB_498_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_532_NBP_before = new TGraphErrors(vRunFloat_532_NBP_before.size(),&vRunFloat_532_NBP_before[0],&vEffNBPLB_532_before[0],NULL,&vErrEffNBPLB_532_before[0]); //MT11
    gEffLB_532_NBP_before->SetName("gEffLB_532_NBP_before");
    gEffLB_532_NBP_before->SetTitle("gEffLB_532_NBP_before");
    gEffLB_532_NBP_before->SetMarkerSize(1.4);
    gEffLB_532_NBP_before->SetMarkerColor(kMagenta);
    gEffLB_532_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_532_NBP_after = new TGraphErrors(vRunFloat_532_NBP_after.size(),&vRunFloat_532_NBP_after[0],&vEffNBPLB_532_after[0],NULL,&vErrEffNBPLB_532_after[0]); //MT11
    gEffLB_532_NBP_after->SetName("gEffLB_532_NBP_before");
    gEffLB_532_NBP_after->SetTitle("gEffLB_532_NBP_before");
    gEffLB_532_NBP_after->SetMarkerSize(1.4);
    gEffLB_532_NBP_after->SetMarkerColor(kOrange);
    gEffLB_532_NBP_after->SetMarkerStyle(24);

    //NBP planes MT22 927, 732, 766
    TGraphErrors *gEffLB_927_NBP_before = new TGraphErrors(vRunFloat_927_NBP_before.size(),&vRunFloat_927_NBP_before[0],&vEffNBPLB_927_before[0],NULL,&vErrEffNBPLB_927_before[0]); //MT11
    gEffLB_927_NBP_before->SetName("gEffLB_927_NBP_before");
    gEffLB_927_NBP_before->SetTitle("gEffLB_927_NBP_before");
    gEffLB_927_NBP_before->SetMarkerSize(1.4);
    gEffLB_927_NBP_before->SetMarkerColor(kBlack);
    gEffLB_927_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_927_NBP_after = new TGraphErrors(vRunFloat_927_NBP_after.size(),&vRunFloat_927_NBP_after[0],&vEffNBPLB_927_after[0],NULL,&vErrEffNBPLB_927_after[0]); //MT11
    gEffLB_927_NBP_after->SetName("gEffLB_927_NBP_after");
    gEffLB_927_NBP_after->SetTitle("gEffLB_927_NBP_after");
    gEffLB_927_NBP_after->SetMarkerSize(1.4);
    gEffLB_927_NBP_after->SetMarkerColor(kRed);
    gEffLB_927_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_732_NBP_before = new TGraphErrors(vRunFloat_732_NBP_before.size(),&vRunFloat_732_NBP_before[0],&vEffNBPLB_732_before[0],NULL,&vErrEffNBPLB_732_before[0]); //MT11
    gEffLB_732_NBP_before->SetName("gEffLB_732_NBP_before");
    gEffLB_732_NBP_before->SetTitle("gEffLB_732_NBP_before");
    gEffLB_732_NBP_before->SetMarkerSize(1.4);
    gEffLB_732_NBP_before->SetMarkerColor(kGreen+3);
    gEffLB_732_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_732_NBP_after = new TGraphErrors(vRunFloat_732_NBP_after.size(),&vRunFloat_732_NBP_after[0],&vEffNBPLB_732_after[0],NULL,&vErrEffNBPLB_732_after[0]); //MT11
    gEffLB_732_NBP_after->SetName("gEffLB_732_NBP_before");
    gEffLB_732_NBP_after->SetTitle("gEffLB_732_NBP_before");
    gEffLB_732_NBP_after->SetMarkerSize(1.4);
    gEffLB_732_NBP_after->SetMarkerColor(kBlue);
    gEffLB_732_NBP_after->SetMarkerStyle(24);
    //
    TGraphErrors *gEffLB_766_NBP_before = new TGraphErrors(vRunFloat_766_NBP_before.size(),&vRunFloat_766_NBP_before[0],&vEffNBPLB_766_before[0],NULL,&vErrEffNBPLB_766_before[0]); //MT11
    gEffLB_766_NBP_before->SetName("gEffLB_766_NBP_before");
    gEffLB_766_NBP_before->SetTitle("gEffLB_766_NBP_before");
    gEffLB_766_NBP_before->SetMarkerSize(1.4);
    gEffLB_766_NBP_before->SetMarkerColor(kMagenta);
    gEffLB_766_NBP_before->SetMarkerStyle(8);
    TGraphErrors *gEffLB_766_NBP_after = new TGraphErrors(vRunFloat_766_NBP_after.size(),&vRunFloat_766_NBP_after[0],&vEffNBPLB_766_after[0],NULL,&vErrEffNBPLB_766_after[0]); //MT11
    gEffLB_766_NBP_after->SetName("gEffLB_766_NBP_before");
    gEffLB_766_NBP_after->SetTitle("gEffLB_766_NBP_before");
    gEffLB_766_NBP_after->SetMarkerSize(1.4);
    gEffLB_766_NBP_after->SetMarkerColor(kOrange);
    gEffLB_766_NBP_after->SetMarkerStyle(24);

    //-------------------//
    //
    // TOT COUNTS PER LB //
    //
    //-------------------//

    //Both planes MT11: 225, 30, 64
    //vTotCountsLB_766_after_BP
    TGraphErrors *gTotCounts_225_Both_before = new TGraphErrors(vRunFloat_225_Both_before.size(),&vRunFloat_225_Both_before[0],&vTotCountsLB_225_before_both[0],NULL,NULL); //MT11
    gTotCounts_225_Both_before->SetName("gTotCounts_225_Both_before");
    gTotCounts_225_Both_before->SetTitle("gTotCounts_225_Both_before");
    gTotCounts_225_Both_before->SetMarkerSize(1.4);
    gTotCounts_225_Both_before->SetMarkerColor(kBlack);
    gTotCounts_225_Both_before->SetMarkerStyle(22);
    lLB225->AddEntry(gTotCounts_225_Both_before,"Before p_{t} cut","p");
    TGraphErrors *gTotCounts_225_Both_after = new TGraphErrors(vRunFloat_225_Both_after.size(),&vRunFloat_225_Both_after[0],&vTotCountsLB_225_after_both[0],NULL,NULL); //MT11
    gTotCounts_225_Both_after->SetName("gTotCounts_225_Both_before");
    gTotCounts_225_Both_after->SetTitle("gTotCounts_225_Both_before");
    gTotCounts_225_Both_after->SetMarkerSize(1.4);
    gTotCounts_225_Both_after->SetMarkerColor(kRed);
    gTotCounts_225_Both_after->SetMarkerStyle(23);
    lLB225->AddEntry(gTotCounts_225_Both_after,"After p_{t} cut","p");
    //
    TGraphErrors *gTotCounts_30_Both_before = new TGraphErrors(vRunFloat_30_Both_before.size(),&vRunFloat_30_Both_before[0],&vTotCountsLB_30_before_both[0],NULL,NULL); //MT11
    gTotCounts_30_Both_before->SetName("gTotCounts_30_Both_before");
    gTotCounts_30_Both_before->SetTitle("gTotCounts_30_Both_before");
    gTotCounts_30_Both_before->SetMarkerSize(1.4);
    gTotCounts_30_Both_before->SetMarkerColor(kGreen+3);
    gTotCounts_30_Both_before->SetMarkerStyle(22);
    lLB30->AddEntry(gTotCounts_30_Both_before,"Before p_{t} cut","p");
    TGraphErrors *gTotCounts_30_Both_after = new TGraphErrors(vRunFloat_30_Both_after.size(),&vRunFloat_30_Both_after[0],&vTotCountsLB_30_after_both[0],NULL,NULL); //MT11
    gTotCounts_30_Both_after->SetName("gTotCounts_30_Both_before");
    gTotCounts_30_Both_after->SetTitle("gTotCounts_30_Both_before");
    gTotCounts_30_Both_after->SetMarkerSize(1.4);
    gTotCounts_30_Both_after->SetMarkerColor(kBlue);
    gTotCounts_30_Both_after->SetMarkerStyle(23);
    lLB30->AddEntry(gTotCounts_30_Both_after,"After p_{t} cut","p");
    //
    TGraphErrors *gTotCounts_64_Both_before = new TGraphErrors(vRunFloat_64_Both_before.size(),&vRunFloat_64_Both_before[0],&vTotCountsLB_64_before_both[0],NULL,NULL); //MT11
    gTotCounts_64_Both_before->SetName("gTotCounts_64_Both_before");
    gTotCounts_64_Both_before->SetTitle("gTotCounts_64_Both_before");
    gTotCounts_64_Both_before->SetMarkerSize(1.4);
    gTotCounts_64_Both_before->SetMarkerColor(kMagenta);
    gTotCounts_64_Both_before->SetMarkerStyle(22);
    lLB64->AddEntry(gTotCounts_64_Both_before,"Before p_{t} cut","p");
    TGraphErrors *gTotCounts_64_Both_after = new TGraphErrors(vRunFloat_64_Both_after.size(),&vRunFloat_64_Both_after[0],&vTotCountsLB_64_after_both[0],NULL,NULL); //MT11
    gTotCounts_64_Both_after->SetName("gTotCounts_64_Both_before");
    gTotCounts_64_Both_after->SetTitle("gTotCounts_64_Both_before");
    gTotCounts_64_Both_after->SetMarkerSize(1.4);
    gTotCounts_64_Both_after->SetMarkerColor(kOrange);
    gTotCounts_64_Both_after->SetMarkerStyle(23);
    lLB64->AddEntry(gTotCounts_64_Both_after,"After p_{t} cut","p");

    //Both planes MT12 459, 264, 298
    TGraphErrors *gTotCounts_459_Both_before = new TGraphErrors(vRunFloat_459_Both_before.size(),&vRunFloat_459_Both_before[0],&vTotCountsLB_459_before_both[0],NULL,NULL); //MT11
    gTotCounts_459_Both_before->SetName("gTotCounts_459_Both_before");
    gTotCounts_459_Both_before->SetTitle("gTotCounts_459_Both_before");
    gTotCounts_459_Both_before->SetMarkerSize(1.4);
    gTotCounts_459_Both_before->SetMarkerColor(kBlack);
    gTotCounts_459_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_459_Both_after = new TGraphErrors(vRunFloat_459_Both_after.size(),&vRunFloat_459_Both_after[0],&vTotCountsLB_459_after_both[0],NULL,NULL); //MT11
    gTotCounts_459_Both_after->SetName("gTotCounts_459_Both_after");
    gTotCounts_459_Both_after->SetTitle("gTotCounts_459_Both_after");
    gTotCounts_459_Both_after->SetMarkerSize(1.4);
    gTotCounts_459_Both_after->SetMarkerColor(kRed);
    gTotCounts_459_Both_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_264_Both_before = new TGraphErrors(vRunFloat_264_Both_before.size(),&vRunFloat_264_Both_before[0],&vTotCountsLB_264_before_both[0],NULL,NULL); //MT11
    gTotCounts_264_Both_before->SetName("gTotCounts_264_Both_before");
    gTotCounts_264_Both_before->SetTitle("gTotCounts_264_Both_before");
    gTotCounts_264_Both_before->SetMarkerSize(1.4);
    gTotCounts_264_Both_before->SetMarkerColor(kGreen+3);
    gTotCounts_264_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_264_Both_after = new TGraphErrors(vRunFloat_264_Both_after.size(),&vRunFloat_264_Both_after[0],&vTotCountsLB_264_after_both[0],NULL,NULL); //MT11
    gTotCounts_264_Both_after->SetName("gTotCounts_264_Both_before");
    gTotCounts_264_Both_after->SetTitle("gTotCounts_264_Both_before");
    gTotCounts_264_Both_after->SetMarkerSize(1.4);
    gTotCounts_264_Both_after->SetMarkerColor(kBlue);
    gTotCounts_264_Both_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_298_Both_before = new TGraphErrors(vRunFloat_298_Both_before.size(),&vRunFloat_298_Both_before[0],&vTotCountsLB_298_before_both[0],NULL,NULL); //MT11
    gTotCounts_298_Both_before->SetName("gTotCounts_298_Both_before");
    gTotCounts_298_Both_before->SetTitle("gTotCounts_298_Both_before");
    gTotCounts_298_Both_before->SetMarkerSize(1.4);
    gTotCounts_298_Both_before->SetMarkerColor(kMagenta);
    gTotCounts_298_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_298_Both_after = new TGraphErrors(vRunFloat_298_Both_after.size(),&vRunFloat_298_Both_after[0],&vTotCountsLB_298_after_both[0],NULL,NULL); //MT11
    gTotCounts_298_Both_after->SetName("gTotCounts_298_Both_before");
    gTotCounts_298_Both_after->SetTitle("gTotCounts_298_Both_before");
    gTotCounts_298_Both_after->SetMarkerSize(1.4);
    gTotCounts_298_Both_after->SetMarkerColor(kOrange);
    gTotCounts_298_Both_after->SetMarkerStyle(23);

    //Both planes MT21 693, 498, 532
    TGraphErrors *gTotCounts_693_Both_before = new TGraphErrors(vRunFloat_693_Both_before.size(),&vRunFloat_693_Both_before[0],&vTotCountsLB_693_before_both[0],NULL,NULL); //MT11
    gTotCounts_693_Both_before->SetName("gTotCounts_693_Both_before");
    gTotCounts_693_Both_before->SetTitle("gTotCounts_693_Both_before");
    gTotCounts_693_Both_before->SetMarkerSize(1.4);
    gTotCounts_693_Both_before->SetMarkerColor(kBlack);
    gTotCounts_693_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_693_Both_after = new TGraphErrors(vRunFloat_693_Both_after.size(),&vRunFloat_693_Both_after[0],&vTotCountsLB_693_after_both[0],NULL,NULL); //MT11
    gTotCounts_693_Both_after->SetName("gTotCounts_693_Both_after");
    gTotCounts_693_Both_after->SetTitle("gTotCounts_693_Both_after");
    gTotCounts_693_Both_after->SetMarkerSize(1.4);
    gTotCounts_693_Both_after->SetMarkerColor(kRed);
    gTotCounts_693_Both_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_498_Both_before = new TGraphErrors(vRunFloat_498_Both_before.size(),&vRunFloat_498_Both_before[0],&vTotCountsLB_498_before_both[0],NULL,NULL); //MT11
    gTotCounts_498_Both_before->SetName("gTotCounts_498_Both_before");
    gTotCounts_498_Both_before->SetTitle("gTotCounts_498_Both_before");
    gTotCounts_498_Both_before->SetMarkerSize(1.4);
    gTotCounts_498_Both_before->SetMarkerColor(kGreen+3);
    gTotCounts_498_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_498_Both_after = new TGraphErrors(vRunFloat_498_Both_after.size(),&vRunFloat_498_Both_after[0],&vTotCountsLB_498_after_both[0],NULL,NULL); //MT11
    gTotCounts_498_Both_after->SetName("gTotCounts_498_Both_before");
    gTotCounts_498_Both_after->SetTitle("gTotCounts_498_Both_before");
    gTotCounts_498_Both_after->SetMarkerSize(1.4);
    gTotCounts_498_Both_after->SetMarkerColor(kBlue);
    gTotCounts_498_Both_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_532_Both_before = new TGraphErrors(vRunFloat_532_Both_before.size(),&vRunFloat_532_Both_before[0],&vTotCountsLB_532_before_both[0],NULL,NULL); //MT11
    gTotCounts_532_Both_before->SetName("gTotCounts_532_Both_before");
    gTotCounts_532_Both_before->SetTitle("gTotCounts_532_Both_before");
    gTotCounts_532_Both_before->SetMarkerSize(1.4);
    gTotCounts_532_Both_before->SetMarkerColor(kMagenta);
    gTotCounts_532_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_532_Both_after = new TGraphErrors(vRunFloat_532_Both_after.size(),&vRunFloat_532_Both_after[0],&vTotCountsLB_532_after_both[0],NULL,NULL); //MT11
    gTotCounts_532_Both_after->SetName("gTotCounts_532_Both_before");
    gTotCounts_532_Both_after->SetTitle("gTotCounts_532_Both_before");
    gTotCounts_532_Both_after->SetMarkerSize(1.4);
    gTotCounts_532_Both_after->SetMarkerColor(kOrange);
    gTotCounts_532_Both_after->SetMarkerStyle(23);

    //Both planes MT22 927, 732, 766
    TGraphErrors *gTotCounts_927_Both_before = new TGraphErrors(vRunFloat_927_Both_before.size(),&vRunFloat_927_Both_before[0],&vTotCountsLB_927_before_both[0],NULL,NULL); //MT11
    gTotCounts_927_Both_before->SetName("gTotCounts_927_Both_before");
    gTotCounts_927_Both_before->SetTitle("gTotCounts_927_Both_before");
    gTotCounts_927_Both_before->SetMarkerSize(1.4);
    gTotCounts_927_Both_before->SetMarkerColor(kBlack);
    gTotCounts_927_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_927_Both_after = new TGraphErrors(vRunFloat_927_Both_after.size(),&vRunFloat_927_Both_after[0],&vTotCountsLB_927_after_both[0],NULL,NULL); //MT11
    gTotCounts_927_Both_after->SetName("gTotCounts_927_Both_after");
    gTotCounts_927_Both_after->SetTitle("gTotCounts_927_Both_after");
    gTotCounts_927_Both_after->SetMarkerSize(1.4);
    gTotCounts_927_Both_after->SetMarkerColor(kRed);
    gTotCounts_927_Both_after->SetMarkerStyle(23);

    //
    TGraphErrors *gTotCounts_732_Both_before = new TGraphErrors(vRunFloat_732_Both_before.size(),&vRunFloat_732_Both_before[0],&vTotCountsLB_732_before_both[0],NULL,NULL); //MT11
    gTotCounts_732_Both_before->SetName("gTotCounts_732_Both_before");
    gTotCounts_732_Both_before->SetTitle("gTotCounts_732_Both_before");
    gTotCounts_732_Both_before->SetMarkerSize(1.4);
    gTotCounts_732_Both_before->SetMarkerColor(kGreen+3);
    gTotCounts_732_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_732_Both_after = new TGraphErrors(vRunFloat_732_Both_after.size(),&vRunFloat_732_Both_after[0],&vTotCountsLB_732_after_both[0],NULL,NULL); //MT11
    gTotCounts_732_Both_after->SetName("gTotCounts_732_Both_before");
    gTotCounts_732_Both_after->SetTitle("gTotCounts_732_Both_before");
    gTotCounts_732_Both_after->SetMarkerSize(1.4);
    gTotCounts_732_Both_after->SetMarkerColor(kBlue);
    gTotCounts_732_Both_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_766_Both_before = new TGraphErrors(vRunFloat_766_Both_before.size(),&vRunFloat_766_Both_before[0],&vTotCountsLB_766_before_both[0],NULL,NULL); //MT11
    gTotCounts_766_Both_before->SetName("gTotCounts_766_Both_before");
    gTotCounts_766_Both_before->SetTitle("gTotCounts_766_Both_before");
    gTotCounts_766_Both_before->SetMarkerSize(1.4);
    gTotCounts_766_Both_before->SetMarkerColor(kMagenta);
    gTotCounts_766_Both_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_766_Both_after = new TGraphErrors(vRunFloat_766_Both_after.size(),&vRunFloat_766_Both_after[0],&vTotCountsLB_766_after_both[0],NULL,NULL); //MT11
    gTotCounts_766_Both_after->SetName("gTotCounts_766_Both_before");
    gTotCounts_766_Both_after->SetTitle("gTotCounts_766_Both_before");
    gTotCounts_766_Both_after->SetMarkerSize(1.4);
    gTotCounts_766_Both_after->SetMarkerColor(kOrange);
    gTotCounts_766_Both_after->SetMarkerStyle(23);

    //BP MT11: 225, 30, 64
    TGraphErrors *gTotCounts_225_BP_before = new TGraphErrors(vRunFloat_225_BP_before.size(),&vRunFloat_225_BP_before[0],&vTotCountsLB_225_before_BP[0],NULL,NULL); //MT11
    gTotCounts_225_BP_before->SetName("gTotCounts_225_BP_before");
    gTotCounts_225_BP_before->SetTitle("gTotCounts_225_BP_before");
    gTotCounts_225_BP_before->SetMarkerSize(1.4);
    gTotCounts_225_BP_before->SetMarkerColor(kBlack);
    gTotCounts_225_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_225_BP_after = new TGraphErrors(vRunFloat_225_BP_after.size(),&vRunFloat_225_BP_after[0],&vTotCountsLB_225_after_BP[0],NULL,NULL); //MT11
    gTotCounts_225_BP_after->SetName("gTotCounts_225_BP_before");
    gTotCounts_225_BP_after->SetTitle("gTotCounts_225_BP_before");
    gTotCounts_225_BP_after->SetMarkerSize(1.4);
    gTotCounts_225_BP_after->SetMarkerColor(kRed);
    gTotCounts_225_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_30_BP_before = new TGraphErrors(vRunFloat_30_BP_before.size(),&vRunFloat_30_BP_before[0],&vTotCountsLB_30_before_BP[0],NULL,NULL); //MT11
    gTotCounts_30_BP_before->SetName("gTotCounts_30_BP_before");
    gTotCounts_30_BP_before->SetTitle("gTotCounts_30_BP_before");
    gTotCounts_30_BP_before->SetMarkerSize(1.4);
    gTotCounts_30_BP_before->SetMarkerColor(kGreen+3);
    gTotCounts_30_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_30_BP_after = new TGraphErrors(vRunFloat_30_BP_after.size(),&vRunFloat_30_BP_after[0],&vTotCountsLB_30_after_BP[0],NULL,NULL); //MT11
    gTotCounts_30_BP_after->SetName("gTotCounts_30_BP_before");
    gTotCounts_30_BP_after->SetTitle("gTotCounts_30_BP_before");
    gTotCounts_30_BP_after->SetMarkerSize(1.4);
    gTotCounts_30_BP_after->SetMarkerColor(kBlue);
    gTotCounts_30_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_64_BP_before = new TGraphErrors(vRunFloat_64_BP_before.size(),&vRunFloat_64_BP_before[0],&vTotCountsLB_64_before_BP[0],NULL,NULL); //MT11
    gTotCounts_64_BP_before->SetName("gTotCounts_64_BP_before");
    gTotCounts_64_BP_before->SetTitle("gTotCounts_64_BP_before");
    gTotCounts_64_BP_before->SetMarkerSize(1.4);
    gTotCounts_64_BP_before->SetMarkerColor(kMagenta);
    gTotCounts_64_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_64_BP_after = new TGraphErrors(vRunFloat_64_BP_after.size(),&vRunFloat_64_BP_after[0],&vTotCountsLB_64_after_BP[0],NULL,NULL); //MT11
    gTotCounts_64_BP_after->SetName("gTotCounts_64_BP_before");
    gTotCounts_64_BP_after->SetTitle("gTotCounts_64_BP_before");
    gTotCounts_64_BP_after->SetMarkerSize(1.4);
    gTotCounts_64_BP_after->SetMarkerColor(kOrange);
    gTotCounts_64_BP_after->SetMarkerStyle(23);

    //BP planes MT12 459, 264, 298
    TGraphErrors *gTotCounts_459_BP_before = new TGraphErrors(vRunFloat_459_BP_before.size(),&vRunFloat_459_BP_before[0],&vTotCountsLB_459_before_BP[0],NULL,NULL); //MT11
    gTotCounts_459_BP_before->SetName("gTotCounts_459_BP_before");
    gTotCounts_459_BP_before->SetTitle("gTotCounts_459_BP_before");
    gTotCounts_459_BP_before->SetMarkerSize(1.4);
    gTotCounts_459_BP_before->SetMarkerColor(kBlack);
    gTotCounts_459_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_459_BP_after = new TGraphErrors(vRunFloat_459_BP_after.size(),&vRunFloat_459_BP_after[0],&vTotCountsLB_459_after_BP[0],NULL,NULL); //MT11
    gTotCounts_459_BP_after->SetName("gTotCounts_459_BP_after");
    gTotCounts_459_BP_after->SetTitle("gTotCounts_459_BP_after");
    gTotCounts_459_BP_after->SetMarkerSize(1.4);
    gTotCounts_459_BP_after->SetMarkerColor(kRed);
    gTotCounts_459_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_264_BP_before = new TGraphErrors(vRunFloat_264_BP_before.size(),&vRunFloat_264_BP_before[0],&vTotCountsLB_264_before_BP[0],NULL,NULL); //MT11
    gTotCounts_264_BP_before->SetName("gTotCounts_264_BP_before");
    gTotCounts_264_BP_before->SetTitle("gTotCounts_264_BP_before");
    gTotCounts_264_BP_before->SetMarkerSize(1.4);
    gTotCounts_264_BP_before->SetMarkerColor(kGreen+3);
    gTotCounts_264_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_264_BP_after = new TGraphErrors(vRunFloat_264_BP_after.size(),&vRunFloat_264_BP_after[0],&vTotCountsLB_264_after_BP[0],NULL,NULL); //MT11
    gTotCounts_264_BP_after->SetName("gTotCounts_264_BP_before");
    gTotCounts_264_BP_after->SetTitle("gTotCounts_264_BP_before");
    gTotCounts_264_BP_after->SetMarkerSize(1.4);
    gTotCounts_264_BP_after->SetMarkerColor(kBlue);
    gTotCounts_264_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_298_BP_before = new TGraphErrors(vRunFloat_298_BP_before.size(),&vRunFloat_298_BP_before[0],&vTotCountsLB_298_before_BP[0],NULL,NULL); //MT11
    gTotCounts_298_BP_before->SetName("gTotCounts_298_BP_before");
    gTotCounts_298_BP_before->SetTitle("gTotCounts_298_BP_before");
    gTotCounts_298_BP_before->SetMarkerSize(1.4);
    gTotCounts_298_BP_before->SetMarkerColor(kMagenta);
    gTotCounts_298_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_298_BP_after = new TGraphErrors(vRunFloat_298_BP_after.size(),&vRunFloat_298_BP_after[0],&vTotCountsLB_298_after_BP[0],NULL,NULL); //MT11
    gTotCounts_298_BP_after->SetName("gTotCounts_298_BP_before");
    gTotCounts_298_BP_after->SetTitle("gTotCounts_298_BP_before");
    gTotCounts_298_BP_after->SetMarkerSize(1.4);
    gTotCounts_298_BP_after->SetMarkerColor(kOrange);
    gTotCounts_298_BP_after->SetMarkerStyle(23);

    //BP planes MT21 693, 498, 532
    TGraphErrors *gTotCounts_693_BP_before = new TGraphErrors(vRunFloat_693_BP_before.size(),&vRunFloat_693_BP_before[0],&vTotCountsLB_693_before_BP[0],NULL,NULL); //MT11
    gTotCounts_693_BP_before->SetName("gTotCounts_693_BP_before");
    gTotCounts_693_BP_before->SetTitle("gTotCounts_693_BP_before");
    gTotCounts_693_BP_before->SetMarkerSize(1.4);
    gTotCounts_693_BP_before->SetMarkerColor(kBlack);
    gTotCounts_693_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_693_BP_after = new TGraphErrors(vRunFloat_693_BP_after.size(),&vRunFloat_693_BP_after[0],&vTotCountsLB_693_after_BP[0],NULL,NULL); //MT11
    gTotCounts_693_BP_after->SetName("gTotCounts_693_BP_after");
    gTotCounts_693_BP_after->SetTitle("gTotCounts_693_BP_after");
    gTotCounts_693_BP_after->SetMarkerSize(1.4);
    gTotCounts_693_BP_after->SetMarkerColor(kRed);
    gTotCounts_693_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_498_BP_before = new TGraphErrors(vRunFloat_498_BP_before.size(),&vRunFloat_498_BP_before[0],&vTotCountsLB_498_before_BP[0],NULL,NULL); //MT11
    gTotCounts_498_BP_before->SetName("gTotCounts_498_BP_before");
    gTotCounts_498_BP_before->SetTitle("gTotCounts_498_BP_before");
    gTotCounts_498_BP_before->SetMarkerSize(1.4);
    gTotCounts_498_BP_before->SetMarkerColor(kGreen+3);
    gTotCounts_498_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_498_BP_after = new TGraphErrors(vRunFloat_498_BP_after.size(),&vRunFloat_498_BP_after[0],&vTotCountsLB_498_after_BP[0],NULL,NULL); //MT11
    gTotCounts_498_BP_after->SetName("gTotCounts_498_BP_before");
    gTotCounts_498_BP_after->SetTitle("gTotCounts_498_BP_before");
    gTotCounts_498_BP_after->SetMarkerSize(1.4);
    gTotCounts_498_BP_after->SetMarkerColor(kBlue);
    gTotCounts_498_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_532_BP_before = new TGraphErrors(vRunFloat_532_BP_before.size(),&vRunFloat_532_BP_before[0],&vTotCountsLB_532_before_BP[0],NULL,NULL); //MT11
    gTotCounts_532_BP_before->SetName("gTotCounts_532_BP_before");
    gTotCounts_532_BP_before->SetTitle("gTotCounts_532_BP_before");
    gTotCounts_532_BP_before->SetMarkerSize(1.4);
    gTotCounts_532_BP_before->SetMarkerColor(kMagenta);
    gTotCounts_532_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_532_BP_after = new TGraphErrors(vRunFloat_532_BP_after.size(),&vRunFloat_532_BP_after[0],&vTotCountsLB_532_after_BP[0],NULL,NULL); //MT11
    gTotCounts_532_BP_after->SetName("gTotCounts_532_BP_before");
    gTotCounts_532_BP_after->SetTitle("gTotCounts_532_BP_before");
    gTotCounts_532_BP_after->SetMarkerSize(1.4);
    gTotCounts_532_BP_after->SetMarkerColor(kOrange);
    gTotCounts_532_BP_after->SetMarkerStyle(23);

    //BP planes MT22 927, 732, 766
    TGraphErrors *gTotCounts_927_BP_before = new TGraphErrors(vRunFloat_927_BP_before.size(),&vRunFloat_927_BP_before[0],&vTotCountsLB_927_before_BP[0],NULL,NULL); //MT11
    gTotCounts_927_BP_before->SetName("gTotCounts_927_BP_before");
    gTotCounts_927_BP_before->SetTitle("gTotCounts_927_BP_before");
    gTotCounts_927_BP_before->SetMarkerSize(1.4);
    gTotCounts_927_BP_before->SetMarkerColor(kBlack);
    gTotCounts_927_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_927_BP_after = new TGraphErrors(vRunFloat_927_BP_after.size(),&vRunFloat_927_BP_after[0],&vTotCountsLB_927_after_BP[0],NULL,NULL); //MT11
    gTotCounts_927_BP_after->SetName("gTotCounts_927_BP_after");
    gTotCounts_927_BP_after->SetTitle("gTotCounts_927_BP_after");
    gTotCounts_927_BP_after->SetMarkerSize(1.4);
    gTotCounts_927_BP_after->SetMarkerColor(kRed);
    gTotCounts_927_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_732_BP_before = new TGraphErrors(vRunFloat_732_BP_before.size(),&vRunFloat_732_BP_before[0],&vTotCountsLB_732_before_BP[0],NULL,NULL); //MT11
    gTotCounts_732_BP_before->SetName("gTotCounts_732_BP_before");
    gTotCounts_732_BP_before->SetTitle("gTotCounts_732_BP_before");
    gTotCounts_732_BP_before->SetMarkerSize(1.4);
    gTotCounts_732_BP_before->SetMarkerColor(kGreen+3);
    gTotCounts_732_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_732_BP_after = new TGraphErrors(vRunFloat_732_BP_after.size(),&vRunFloat_732_BP_after[0],&vTotCountsLB_732_after_BP[0],NULL,NULL); //MT11
    gTotCounts_732_BP_after->SetName("gTotCounts_732_BP_before");
    gTotCounts_732_BP_after->SetTitle("gTotCounts_732_BP_before");
    gTotCounts_732_BP_after->SetMarkerSize(1.4);
    gTotCounts_732_BP_after->SetMarkerColor(kBlue);
    gTotCounts_732_BP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_766_BP_before = new TGraphErrors(vRunFloat_766_BP_before.size(),&vRunFloat_766_BP_before[0],&vTotCountsLB_766_before_BP[0],NULL,NULL); //MT11
    gTotCounts_766_BP_before->SetName("gTotCounts_766_BP_before");
    gTotCounts_766_BP_before->SetTitle("gTotCounts_766_BP_before");
    gTotCounts_766_BP_before->SetMarkerSize(1.4);
    gTotCounts_766_BP_before->SetMarkerColor(kMagenta);
    gTotCounts_766_BP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_766_BP_after = new TGraphErrors(vRunFloat_766_BP_after.size(),&vRunFloat_766_BP_after[0],&vTotCountsLB_766_after_BP[0],NULL,NULL); //MT11
    gTotCounts_766_BP_after->SetName("gTotCounts_766_BP_before");
    gTotCounts_766_BP_after->SetTitle("gTotCounts_766_BP_before");
    gTotCounts_766_BP_after->SetMarkerSize(1.4);
    gTotCounts_766_BP_after->SetMarkerColor(kOrange);
    gTotCounts_766_BP_after->SetMarkerStyle(23);

    //NBP MT11: 225, 30, 64
    TGraphErrors *gTotCounts_225_NBP_before = new TGraphErrors(vRunFloat_225_NBP_before.size(),&vRunFloat_225_NBP_before[0],&vTotCountsLB_225_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_225_NBP_before->SetName("gTotCounts_225_NBP_before");
    gTotCounts_225_NBP_before->SetTitle("gTotCounts_225_NBP_before");
    gTotCounts_225_NBP_before->SetMarkerSize(1.4);
    gTotCounts_225_NBP_before->SetMarkerColor(kBlack);
    gTotCounts_225_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_225_NBP_after = new TGraphErrors(vRunFloat_225_NBP_after.size(),&vRunFloat_225_NBP_after[0],&vTotCountsLB_225_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_225_NBP_after->SetName("gTotCounts_225_NBP_before");
    gTotCounts_225_NBP_after->SetTitle("gTotCounts_225_NBP_before");
    gTotCounts_225_NBP_after->SetMarkerSize(1.4);
    gTotCounts_225_NBP_after->SetMarkerColor(kRed);
    gTotCounts_225_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_30_NBP_before = new TGraphErrors(vRunFloat_30_NBP_before.size(),&vRunFloat_30_NBP_before[0],&vTotCountsLB_30_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_30_NBP_before->SetName("gTotCounts_30_NBP_before");
    gTotCounts_30_NBP_before->SetTitle("gTotCounts_30_NBP_before");
    gTotCounts_30_NBP_before->SetMarkerSize(1.4);
    gTotCounts_30_NBP_before->SetMarkerColor(kGreen+3);
    gTotCounts_30_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_30_NBP_after = new TGraphErrors(vRunFloat_30_NBP_after.size(),&vRunFloat_30_NBP_after[0],&vTotCountsLB_30_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_30_NBP_after->SetName("gTotCounts_30_NBP_before");
    gTotCounts_30_NBP_after->SetTitle("gTotCounts_30_NBP_before");
    gTotCounts_30_NBP_after->SetMarkerSize(1.4);
    gTotCounts_30_NBP_after->SetMarkerColor(kBlue);
    gTotCounts_30_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_64_NBP_before = new TGraphErrors(vRunFloat_64_NBP_before.size(),&vRunFloat_64_NBP_before[0],&vTotCountsLB_64_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_64_NBP_before->SetName("gTotCounts_64_NBP_before");
    gTotCounts_64_NBP_before->SetTitle("gTotCounts_64_NBP_before");
    gTotCounts_64_NBP_before->SetMarkerSize(1.4);
    gTotCounts_64_NBP_before->SetMarkerColor(kMagenta);
    gTotCounts_64_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_64_NBP_after = new TGraphErrors(vRunFloat_64_NBP_after.size(),&vRunFloat_64_NBP_after[0],&vTotCountsLB_64_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_64_NBP_after->SetName("gTotCounts_64_NBP_before");
    gTotCounts_64_NBP_after->SetTitle("gTotCounts_64_NBP_before");
    gTotCounts_64_NBP_after->SetMarkerSize(1.4);
    gTotCounts_64_NBP_after->SetMarkerColor(kOrange);
    gTotCounts_64_NBP_after->SetMarkerStyle(23);

    //NBP planes MT12 459, 264, 298
    TGraphErrors *gTotCounts_459_NBP_before = new TGraphErrors(vRunFloat_459_NBP_before.size(),&vRunFloat_459_NBP_before[0],&vTotCountsLB_459_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_459_NBP_before->SetName("gTotCounts_459_NBP_before");
    gTotCounts_459_NBP_before->SetTitle("gTotCounts_459_NBP_before");
    gTotCounts_459_NBP_before->SetMarkerSize(1.4);
    gTotCounts_459_NBP_before->SetMarkerColor(kBlack);
    gTotCounts_459_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_459_NBP_after = new TGraphErrors(vRunFloat_459_NBP_after.size(),&vRunFloat_459_NBP_after[0],&vTotCountsLB_459_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_459_NBP_after->SetName("gTotCounts_459_NBP_after");
    gTotCounts_459_NBP_after->SetTitle("gTotCounts_459_NBP_after");
    gTotCounts_459_NBP_after->SetMarkerSize(1.4);
    gTotCounts_459_NBP_after->SetMarkerColor(kRed);
    gTotCounts_459_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_264_NBP_before = new TGraphErrors(vRunFloat_264_NBP_before.size(),&vRunFloat_264_NBP_before[0],&vTotCountsLB_264_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_264_NBP_before->SetName("gTotCounts_264_NBP_before");
    gTotCounts_264_NBP_before->SetTitle("gTotCounts_264_NBP_before");
    gTotCounts_264_NBP_before->SetMarkerSize(1.4);
    gTotCounts_264_NBP_before->SetMarkerColor(kGreen+3);
    gTotCounts_264_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_264_NBP_after = new TGraphErrors(vRunFloat_264_NBP_after.size(),&vRunFloat_264_NBP_after[0],&vTotCountsLB_264_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_264_NBP_after->SetName("gTotCounts_264_NBP_before");
    gTotCounts_264_NBP_after->SetTitle("gTotCounts_264_NBP_before");
    gTotCounts_264_NBP_after->SetMarkerSize(1.4);
    gTotCounts_264_NBP_after->SetMarkerColor(kBlue);
    gTotCounts_264_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_298_NBP_before = new TGraphErrors(vRunFloat_298_NBP_before.size(),&vRunFloat_298_NBP_before[0],&vTotCountsLB_298_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_298_NBP_before->SetName("gTotCounts_298_NBP_before");
    gTotCounts_298_NBP_before->SetTitle("gTotCounts_298_NBP_before");
    gTotCounts_298_NBP_before->SetMarkerSize(1.4);
    gTotCounts_298_NBP_before->SetMarkerColor(kMagenta);
    gTotCounts_298_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_298_NBP_after = new TGraphErrors(vRunFloat_298_NBP_after.size(),&vRunFloat_298_NBP_after[0],&vTotCountsLB_298_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_298_NBP_after->SetName("gTotCounts_298_NBP_before");
    gTotCounts_298_NBP_after->SetTitle("gTotCounts_298_NBP_before");
    gTotCounts_298_NBP_after->SetMarkerSize(1.4);
    gTotCounts_298_NBP_after->SetMarkerColor(kOrange);
    gTotCounts_298_NBP_after->SetMarkerStyle(23);

    //NBP planes MT21 693, 498, 532
    TGraphErrors *gTotCounts_693_NBP_before = new TGraphErrors(vRunFloat_693_NBP_before.size(),&vRunFloat_693_NBP_before[0],&vTotCountsLB_693_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_693_NBP_before->SetName("gTotCounts_693_NBP_before");
    gTotCounts_693_NBP_before->SetTitle("gTotCounts_693_NBP_before");
    gTotCounts_693_NBP_before->SetMarkerSize(1.4);
    gTotCounts_693_NBP_before->SetMarkerColor(kBlack);
    gTotCounts_693_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_693_NBP_after = new TGraphErrors(vRunFloat_693_NBP_after.size(),&vRunFloat_693_NBP_after[0],&vTotCountsLB_693_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_693_NBP_after->SetName("gTotCounts_693_NBP_after");
    gTotCounts_693_NBP_after->SetTitle("gTotCounts_693_NBP_after");
    gTotCounts_693_NBP_after->SetMarkerSize(1.4);
    gTotCounts_693_NBP_after->SetMarkerColor(kRed);
    gTotCounts_693_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_498_NBP_before = new TGraphErrors(vRunFloat_498_NBP_before.size(),&vRunFloat_498_NBP_before[0],&vTotCountsLB_498_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_498_NBP_before->SetName("gTotCounts_498_NBP_before");
    gTotCounts_498_NBP_before->SetTitle("gTotCounts_498_NBP_before");
    gTotCounts_498_NBP_before->SetMarkerSize(1.4);
    gTotCounts_498_NBP_before->SetMarkerColor(kGreen+3);
    gTotCounts_498_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_498_NBP_after = new TGraphErrors(vRunFloat_498_NBP_after.size(),&vRunFloat_498_NBP_after[0],&vTotCountsLB_498_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_498_NBP_after->SetName("gTotCounts_498_NBP_before");
    gTotCounts_498_NBP_after->SetTitle("gTotCounts_498_NBP_before");
    gTotCounts_498_NBP_after->SetMarkerSize(1.4);
    gTotCounts_498_NBP_after->SetMarkerColor(kBlue);
    gTotCounts_498_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_532_NBP_before = new TGraphErrors(vRunFloat_532_NBP_before.size(),&vRunFloat_532_NBP_before[0],&vTotCountsLB_532_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_532_NBP_before->SetName("gTotCounts_532_NBP_before");
    gTotCounts_532_NBP_before->SetTitle("gTotCounts_532_NBP_before");
    gTotCounts_532_NBP_before->SetMarkerSize(1.4);
    gTotCounts_532_NBP_before->SetMarkerColor(kMagenta);
    gTotCounts_532_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_532_NBP_after = new TGraphErrors(vRunFloat_532_NBP_after.size(),&vRunFloat_532_NBP_after[0],&vTotCountsLB_532_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_532_NBP_after->SetName("gTotCounts_532_NBP_before");
    gTotCounts_532_NBP_after->SetTitle("gTotCounts_532_NBP_before");
    gTotCounts_532_NBP_after->SetMarkerSize(1.4);
    gTotCounts_532_NBP_after->SetMarkerColor(kOrange);
    gTotCounts_532_NBP_after->SetMarkerStyle(23);

    //NBP planes MT22 927, 732, 766
    TGraphErrors *gTotCounts_927_NBP_before = new TGraphErrors(vRunFloat_927_NBP_before.size(),&vRunFloat_927_NBP_before[0],&vTotCountsLB_927_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_927_NBP_before->SetName("gTotCounts_927_NBP_before");
    gTotCounts_927_NBP_before->SetTitle("gTotCounts_927_NBP_before");
    gTotCounts_927_NBP_before->SetMarkerSize(1.4);
    gTotCounts_927_NBP_before->SetMarkerColor(kBlack);
    gTotCounts_927_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_927_NBP_after = new TGraphErrors(vRunFloat_927_NBP_after.size(),&vRunFloat_927_NBP_after[0],&vTotCountsLB_927_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_927_NBP_after->SetName("gTotCounts_927_NBP_after");
    gTotCounts_927_NBP_after->SetTitle("gTotCounts_927_NBP_after");
    gTotCounts_927_NBP_after->SetMarkerSize(1.4);
    gTotCounts_927_NBP_after->SetMarkerColor(kRed);
    gTotCounts_927_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_732_NBP_before = new TGraphErrors(vRunFloat_732_NBP_before.size(),&vRunFloat_732_NBP_before[0],&vTotCountsLB_732_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_732_NBP_before->SetName("gTotCounts_732_NBP_before");
    gTotCounts_732_NBP_before->SetTitle("gTotCounts_732_NBP_before");
    gTotCounts_732_NBP_before->SetMarkerSize(1.4);
    gTotCounts_732_NBP_before->SetMarkerColor(kGreen+3);
    gTotCounts_732_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_732_NBP_after = new TGraphErrors(vRunFloat_732_NBP_after.size(),&vRunFloat_732_NBP_after[0],&vTotCountsLB_732_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_732_NBP_after->SetName("gTotCounts_732_NBP_before");
    gTotCounts_732_NBP_after->SetTitle("gTotCounts_732_NBP_before");
    gTotCounts_732_NBP_after->SetMarkerSize(1.4);
    gTotCounts_732_NBP_after->SetMarkerColor(kBlue);
    gTotCounts_732_NBP_after->SetMarkerStyle(23);
    //
    TGraphErrors *gTotCounts_766_NBP_before = new TGraphErrors(vRunFloat_766_NBP_before.size(),&vRunFloat_766_NBP_before[0],&vTotCountsLB_766_before_NBP[0],NULL,NULL); //MT11
    gTotCounts_766_NBP_before->SetName("gTotCounts_766_NBP_before");
    gTotCounts_766_NBP_before->SetTitle("gTotCounts_766_NBP_before");
    gTotCounts_766_NBP_before->SetMarkerSize(1.4);
    gTotCounts_766_NBP_before->SetMarkerColor(kMagenta);
    gTotCounts_766_NBP_before->SetMarkerStyle(22);
    TGraphErrors *gTotCounts_766_NBP_after = new TGraphErrors(vRunFloat_766_NBP_after.size(),&vRunFloat_766_NBP_after[0],&vTotCountsLB_766_after_NBP[0],NULL,NULL); //MT11
    gTotCounts_766_NBP_after->SetName("gTotCounts_766_NBP_before");
    gTotCounts_766_NBP_after->SetTitle("gTotCounts_766_NBP_before");
    gTotCounts_766_NBP_after->SetMarkerSize(1.4);
    gTotCounts_766_NBP_after->SetMarkerColor(kOrange);
    gTotCounts_766_NBP_after->SetMarkerStyle(23);

    //Multigraph both
    //MT11 225, 30, 64
    TMultiGraph *mEffLB_225_Both = new TMultiGraph();
    mEffLB_225_Both->Add(gEffLB_225_Both_before);
    mEffLB_225_Both->Add(gEffLB_225_Both_after);
    mEffLB_225_Both->Add(gTotCounts_225_Both_before);
    mEffLB_225_Both->Add(gTotCounts_225_Both_after);

    TCanvas *cEffLB_225_Both = new TCanvas();
    cEffLB_225_Both->cd();
    mEffLB_225_Both->SetTitle("LB 225 both planes MT11");
    mEffLB_225_Both->Draw("AP");
    mEffLB_225_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_225_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_225_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_225_Both->GetYaxis()->SetTitleFont(62);
    mEffLB_225_Both->GetXaxis()->SetLabelFont(62);
    mEffLB_225_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_225_Both->GetXaxis()->CenterTitle(true);
    mEffLB_225_Both->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_30_Both = new TMultiGraph();
    mEffLB_30_Both->Add(gEffLB_30_Both_before);
    mEffLB_30_Both->Add(gEffLB_30_Both_after);
    mEffLB_30_Both->Add(gTotCounts_30_Both_before);
    mEffLB_30_Both->Add(gTotCounts_30_Both_after);

    TCanvas *cEffLB_30_Both = new TCanvas();
    cEffLB_30_Both->cd();
    mEffLB_30_Both->SetTitle("LB 30 both planes MT11");
    mEffLB_30_Both->Draw("AP");
    mEffLB_30_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_30_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_30_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_30_Both->GetYaxis()->SetTitleFont(62);
    mEffLB_30_Both->GetXaxis()->SetLabelFont(62);
    mEffLB_30_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_30_Both->GetXaxis()->CenterTitle(true);
    mEffLB_30_Both->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_64_Both = new TMultiGraph();
    mEffLB_64_Both->Add(gEffLB_64_Both_before);
    mEffLB_64_Both->Add(gEffLB_64_Both_after);
    mEffLB_64_Both->Add(gTotCounts_64_Both_before);
    mEffLB_64_Both->Add(gTotCounts_64_Both_after);

    TCanvas *cEffLB_64_Both = new TCanvas();
    cEffLB_64_Both->cd();
    mEffLB_64_Both->SetTitle("LB 64 both planes MT11");
    mEffLB_64_Both->Draw("AP");
    mEffLB_64_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_64_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_64_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_64_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_64_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_64_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_64_Both->GetXaxis()->CenterTitle(true);
    mEffLB_64_Both->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");

    //MT12 459, 264, 298
    TMultiGraph *mEffLB_459_Both = new TMultiGraph();
    mEffLB_459_Both->Add(gEffLB_459_Both_before);
    mEffLB_459_Both->Add(gEffLB_459_Both_after);
    mEffLB_459_Both->Add(gTotCounts_459_Both_before);
    mEffLB_459_Both->Add(gTotCounts_459_Both_after);

    TCanvas *cEffLB_459_Both = new TCanvas();
    cEffLB_459_Both->cd();
    mEffLB_459_Both->SetTitle("LB 459 both planes MT12");
    mEffLB_459_Both->Draw("AP");
    mEffLB_459_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_459_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_459_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_459_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_459_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_459_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_459_Both->GetXaxis()->CenterTitle(true);
    mEffLB_459_Both->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_264_Both = new TMultiGraph();
    mEffLB_264_Both->Add(gEffLB_264_Both_before);
    mEffLB_264_Both->Add(gEffLB_264_Both_after);
    mEffLB_264_Both->Add(gTotCounts_264_Both_before);
    mEffLB_264_Both->Add(gTotCounts_264_Both_after);

    TCanvas *cEffLB_264_Both = new TCanvas();
    cEffLB_264_Both->cd();
    mEffLB_264_Both->SetTitle("LB 264 both planes MT12");
    mEffLB_264_Both->Draw("AP");
    mEffLB_264_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_264_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_264_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_264_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_264_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_264_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_264_Both->GetXaxis()->CenterTitle(true);
    mEffLB_264_Both->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_298_Both = new TMultiGraph();
    mEffLB_298_Both->Add(gEffLB_298_Both_before);
    mEffLB_298_Both->Add(gEffLB_298_Both_after);
    mEffLB_298_Both->Add(gTotCounts_298_Both_before);
    mEffLB_298_Both->Add(gTotCounts_298_Both_after);

    TCanvas *cEffLB_298_Both = new TCanvas();
    cEffLB_298_Both->cd();
    mEffLB_298_Both->SetTitle("LB 298 both planes MT12");
    mEffLB_298_Both->Draw("AP");
    mEffLB_298_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_298_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_298_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_298_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_298_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_298_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_298_Both->GetXaxis()->CenterTitle(true);
    mEffLB_298_Both->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");


    //MT21 693, 498, 532
    TMultiGraph *mEffLB_693_Both = new TMultiGraph();
    mEffLB_693_Both->Add(gEffLB_693_Both_before);
    mEffLB_693_Both->Add(gEffLB_693_Both_after);
    mEffLB_693_Both->Add(gTotCounts_693_Both_before);
    mEffLB_693_Both->Add(gTotCounts_693_Both_after);

    TCanvas *cEffLB_693_Both = new TCanvas();
    cEffLB_693_Both->cd();
    mEffLB_693_Both->SetTitle("LB 693 both planes MT21");
    mEffLB_693_Both->Draw("AP");
    mEffLB_693_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_693_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_693_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_693_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_693_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_693_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_693_Both->GetXaxis()->CenterTitle(true);
    mEffLB_693_Both->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_498_Both = new TMultiGraph();
    mEffLB_498_Both->Add(gEffLB_498_Both_before);
    mEffLB_498_Both->Add(gEffLB_498_Both_after);
    mEffLB_498_Both->Add(gTotCounts_498_Both_before);
    mEffLB_498_Both->Add(gTotCounts_498_Both_after);

    TCanvas *cEffLB_498_Both = new TCanvas();
    cEffLB_498_Both->cd();
    mEffLB_498_Both->SetTitle("LB 498 both planes MT21");
    mEffLB_498_Both->Draw("AP");
    mEffLB_498_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_498_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_498_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_498_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_498_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_498_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_498_Both->GetXaxis()->CenterTitle(true);
    mEffLB_498_Both->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_532_Both = new TMultiGraph();
    mEffLB_532_Both->Add(gEffLB_532_Both_before);
    mEffLB_532_Both->Add(gEffLB_532_Both_after);
    mEffLB_532_Both->Add(gTotCounts_532_Both_before);
    mEffLB_532_Both->Add(gTotCounts_532_Both_after);

    TCanvas *cEffLB_532_Both = new TCanvas();
    cEffLB_532_Both->cd();
    mEffLB_532_Both->SetTitle("LB 532 both planes MT21");
    mEffLB_532_Both->Draw("AP");
    mEffLB_532_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_532_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_532_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_532_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_532_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_532_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_532_Both->GetXaxis()->CenterTitle(true);
    mEffLB_532_Both->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");

    //MT22 927, 732, 766
    TMultiGraph *mEffLB_927_Both = new TMultiGraph();
    mEffLB_927_Both->Add(gEffLB_927_Both_before);
    mEffLB_927_Both->Add(gEffLB_927_Both_after);
    mEffLB_927_Both->Add(gTotCounts_927_Both_before);
    mEffLB_927_Both->Add(gTotCounts_927_Both_after);
    
    TCanvas *cEffLB_927_Both = new TCanvas();
    cEffLB_927_Both->cd();
    mEffLB_927_Both->SetTitle("LB 927 both planes MT22");
    mEffLB_927_Both->Draw("AP");
    mEffLB_927_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_927_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_927_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_927_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_927_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_927_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_927_Both->GetXaxis()->CenterTitle(true);
    mEffLB_927_Both->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_732_Both = new TMultiGraph();
    mEffLB_732_Both->Add(gEffLB_732_Both_before);
    mEffLB_732_Both->Add(gEffLB_732_Both_after);
    mEffLB_732_Both->Add(gTotCounts_732_Both_before);
    mEffLB_732_Both->Add(gTotCounts_732_Both_after);

    TCanvas *cEffLB_732_Both = new TCanvas();
    cEffLB_732_Both->cd();
    mEffLB_732_Both->SetTitle("LB 732 both planes MT22");
    mEffLB_732_Both->Draw("AP");
    mEffLB_732_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_732_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_732_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_732_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_732_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_732_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_732_Both->GetXaxis()->CenterTitle(true);
    mEffLB_732_Both->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_766_Both = new TMultiGraph();
    mEffLB_766_Both->Add(gEffLB_766_Both_before);
    mEffLB_766_Both->Add(gEffLB_766_Both_after);
    mEffLB_766_Both->Add(gTotCounts_766_Both_before);
    mEffLB_766_Both->Add(gTotCounts_766_Both_after);

    TCanvas *cEffLB_766_Both = new TCanvas();
    cEffLB_766_Both->cd();
    mEffLB_766_Both->SetTitle("LB 766 both planes MT22");
    mEffLB_766_Both->Draw("AP");
    mEffLB_766_Both->GetXaxis()->SetTitle("Run #");
    mEffLB_766_Both->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_766_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_766_Both->GetXaxis()->SetTitleFont(62);
    mEffLB_766_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_766_Both->GetYaxis()->SetLabelFont(62);
    mEffLB_766_Both->GetXaxis()->CenterTitle(true);
    mEffLB_766_Both->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");

    //Multigraph BP
    //MT11 225, 30, 64
    TMultiGraph *mEffLB_225_BP = new TMultiGraph();
    mEffLB_225_BP->Add(gEffLB_225_BP_before);
    mEffLB_225_BP->Add(gEffLB_225_BP_after);
    mEffLB_225_BP->Add(gTotCounts_225_BP_before);
    mEffLB_225_BP->Add(gTotCounts_225_BP_after);

    TCanvas *cEffLB_225_BP = new TCanvas();
    cEffLB_225_BP->cd();
    mEffLB_225_BP->SetTitle("LB 225 BP planes MT11");
    mEffLB_225_BP->Draw("AP");
    mEffLB_225_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_225_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_225_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_225_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_225_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_225_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_225_BP->GetXaxis()->CenterTitle(true);
    mEffLB_225_BP->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mTotCountsLB_30_BP = new TMultiGraph();
    mTotCountsLB_30_BP->Add(gTotCounts_30_BP_before);
    mTotCountsLB_30_BP->Add(gTotCounts_30_BP_after);
    TMultiGraph *mEffLB_30_BP = new TMultiGraph();
    mEffLB_30_BP->Add(gEffLB_30_BP_before);
    mEffLB_30_BP->Add(gEffLB_30_BP_after);
    //mEffLB_30_BP->Add(gTotCounts_30_BP_before);
    //mEffLB_30_BP->Add(gTotCounts_30_BP_after);
    mEffLB_30_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_30_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_30_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_30_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_30_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_30_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_30_BP->GetXaxis()->CenterTitle(true);
    mEffLB_30_BP->GetYaxis()->CenterTitle(true);
    mEffLB_30_BP->SetTitle("LB 30 BP planes MT11");
    //mEffLB_30_BP->Draw("AP");
    
    //lLB30->Draw("SAME");

    TCanvas *cEffLB_30_BP = new TCanvas();
    cEffLB_30_BP->cd();
    mTotCountsLB_30_BP->Draw("AP");
    lLB30->Draw("SAME");

    /*TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
  	p2->SetGrid();
  	TPad *p3 = new TPad("p3", "", 0, 0, 1, 1);
 	p3->SetFillStyle(4000); // will be transparent
  
  	p2->Draw();
  	p2->cd();
    mEffLB_30_BP->Draw("AP");
  	gPad->Update();
  
	double xmin = p2->GetUxmin();
	double xmax = p2->GetUxmax();
	double dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
    mTotCountsLB_30_BP->GetYaxis()->SetRangeUser(0,200);
	double ymin = mTotCountsLB_30_BP->GetHistogram()->GetMinimum();
	double ymax = mTotCountsLB_30_BP->GetHistogram()->GetMaximum() + 0.5*(mTotCountsLB_30_BP->GetHistogram()->GetMaximum());
	double dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p3->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p3->Draw();
	p3->cd();
	mTotCountsLB_30_BP->Draw("P");
	gPad->Update();
  
    TGaxis *axisCounts_LB_30_BP = new TGaxis(xmax, 0, xmax, 200, 0, 200, 510, "+L");

    axisCounts_LB_30_BP->SetLineColor(kBlack);
  	axisCounts_LB_30_BP->SetLabelColor(kBlack);
	axisCounts_LB_30_BP->SetTitle("Total counts");
	axisCounts_LB_30_BP->SetTitleOffset(1.2);
	axisCounts_LB_30_BP->SetLabelFont(62);
	axisCounts_LB_30_BP->SetTextFont(62);
	axisCounts_LB_30_BP->CenterTitle(true);
  	axisCounts_LB_30_BP->Draw();
	lLB30->Draw("SAME");
  	gPad->Update();*/

    //
    TMultiGraph *mEffLB_64_BP = new TMultiGraph();
    mEffLB_64_BP->Add(gEffLB_64_BP_before);
    mEffLB_64_BP->Add(gEffLB_64_BP_after);
    mEffLB_64_BP->Add(gTotCounts_64_BP_before);
    mEffLB_64_BP->Add(gTotCounts_64_BP_after);

    TCanvas *cEffLB_64_BP = new TCanvas();
    cEffLB_64_BP->cd();
    mEffLB_64_BP->SetTitle("LB 64 BP planes MT11");
    mEffLB_64_BP->Draw("AP");
    mEffLB_64_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_64_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_64_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_64_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_64_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_64_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_64_BP->GetXaxis()->CenterTitle(true);
    mEffLB_64_BP->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");

    //MT12 459, 264, 298
    TMultiGraph *mEffLB_459_BP = new TMultiGraph();
    mEffLB_459_BP->Add(gEffLB_459_BP_before);
    mEffLB_459_BP->Add(gEffLB_459_BP_after);
    mEffLB_459_BP->Add(gTotCounts_459_BP_before);
    mEffLB_459_BP->Add(gTotCounts_459_BP_after);

    TCanvas *cEffLB_459_BP = new TCanvas();
    cEffLB_459_BP->cd();
    mEffLB_459_BP->SetTitle("LB 459 BP planes MT12");
    mEffLB_459_BP->Draw("AP");
    mEffLB_459_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_459_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_459_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_459_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_459_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_459_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_459_BP->GetXaxis()->CenterTitle(true);
    mEffLB_459_BP->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_264_BP = new TMultiGraph();
    mEffLB_264_BP->Add(gEffLB_264_BP_before);
    mEffLB_264_BP->Add(gEffLB_264_BP_after);
    mEffLB_264_BP->Add(gTotCounts_264_BP_before);
    mEffLB_264_BP->Add(gTotCounts_264_BP_after);

    TCanvas *cEffLB_264_BP = new TCanvas();
    cEffLB_264_BP->cd();
    mEffLB_264_BP->SetTitle("LB 264 BP planes MT12");
    mEffLB_264_BP->Draw("AP");
    mEffLB_264_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_264_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_264_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_264_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_264_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_264_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_264_BP->GetXaxis()->CenterTitle(true);
    mEffLB_264_BP->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_298_BP = new TMultiGraph();
    mEffLB_298_BP->Add(gEffLB_298_BP_before);
    mEffLB_298_BP->Add(gEffLB_298_BP_after);
    mEffLB_298_BP->Add(gTotCounts_298_BP_before);
    mEffLB_298_BP->Add(gTotCounts_298_BP_after);

    TCanvas *cEffLB_298_BP = new TCanvas();
    cEffLB_298_BP->cd();
    mEffLB_298_BP->SetTitle("LB 298 BP planes MT12");
    mEffLB_298_BP->Draw("AP");
    mEffLB_298_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_298_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_298_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_298_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_298_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_298_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_298_BP->GetXaxis()->CenterTitle(true);
    mEffLB_298_BP->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");


    //MT21 693, 498, 532
    TMultiGraph *mEffLB_693_BP = new TMultiGraph();
    mEffLB_693_BP->Add(gEffLB_693_BP_before);
    mEffLB_693_BP->Add(gEffLB_693_BP_after);
    mEffLB_693_BP->Add(gTotCounts_693_BP_before);
    mEffLB_693_BP->Add(gTotCounts_693_BP_after);

    TCanvas *cEffLB_693_BP = new TCanvas();
    cEffLB_693_BP->cd();
    mEffLB_693_BP->SetTitle("LB 693 BP planes MT21");
    mEffLB_693_BP->Draw("AP");
    mEffLB_693_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_693_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_693_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_693_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_693_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_693_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_693_BP->GetXaxis()->CenterTitle(true);
    mEffLB_693_BP->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_498_BP = new TMultiGraph();
    mEffLB_498_BP->Add(gEffLB_498_BP_before);
    mEffLB_498_BP->Add(gEffLB_498_BP_after);
    mEffLB_498_BP->Add(gTotCounts_498_BP_before);
    mEffLB_498_BP->Add(gTotCounts_498_BP_after);

    TCanvas *cEffLB_498_BP = new TCanvas();
    cEffLB_498_BP->cd();
    mEffLB_498_BP->SetTitle("LB 498 BP planes MT21");
    mEffLB_498_BP->Draw("AP");
    mEffLB_498_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_498_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_498_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_498_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_498_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_498_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_498_BP->GetXaxis()->CenterTitle(true);
    mEffLB_498_BP->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_532_BP = new TMultiGraph();
    mEffLB_532_BP->Add(gEffLB_532_BP_before);
    mEffLB_532_BP->Add(gEffLB_532_BP_after);
    mEffLB_532_BP->Add(gTotCounts_532_BP_before);
    mEffLB_532_BP->Add(gTotCounts_532_BP_after);

    TCanvas *cEffLB_532_BP = new TCanvas();
    cEffLB_532_BP->cd();
    mEffLB_532_BP->SetTitle("LB 532 BP planes MT21");
    mEffLB_532_BP->Draw("AP");
    mEffLB_532_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_532_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_532_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_532_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_532_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_532_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_532_BP->GetXaxis()->CenterTitle(true);
    mEffLB_532_BP->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");

    //MT22 927, 732, 766
    TMultiGraph *mEffLB_927_BP = new TMultiGraph();
    mEffLB_927_BP->Add(gEffLB_927_BP_before);
    mEffLB_927_BP->Add(gEffLB_927_BP_after);
    mEffLB_927_BP->Add(gTotCounts_927_BP_before);
    mEffLB_927_BP->Add(gTotCounts_927_BP_after);

    TCanvas *cEffLB_927_BP = new TCanvas();
    cEffLB_927_BP->cd();
    mEffLB_927_BP->SetTitle("LB 927 BP planes MT22");
    mEffLB_927_BP->Draw("AP");
    mEffLB_927_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_927_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_927_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_927_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_927_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_927_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_927_BP->GetXaxis()->CenterTitle(true);
    mEffLB_927_BP->GetYaxis()->CenterTitle(true);
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_732_BP = new TMultiGraph();
    mEffLB_732_BP->Add(gEffLB_732_BP_before);
    mEffLB_732_BP->Add(gEffLB_732_BP_after);
    mEffLB_732_BP->Add(gTotCounts_732_BP_before);
    mEffLB_732_BP->Add(gTotCounts_732_BP_after);

    TCanvas *cEffLB_732_BP = new TCanvas();
    cEffLB_732_BP->cd();
    mEffLB_732_BP->SetTitle("LB 732 BP planes MT22");
    mEffLB_732_BP->Draw("AP");
    mEffLB_732_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_732_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_732_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_732_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_732_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_732_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_732_BP->GetXaxis()->CenterTitle(true);
    mEffLB_732_BP->GetYaxis()->CenterTitle(true);
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_766_BP = new TMultiGraph();
    mEffLB_766_BP->Add(gEffLB_766_BP_before);
    mEffLB_766_BP->Add(gEffLB_766_BP_after);
    mEffLB_766_BP->Add(gTotCounts_766_BP_before);
    mEffLB_766_BP->Add(gTotCounts_766_BP_after);

    TCanvas *cEffLB_766_BP = new TCanvas();
    cEffLB_766_BP->cd();
    mEffLB_766_BP->SetTitle("LB 766 BP planes MT22");
    mEffLB_766_BP->Draw("AP");
    mEffLB_766_BP->GetXaxis()->SetTitle("Run #");
    mEffLB_766_BP->GetYaxis()->SetTitle("Efficiency [%]");
    mEffLB_766_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_766_BP->GetXaxis()->SetTitleFont(62);
    mEffLB_766_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_766_BP->GetYaxis()->SetLabelFont(62);
    mEffLB_766_BP->GetXaxis()->CenterTitle(true);
    mEffLB_766_BP->GetYaxis()->CenterTitle(true);
    lLB64->Draw("SAME");

    //Multigraph NBP
    //MT11 225, 30, 64
    TMultiGraph *mEffLB_225_NBP = new TMultiGraph();
    mEffLB_225_NBP->Add(gEffLB_225_NBP_before);
    mEffLB_225_NBP->Add(gEffLB_225_NBP_after);
    mEffLB_225_NBP->Add(gTotCounts_225_NBP_before);
    mEffLB_225_NBP->Add(gTotCounts_225_NBP_after);

    TCanvas *cEffLB_225_NBP = new TCanvas();
    cEffLB_225_NBP->cd();
    mEffLB_225_NBP->SetTitle("LB 225 NBP planes MT11");
    mEffLB_225_NBP->Draw("AP");
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_30_NBP = new TMultiGraph();
    mEffLB_30_NBP->Add(gEffLB_30_NBP_before);
    mEffLB_30_NBP->Add(gEffLB_30_NBP_after);
    mEffLB_30_NBP->Add(gTotCounts_30_NBP_before);
    mEffLB_30_NBP->Add(gTotCounts_30_NBP_after);

    TCanvas *cEffLB_30_NBP = new TCanvas();
    cEffLB_30_NBP->cd();
    mEffLB_30_NBP->SetTitle("LB 30 NBP planes MT11");
    mEffLB_30_NBP->Draw("AP");
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_64_NBP = new TMultiGraph();
    mEffLB_64_NBP->Add(gEffLB_64_NBP_before);
    mEffLB_64_NBP->Add(gEffLB_64_NBP_after);
    mEffLB_64_NBP->Add(gTotCounts_64_NBP_before);
    mEffLB_64_NBP->Add(gTotCounts_64_NBP_after);

    TCanvas *cEffLB_64_NBP = new TCanvas();
    cEffLB_64_NBP->cd();
    mEffLB_64_NBP->SetTitle("LB 64 NBP planes MT11");
    mEffLB_64_NBP->Draw("AP");
    lLB64->Draw("SAME");

    //MT12 459, 264, 298
    TMultiGraph *mEffLB_459_NBP = new TMultiGraph();
    mEffLB_459_NBP->Add(gEffLB_459_NBP_before);
    mEffLB_459_NBP->Add(gEffLB_459_NBP_after);
    mEffLB_459_NBP->Add(gTotCounts_459_NBP_before);
    mEffLB_459_NBP->Add(gTotCounts_459_NBP_after);

    TCanvas *cEffLB_459_NBP = new TCanvas();
    cEffLB_459_NBP->cd();
    mEffLB_459_NBP->SetTitle("LB 459 NBP planes MT12");
    mEffLB_459_NBP->Draw("AP");
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_264_NBP = new TMultiGraph();
    mEffLB_264_NBP->Add(gEffLB_264_NBP_before);
    mEffLB_264_NBP->Add(gEffLB_264_NBP_after);
    mEffLB_264_NBP->Add(gTotCounts_264_NBP_before);
    mEffLB_264_NBP->Add(gTotCounts_264_NBP_after);

    TCanvas *cEffLB_264_NBP = new TCanvas();
    cEffLB_264_NBP->cd();
    mEffLB_264_NBP->SetTitle("LB 264 NBP planes MT12");
    mEffLB_264_NBP->Draw("AP");
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_298_NBP = new TMultiGraph();
    mEffLB_298_NBP->Add(gEffLB_298_NBP_before);
    mEffLB_298_NBP->Add(gEffLB_298_NBP_after);
    mEffLB_298_NBP->Add(gTotCounts_298_NBP_before);
    mEffLB_298_NBP->Add(gTotCounts_298_NBP_after);

    TCanvas *cEffLB_298_NBP = new TCanvas();
    cEffLB_298_NBP->cd();
    mEffLB_298_NBP->SetTitle("LB 298 NBP planes MT12");
    mEffLB_298_NBP->Draw("AP");
    lLB64->Draw("SAME");

    //MT21 693, 498, 532
    TMultiGraph *mEffLB_693_NBP = new TMultiGraph();
    mEffLB_693_NBP->Add(gEffLB_693_NBP_before);
    mEffLB_693_NBP->Add(gEffLB_693_NBP_after);
    mEffLB_693_NBP->Add(gTotCounts_693_NBP_before);
    mEffLB_693_NBP->Add(gTotCounts_693_NBP_after);

    TCanvas *cEffLB_693_NBP = new TCanvas();
    cEffLB_693_NBP->cd();
    mEffLB_693_NBP->SetTitle("LB 693 NBP planes MT21");
    mEffLB_693_NBP->Draw("AP");
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_498_NBP = new TMultiGraph();
    mEffLB_498_NBP->Add(gEffLB_498_NBP_before);
    mEffLB_498_NBP->Add(gEffLB_498_NBP_after);
    mEffLB_498_NBP->Add(gTotCounts_498_NBP_before);
    mEffLB_498_NBP->Add(gTotCounts_498_NBP_after);

    TCanvas *cEffLB_498_NBP = new TCanvas();
    cEffLB_498_NBP->cd();
    mEffLB_498_NBP->SetTitle("LB 498 NBP planes MT21");
    mEffLB_498_NBP->Draw("AP");
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_532_NBP = new TMultiGraph();
    mEffLB_532_NBP->Add(gEffLB_532_NBP_before);
    mEffLB_532_NBP->Add(gEffLB_532_NBP_after);
    mEffLB_532_NBP->Add(gTotCounts_532_NBP_before);
    mEffLB_532_NBP->Add(gTotCounts_532_NBP_after);

    TCanvas *cEffLB_532_NBP = new TCanvas();
    cEffLB_532_NBP->cd();
    mEffLB_532_NBP->SetTitle("LB 532 NBP planes MT21");
    mEffLB_532_NBP->Draw("AP");
    lLB64->Draw("SAME");

    //MT22 927, 732, 766
    TMultiGraph *mEffLB_927_NBP = new TMultiGraph();
    mEffLB_927_NBP->Add(gEffLB_927_NBP_before);
    mEffLB_927_NBP->Add(gEffLB_927_NBP_after);
    mEffLB_927_NBP->Add(gTotCounts_927_NBP_before);
    mEffLB_927_NBP->Add(gTotCounts_927_NBP_after);

    TCanvas *cEffLB_927_NBP = new TCanvas();
    cEffLB_927_NBP->cd();
    mEffLB_927_NBP->SetTitle("LB 927 NBP planes MT22");
    mEffLB_927_NBP->Draw("AP");
    lLB225->Draw("SAME");
    //
    TMultiGraph *mEffLB_732_NBP = new TMultiGraph();
    mEffLB_732_NBP->Add(gEffLB_732_NBP_before);
    mEffLB_732_NBP->Add(gEffLB_732_NBP_after);
    mEffLB_732_NBP->Add(gTotCounts_732_NBP_before);
    mEffLB_732_NBP->Add(gTotCounts_732_NBP_after);

    TCanvas *cEffLB_732_NBP = new TCanvas();
    cEffLB_732_NBP->cd();
    mEffLB_732_NBP->SetTitle("LB 732 NBP planes MT22");
    mEffLB_732_NBP->Draw("AP");
    lLB30->Draw("SAME");
    //
    TMultiGraph *mEffLB_766_NBP = new TMultiGraph();
    mEffLB_766_NBP->Add(gEffLB_766_NBP_before);
    mEffLB_766_NBP->Add(gEffLB_766_NBP_after);
    mEffLB_766_NBP->Add(gTotCounts_766_NBP_before);
    mEffLB_766_NBP->Add(gTotCounts_766_NBP_after);

    TCanvas *cEffLB_766_NBP = new TCanvas();
    cEffLB_766_NBP->cd();
    mEffLB_766_NBP->SetTitle("LB 766 NBP planes MT22");
    mEffLB_766_NBP->Draw("AP");
    lLB64->Draw("SAME");


} //End of main function