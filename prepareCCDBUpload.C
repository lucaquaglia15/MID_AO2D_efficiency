//Helper file used to count the total number of entries in a (reduced)AO2D file

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include "THStack.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TROOT.h"

void prepareCCDBUpload() {

    //File name: 
    //Start and end of run from RCT [STF - 120000:ETF + 120000] 2 minutes buffer, from "run_dates.txt"

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    //string period = "LHC22_pass7_skimmed";
    string period = "LHC23_pass4_skimmed";
    string globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/"; //Where txt file with star/end of run is stored

    string ccdbPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/"; //Where CCDB objects are stored

    ifstream hRun;
    hRun.open((globalPath+"run_IR_Bfield.txt").c_str());

    /*bool isIn;
    int run;
    long int start, end;
    vector<int> vRun;
    vector<long int> vStart, vEnd;

    while (runDates >> isIn >> run >> start >> end) {
        if (isIn) {
            vRun.push_back(run);
            vStart.push_back(start-120000);
            vEnd.push_back(end+120000);
        }
    }*/

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
            vStart.push_back(start-120000);
            vEnd.push_back(end+120000);
        }
    }

    ofstream hOut;
    hOut.open((ccdbPath+"upload.sh").c_str());
    hOut << "#!/bin/bash\n";
    hOut << "\n";
    hOut << "ccdbhost=http://alice-ccdb.cern.ch";
    hOut << "\n\n";

    for (unsigned int i = 0; i < vRun.size(); i++) {
        hOut << "o2-ccdb-upload --host \"$ccdbhost\" -p MID/Calib/ChamberEfficiency -f o2-mid-ChEffCounter_" + to_string(vRun.at(i)) + 
        ".root -k ccdb-object --starttimestamp " + to_string(vStart.at(i)) + " --endtimestamp " + to_string(vEnd.at(i)) + " -m \"runNumber=" + to_string(vRun.at(i)) + ";JIRA=O2-5759;\"\n";
    }
}