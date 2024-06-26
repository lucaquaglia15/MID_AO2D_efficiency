#include <iostream>
#include <vector>
#include <stdlib.h> 

#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

void fetchCCDB() {

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other dataset for the analyses of Nazar
    string globalPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/";

    //Path for the .txt file of the run list of the period
    string runNumbers = globalPath+"run_list.txt"; 
    string runDates = globalPath+"run_dates.txt";

    //Open txt file of runs
    ifstream hRun;
    hRun.open(runNumbers.c_str());

    ofstream hRunDates;
    hRunDates.open(runDates.c_str());

    //Push back to a vector of int (no need to care about size)
    int run;
    vector<int> vRun;

    while(hRun >> run) {
        vRun.push_back(run);
    }
    //sort in ascending order
    sort(vRun.begin(), vRun.end());

    o2::ccdb::CcdbApi api; //CCDB API
    api.init("http://alice-ccdb.cern.ch"); //Open connection to ALICE CCDB

    std::map<std::string, std::string> md; //Metada map
    std::map<std::string, std::string> headers; //Headers map

    for (unsigned iRun = 0; iRun < vRun.size(); iRun++) {
        headers = api.retrieveHeaders("RCT/Info/RunInformation", md, vRun.at(iRun));
        hRunDates << vRun.at(iRun) << "\t" << headers["STF"] << "\t" << headers["ETF"] << "\n";
        cout << vRun.at(iRun) << "\t" << headers["STF"] << "\t" << headers["ETF"] << "\n";
    }
        
    //string runToBeFetched = "RCT/Info/RunInformation/" + to_string(int(vRun[0]));
    //cout << runToBeFetched << endl;    
    
    //cout << headers["STF"] << endl;
    //cout << headers["ETF"] << endl;

    //hRunDates << vRun.at(0) << "\t" << headers["STF"] << "\t" << headers["ETF"] << "\n";
    
    hRun.close();
    hRunDates.close();

}