#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TBox.h"
#include "TEfficiency.h"

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameter
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

//To create/delete folders and/or check they exist
namespace fs = std::filesystem;

const char ext[20] =".root";

void printEff() {

    vector<float> vErrEffBP, vErrEffNBP, vErrEffBoth;
    vector<float> vErrErrEffBP, vErrErrEffNBP, vErrErrEffBoth;
    vector<float> vFirstRun, vLastRun, vRun;

    //For nice plots
    vector<TBox*>  boxesBP, boxesNBP, boxesBoth;  // error boxes

    double yminBP = 1e+9, ymaxBP = 1e-9, yminNBP = 1e+9, ymaxNBP = 1e-9, yminBoth = 1e+9, ymaxBoth = 1e-9;

    //Enter in the scan folder 
    //string period = "LHC23_pass4_skimmed";
    //string period = "LHC23_PbPb_pass4";
    string period = "LHC25ad_pass2";
    //string period = "LHC25ae_pass2";
    //string period = "LHC24_ppref_pass1"; //pp ref 2024

    int trackGoal = 0;
    //int trackGoal = 300000000;

    //Crate sub-folder "merged" within period subfolder
	string fol = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/" + period + "/ccdb/merged/";
    
    if (!fs::exists(fol.c_str() )) {
        fs::create_directories(fol.c_str());
    }
    
    //Add track goal number to folder
    fol += to_string(trackGoal) + "/";
    
    if (!fs::exists(fol.c_str() )) {
        fs::create_directories(fol.c_str());
    }

    //Output files
    TFile *fOut = new TFile((fol + "outFile_merged_" + period + "_" + to_string(trackGoal)+".root").c_str() ,"RECREATE");

	gSystem->cd(fol.c_str()); //Enter folder

	//Count how many CAEN.root files are in the folder
	const char* entry; 
	const char* filename;
    int file_count = 0;
    TString str;

    //Vector of tuple to store the name of the file and first/last run number to order the file
    //otherwise they are not opened in the normal order (time wise) but according to the
    //binary conversion, strange things
    vector<tuple<TString,int,int>> files;
      
    char* dir = gSystem->ExpandPathName(fol.c_str());
    void* dirp = gSystem->OpenDirectory(dir);
    
    while((entry=gSystem->GetDirEntry(dirp))) {  	
        cout << entry << endl;

        TString fname(entry);
        TString tmp(entry);
    
        if (!tmp.EndsWith(".root")) continue;
        if (!tmp.BeginsWith("o2-mid-ChEffCounter_")) continue;

        // Remove prefix and extension
        tmp.ReplaceAll("o2-mid-ChEffCounter_", "");
        tmp.ReplaceAll(".root", "");

        // Split into x and y
        Ssiz_t pos = tmp.First('_');
        if (pos < 0) continue;
        int x = TString(tmp(0, pos)).Atoi();
        int y = TString(tmp(pos + 1, tmp.Length())).Atoi();

        cout << x << "\t" << y << endl;
        files.push_back(make_tuple(fname,x,y));
    }

    //Sort by first run number
    std::sort(files.begin(), files.end(), [](auto &a, auto &b){return std::get<1>(a) < std::get<1>(b);});

    cout << "After ordeing" << endl;

    //loop over sorted files
    for (auto &f : files) {

        TString entry = get<0>(f);
        int first = get<1>(f);
        int last = get<2>(f);

        vFirstRun.push_back(first);
        vLastRun.push_back(last);

        TH1F *hEffBP = new TH1F("hEffBP",("Efficiency error BP " +to_string(first) + "_" + to_string(last)).c_str(),30,0.,5.);
        TH1F *hEffNBP = new TH1F("hEffNBP",("Efficiency error NBP " +to_string(first) + "_" + to_string(last)).c_str(),30,0.,5.);
        TH1F *hEffBoth = new TH1F("hEffBoth",("Efficiency error Both " +to_string(first) + "_" + to_string(last)).c_str(),30,0.,5.);

        TH1F *hBP = new TH1F("hBP",("Efficiency BP " +to_string(first) + "_" + to_string(last)).c_str(),300,0.,105.);
        TH1F *hNBP = new TH1F("hNBP",("Efficiency NBP " +to_string(first) + "_" + to_string(last)).c_str(),300,0.,105.);
        TH1F *hBoth = new TH1F("hBoth",("Efficiency Both " +to_string(first) + "_" + to_string(last)).c_str(),300,0.,105.);

        std::cout << "Opening " << entry << " first: " << first << " last: " << last << std::endl;

        //Placeholder for two runs, done visually later
        vRun.push_back(file_count+0.5);
        double x = file_count+0.5;
        file_count++;

        TFile *fileEff = new TFile(entry,"READ"); //xxxxxx is the run number
        fileEff->cd();
        vector<o2::mid::ChEffCounter> *v = (vector<o2::mid::ChEffCounter>*)fileEff->Get("ccdb-object");

        float BPeff = 0., errBPeff = 0.;
        float NBPeff = 0., errNBPeff = 0.;
        float Botheff = 0., errBotheff = 0.;

        //Print and add to histo
        for (int i = 0; i < 936; i++) {
            //cout << "Board: " << i << endl;

            if  (v->at(i).getCounts(o2::mid::EffCountType::AllTracks) != 0) {
                
                //With TEfficiency class root
                /*double eff = (double)nBend / nAll;
                double effLow = TEfficiency::ClopperPearson(nAll, nBend, 0.6827, false);
                double effUp  = TEfficiency::ClopperPearson(nAll, nBend, 0.6827, true);

                double errLow = eff - effLow;
                double errUp  = effUp - eff;*/

                BPeff = (v->at(i).getCounts(o2::mid::EffCountType::BendPlane)/(float)v->at(i).getCounts(o2::mid::EffCountType::AllTracks))*100;
                errBPeff = TMath::Sqrt((BPeff*(100-BPeff))/(v->at(i).getCounts(o2::mid::EffCountType::AllTracks)));
                
                NBPeff = (v->at(i).getCounts(o2::mid::EffCountType::NonBendPlane)/(float)v->at(i).getCounts(o2::mid::EffCountType::AllTracks))*100;
                errNBPeff = TMath::Sqrt((NBPeff*(100-NBPeff))/(v->at(i).getCounts(o2::mid::EffCountType::AllTracks)));
                
                Botheff = (v->at(i).getCounts(o2::mid::EffCountType::BothPlanes)/(float)v->at(i).getCounts(o2::mid::EffCountType::AllTracks))*100; 
                errBotheff = TMath::Sqrt((Botheff*(100-Botheff))/(v->at(i).getCounts(o2::mid::EffCountType::AllTracks)));

                //Efficiency
                hBP->Fill(BPeff);
                hNBP->Fill(NBPeff);
                hBoth->Fill(Botheff);

                //Error on efficiency
                hEffBP->Fill(errBPeff);
                hEffNBP->Fill(errNBPeff);
                hEffBoth->Fill(errBotheff);
            }

            else {
                continue;
                //cout << "Board: " << i << " has 0 total counts" << endl; 
            }   
        }

        vErrEffBP.push_back(hEffBP->GetMean());
        vErrErrEffBP.push_back(hEffBP->GetMeanError());
        
        vErrEffNBP.push_back(hEffNBP->GetMean());
        vErrErrEffNBP.push_back(hEffNBP->GetMeanError());
        
        vErrEffBoth.push_back(hEffBoth->GetMean());
        vErrErrEffBoth.push_back(hEffBoth->GetMeanError());

        double dx = 0.3;
        //
        double meanBP = hEffBP->GetMean();
        double meanErrBP = hEffBP->GetStdDev();

        double meanNBP = hEffNBP->GetMean();
        double meanErrNBP = hEffNBP->GetStdDev();

        double meanBoth = hEffBoth->GetMean();
        double meanErrBoth = hEffBoth->GetStdDev();
        //
        yminBP = std::min(yminBP, meanBP - meanErrBP);
        ymaxBP = std::max(ymaxBP, meanBP + meanErrBP);

        yminNBP = std::min(yminNBP, meanNBP - meanErrNBP);
        ymaxNBP = std::max(ymaxNBP, meanNBP + meanErrNBP);

        yminBoth = std::min(yminBoth, meanBoth - meanErrBoth);
        ymaxBoth = std::max(ymaxBoth, meanBoth + meanErrBoth);
        //

        // Error box (centered on mean, height = ±meanErr)
        TBox* boxBP = new TBox(x - dx, meanBP - meanErrBP, x + dx, meanBP + meanErrBP);
        boxBP->SetFillColor(kRed - 9);
        boxBP->SetFillStyle(1001); // solid
        boxBP->SetLineColor(kRed);
        boxesBP.push_back(boxBP);

        TBox* boxNBP = new TBox(x - dx, meanNBP - meanErrNBP, x + dx, meanNBP + meanErrNBP);
        boxNBP->SetFillColor(kRed - 9);
        boxNBP->SetFillStyle(1001); // solid
        boxNBP->SetLineColor(kRed);
        boxesNBP.push_back(boxNBP);

        TBox* boxBoth = new TBox(x - dx, meanBoth - meanErrBoth, x + dx, meanBoth + meanErrBoth);
        boxBoth->SetFillColor(kRed - 9);
        boxBoth->SetFillStyle(1001); // solid
        boxBoth->SetLineColor(kRed);
        boxesBoth.push_back(boxBoth);

        fileEff->Close();
        delete fileEff;

        //Efficiency
        TCanvas *cEffBP = new TCanvas();
        cEffBP->cd();
        hBP->Draw("HISTO");
        hBP->GetXaxis()->SetTitle("LB efficiency [%] BP");
        hBP->GetYaxis()->SetTitle("Counts");

        TCanvas *cEffNBP = new TCanvas();
        cEffNBP->cd();
        hNBP->Draw("HISTO");
        hNBP->GetXaxis()->SetTitle("LB efficiency [%] NBP");
        hNBP->GetYaxis()->SetTitle("Counts");

        TCanvas *cEffBoth = new TCanvas();
        cEffBoth->cd();
        hBoth->Draw("HISTO");
        hBoth->GetXaxis()->SetTitle("LB efficiency [%] Both");
        hBoth->GetYaxis()->SetTitle("Counts");

        //Error on efficiency
        TCanvas *cErrEffBP = new TCanvas();
        cErrEffBP->cd();
        hEffBP->Draw("HISTO");
        hEffBP->GetXaxis()->SetTitle("LB error on efficiency [%] BP");
        hEffBP->GetYaxis()->SetTitle("Counts");

        TCanvas *cErrEffNBP = new TCanvas();
        cErrEffNBP->cd();
        hEffNBP->Draw("HISTO");
        hEffNBP->GetXaxis()->SetTitle("LB error on efficiency [%] NBP");
        hEffNBP->GetYaxis()->SetTitle("Counts");

        TCanvas *cErrEffBoth = new TCanvas();
        cErrEffBoth->cd();
        hEffBoth->Draw("HISTO");
        hEffBoth->GetXaxis()->SetTitle("LB error on efficiency [%] Both");
        hEffBoth->GetYaxis()->SetTitle("Counts");

        fOut->cd();
        //Efficiency
        cEffBP->Write(("BP_efficiency_" + to_string(first) + "_" + to_string(last)).c_str());
        cEffNBP->Write(("NBP_efficiency_" + to_string(first) + "_" + to_string(last)).c_str());
        cEffBoth->Write(("Both_efficiency_" + to_string(first) + "_" + to_string(last)).c_str());
        //Error on efficiency
        cErrEffBP->Write(("BP_" + to_string(first) + "_" + to_string(last)).c_str());
        cErrEffNBP->Write(("NBP_" + to_string(first) + "_" + to_string(last)).c_str());
        cErrEffBoth->Write(("Both_" + to_string(first) + "_" + to_string(last)).c_str());

        delete hEffBP;
        delete hEffNBP;
        delete hEffBoth;

        delete hBP;
        delete hNBP;
        delete hBoth;
    }

    //////
    //BP//
    //////
    TCanvas *cErrEffBPtrend = new TCanvas();
    cErrEffBPtrend->cd();
    cErrEffBPtrend->SetGridy();
    TGraphErrors *gErrEffBP = new TGraphErrors(vRun.size(),&vRun[0],&vErrEffBP[0],NULL,NULL);
    gErrEffBP->SetMarkerStyle(8);
    gErrEffBP->SetMarkerSize(1.4);
    gErrEffBP->GetXaxis()->SetTitle("Run #");
    gErrEffBP->GetXaxis()->SetTitleOffset(1.5);
    gErrEffBP->GetYaxis()->SetTitle("Error on efficiency BP");

    //Margin for y-axis
    double marginBP = 0.1 * (ymaxBP - yminBP);
    yminBP -= marginBP;
    ymaxBP += marginBP;
    gErrEffBP->GetYaxis()->SetRangeUser(yminBP, ymaxBP);

    //Draw
    gErrEffBP->Draw("AP");
    //Change x-axis labels
    for (auto box : boxesBP) box->Draw("SAME");

    for (int i = 0; i < gErrEffBP->GetN(); i++) { // Loop over all entries
        gErrEffBP->GetXaxis()->SetBinLabel(gErrEffBP->GetXaxis()->FindBin(i+0.5), Form("%d-%d", (int)vFirstRun[i], (int)vLastRun[i])); // Find out which bin on the x-axis the point corresponds to and set the bin label
    }
    TAxis *axBP = gErrEffBP->GetXaxis();
    //axBP->LabelsOption("d");  
    //Re-Draw graph on top of the error boxes
    gErrEffBP->Draw("P SAME");
    gPad->Update();

    ///////
    //NBP//
    ///////
    TCanvas *cErrEffNBPtrend = new TCanvas();
    cErrEffNBPtrend->cd();
    cErrEffNBPtrend->SetGridy();
    TGraphErrors *gErrEffNBP = new TGraphErrors(vRun.size(),&vRun[0],&vErrEffNBP[0],NULL,NULL);
    gErrEffNBP->SetMarkerStyle(8);
    gErrEffNBP->SetMarkerSize(1.4);
    gErrEffNBP->GetXaxis()->SetTitle("Run #");
    gErrEffNBP->GetXaxis()->SetTitleOffset(1.5);
    gErrEffNBP->GetYaxis()->SetTitle("Error on efficiency NBP");

    //Margin for y-axis
    double marginNBP = 0.1 * (ymaxNBP - yminNBP);
    yminNBP -= marginNBP;
    ymaxNBP += marginNBP;
    gErrEffNBP->GetYaxis()->SetRangeUser(yminNBP, ymaxNBP);

    //Draw
    gErrEffNBP->Draw("AP");
    //Change x-axis labels
    for (auto box : boxesNBP) box->Draw("SAME");

    for (int i = 0; i < gErrEffNBP->GetN(); i++) { // Loop over all entries
        gErrEffNBP->GetXaxis()->SetBinLabel(gErrEffNBP->GetXaxis()->FindBin(i+0.5), Form("%d-%d", (int)vFirstRun[i], (int)vLastRun[i])); // Find out which bin on the x-axis the point corresponds to and set the bin label
    }
    TAxis *axNBP = gErrEffNBP->GetXaxis();
    axNBP->LabelsOption("d");  
    //Re-Draw graph on top of the error boxes
    gErrEffNBP->Draw("P SAME");
    gPad->Update();

    ///////////////
    //Both planes//
    ///////////////
    TCanvas *cErrEffBothtrend = new TCanvas();
    cErrEffBothtrend->cd();
    cErrEffBothtrend->SetGridy();
    TGraphErrors *gErrEffBoth = new TGraphErrors(vRun.size(),&vRun[0],&vErrEffBoth[0],NULL,NULL);
    gErrEffBoth->SetMarkerStyle(8);
    gErrEffBoth->SetMarkerSize(1.4);
    gErrEffBoth->GetXaxis()->SetTitleOffset(1.5);
    gErrEffBoth->GetXaxis()->SetTitle("Run #");
    gErrEffBoth->GetYaxis()->SetTitle("Error on efficiency Both");

    //Margin for y-axis
    double marginBoth = 0.1 * (ymaxBoth - yminBoth);
    yminBoth -= marginBoth;
    ymaxBoth += marginBoth;
    gErrEffBoth->GetYaxis()->SetRangeUser(yminBoth, ymaxBoth);

    //Draw
    gErrEffBoth->Draw("AP");
    //Change x-axis labels
    for (auto box : boxesBoth) box->Draw("SAME");

    for (int i = 0; i < gErrEffBoth->GetN(); i++) { // Loop over all entries
        gErrEffBoth->GetXaxis()->SetBinLabel(gErrEffBoth->GetXaxis()->FindBin(i+0.5), Form("%d-%d", (int)vFirstRun[i], (int)vLastRun[i])); // Find out which bin on the x-axis the point corresponds to and set the bin label
    }
    TAxis *axBoth = gErrEffBoth->GetXaxis();
    axBoth->LabelsOption("d");  
    //Re-Draw graph on top of the error boxes
    gErrEffBoth->Draw("P SAME");
    gPad->Update();

    fOut->cd();
    cErrEffBPtrend->Write(("errEffTrend_BP_" + period).c_str());
    cErrEffNBPtrend->Write(("errEffTrend_NBP_" + period).c_str());
    cErrEffBothtrend->Write(("errEffTrend_Both_" + period).c_str());

    fOut->Close();
    delete fOut;

}