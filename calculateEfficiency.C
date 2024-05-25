#include "TH1F.h"
#include "TH2F.h"
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
#include "TKey.h"

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameters
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

#include "CCDB/CcdbApi.h" //CCDB api library

using namespace std;

void calculateEfficiency() {

    o2::mid::Mapping mapping;

    //string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC22o_pass6_small.root";
    string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC23_pass4_skimmed_QC1.root";

    TFile *fIn = new TFile(inFileName.c_str(),"READ");
    
    //string outFileName = "efficiency.root";
    //TFile *fOut = new TFile(outFileName.c_str(),"RECREATE");

    int nBinsPlane = 4;
    int nBinsRPC = 72;
    int nBinsBoard = 936;

    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    TH1F *hEffLBplanes1D_both[4]; //4 1D histograms for LB (one per plane)
    TH2F *hEffRPCplanes2D_both[4]; //4 2D histograms for RPC efficiency (one per plane)

    for (int i = 0; i < 4; i++) {
        hEffRPCplanes2D_both[i] = new TH2F(("RPC 2D efficiency "+planeName[i]).c_str(),("RPC 2D efficiency "+planeName[i]).c_str(),2,-1,1,9,0.5,9.5);
        hEffLBplanes1D_both[i] = new TH1F(("LB efficiency "+planeName[i]).c_str(),("LB efficiency "+planeName[i]).c_str(),234,0.5,243.5);
    }

    //TIter keyList(fIn->GetListOfKeys());
    //TKey *key;
    //while ((key = (TKey*)keyList())) {
        //key = (TKey*)keyList();
        //TClass *cl = gROOT->GetClass(key->GetClassName());
        //if (!cl->InheritsFrom("TDirectoryFile")) 
        //    continue;
        
        //TDirectoryFile *d = (TDirectoryFile*)key->ReadObj();
        TDirectoryFile *d = (TDirectoryFile*)fIn->Get("mid-efficiency");
        d->cd();
        
        //Get RPC counts
        TH1F *hFiredBothPlanesRPC = (TH1F*)d->Get("nFiredBothperRPC");
        TH1F *hFiredBPRPC= (TH1F*)d->Get("nFiredBPperRPC");
        TH1F *hFiredNBPRPC = (TH1F*)d->Get("nFiredNBPperRPC");
        TH1F *hTotRPC = (TH1F*)d->Get("nTotperRPC");

        //Get LB counts
        TH1F *hFiredBothPlanesLB = (TH1F*)d->Get("nFiredBothperBoard");
        TH1F *hFiredBPLB = (TH1F*)d->Get("nFiredBPperBoard");
        TH1F *hFiredNBPLB = (TH1F*)d->Get("nFiredNBPperBoard");
        TH1F *hTotLB = (TH1F*)d->Get("nTotperBoard");

        //Analyze RPCs
        TH1F *hEffRPC_both = new TH1F("effRPC_both","effRPC_both",nBinsRPC,-0.5,71.5);
        TH1F *hEffRPC_BP = new TH1F("effRPC_BP","effRPC_BP",nBinsRPC,-0.5,71.5);
        TH1F *hEffRPC_NBP = new TH1F("effRPC_NBP","effRPC_NBP",nBinsRPC,-0.5,71.5);

        //Analyze LB
        TH1F *hEffLB_both = new TH1F("effBoard_both","effBoard_both",nBinsBoard,0.5,936.5);
        TH1F *hEffLB_BP = new TH1F("effBoard_BP","effBoard_BP",nBinsBoard,0.5,936.5);
        TH1F *hEffLB_NBP = new TH1F("effBoard_NBP","effBoard_NBP",nBinsBoard,0.5,936.5);

        //fIn->Close();
        
        float effBothRPC, effBPRPC, effNBPRPC;
        float errEffBothRPC, errEffBPRPC, errEffNBPRPC;

        float effBothLB, effBPLB, effNBPLB;
        float errEffBothLB, errEffBPLB, errEffNBPLB;

        for (int i = 1; i <= nBinsRPC; i++) {
            cout << i << "\t" << hFiredBothPlanesRPC->GetBinContent(i) << "\t" << hTotRPC->GetBinContent(i) << endl;

            if (hTotRPC->GetBinContent(i) != 0) {
                effBothRPC = (hFiredBothPlanesRPC->GetBinContent(i)/hTotRPC->GetBinContent(i))*100;
                effBPRPC = (hFiredBPRPC->GetBinContent(i)/hTotRPC->GetBinContent(i))*100;
                effNBPRPC = (hFiredNBPRPC->GetBinContent(i)/hTotRPC->GetBinContent(i))*100;

                errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/hTotRPC->GetBinContent(i));
                errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/hTotRPC->GetBinContent(i));
                errEffNBPRPC = TMath::Sqrt(effNBPRPC*(100-effNBPRPC)/hTotRPC->GetBinContent(i));

                hEffRPC_both->SetBinContent(i,effBothRPC);
                hEffRPC_both->SetBinError(i,errEffBothRPC);
                hEffRPC_BP->SetBinContent(i,effBPRPC);
                hEffRPC_BP->SetBinError(i,errEffBPRPC);
                hEffRPC_NBP->SetBinContent(i,effNBPRPC);
                hEffRPC_NBP->SetBinError(i,errEffNBPRPC);
            }   

            //MT11
            if (i >= 1 && i <= 9) {
                hEffRPCplanes2D_both[0]->SetBinContent(2,i,effBothRPC);
            }
            else if (i >= 37 && i <= 45) {
                hEffRPCplanes2D_both[0]->SetBinContent(1,i-36,effBothRPC);
            }

            //MT12
            else if (i >= 10 && i <= 18) {
                hEffRPCplanes2D_both[1]->SetBinContent(2,i-9,effBothRPC);
            }
            else if (i >= 46 && i <= 54) {
                hEffRPCplanes2D_both[1]->SetBinContent(1,i-45,effBothRPC);
            }

            //MT21
            else if (i >= 19 && i <= 27) {
                hEffRPCplanes2D_both[2]->SetBinContent(2,i-18,effBothRPC);
            }
            else if (i >= 55 && i <= 63) {
                hEffRPCplanes2D_both[2]->SetBinContent(1,i-54,effBothRPC);
            }

            //MT22
            else if (i >= 28 && i <= 36) {
                hEffRPCplanes2D_both[3]->SetBinContent(2,i-27,effBothRPC);
            }
            else if (i >= 64 && i <= 72) {
                hEffRPCplanes2D_both[3]->SetBinContent(1,i-63,effBothRPC);
            }
        }

        for (int i = 1; i <= nBinsBoard; i++) {
            cout << "LB Both planes " <<  i << "\t" << hFiredBothPlanesLB->GetBinContent(i) << "\t" << hTotLB->GetBinContent(i) << endl;

            if (hTotLB->GetBinContent(i) != 0) {

                effBothLB = (hFiredBothPlanesLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effBPLB = (hFiredBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                effNBPLB = (hFiredNBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;

                errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLB->GetBinContent(i));
                errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLB->GetBinContent(i));
                errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLB->GetBinContent(i));

                hEffLB_both->Fill(i,effBothLB);
                hEffLB_both->SetBinError(i,errEffBothLB);
                hEffLB_BP->Fill(i,effBPLB);
                hEffLB_BP->SetBinError(i,errEffBPLB);
                hEffLB_NBP->Fill(i,effNBPLB);
                hEffLB_NBP->SetBinError(i,errEffNBPLB);

                if (i <= 234) {
                    hEffLBplanes1D_both[0]->SetBinContent(i,effBothLB);
                    hEffLBplanes1D_both[0]->SetBinError(i,errEffBothLB);
                }
                else if (i >= 235 && i <= 468) {
                    hEffLBplanes1D_both[1]->SetBinContent(i-234,effBothLB);
                    hEffLBplanes1D_both[1]->SetBinError(i-234,errEffBothLB);
                }
                else if (i>= 469 && i <= 702) {
                    hEffLBplanes1D_both[2]->SetBinContent(i-468,effBothLB);
                    hEffLBplanes1D_both[2]->SetBinError(i-468,errEffBothLB);
                }

                else {
                    hEffLBplanes1D_both[3]->SetBinContent(i-702,effBothLB);
                    hEffLBplanes1D_both[3]->SetBinError(i-702,errEffBothLB);
                }
            }
        }

        TCanvas *cEffBotPlanesRPC = new TCanvas();
        cEffBotPlanesRPC->cd();
        hEffRPC_both->SetStats(0);
        hEffRPC_both->GetXaxis()->SetTitle("RPC");
        hEffRPC_both->GetYaxis()->SetTitle("Efficiency [%]");
        hEffRPC_both->Draw("P");

        TCanvas *cEffBPRPC = new TCanvas();
        cEffBPRPC->cd();
        hEffRPC_BP->SetStats(0);
        hEffRPC_BP->GetXaxis()->SetTitle("RPC");
        hEffRPC_BP->GetYaxis()->SetTitle("Efficiency [%]");
        hEffRPC_BP->Draw("P");

        TCanvas *cEffNBPRPC = new TCanvas();
        cEffNBPRPC->cd();
        hEffRPC_NBP->SetStats(0);
        hEffRPC_NBP->GetXaxis()->SetTitle("RPC");
        hEffRPC_NBP->GetYaxis()->SetTitle("Efficiency [%]");
        hEffRPC_NBP->Draw("P");

        TCanvas *cEffBotPlanesBoard = new TCanvas();
        cEffBotPlanesBoard->cd();
        hEffLB_both->SetStats(0);
        hEffLB_both->GetXaxis()->SetTitle("Local Board");
        hEffLB_both->GetYaxis()->SetTitle("Efficiency [%]");
        hEffLB_both->Draw("HISTO");

        TCanvas *cEffBPBoard = new TCanvas();
        cEffBPBoard->cd();
        hEffLB_BP->SetStats(0);
        hEffLB_BP->GetXaxis()->SetTitle("Local Board");
        hEffLB_BP->GetYaxis()->SetTitle("Efficiency [%]");
        hEffLB_BP->Draw("HISTO");

        TCanvas *cEffNBPBoard = new TCanvas();
        cEffNBPBoard->cd();
        hEffLB_NBP->SetStats(0);
        hEffLB_NBP->GetXaxis()->SetTitle("Local Board");
        hEffLB_NBP->GetYaxis()->SetTitle("Efficiency [%]");
        hEffLB_NBP->Draw("HISTO");

        TCanvas *cEffLB_plane_both = new TCanvas();
        cEffLB_plane_both->Divide(2,2);
        cEffLB_plane_both->cd(1);
        hEffLBplanes1D_both[0]->SetStats(0);
        hEffLBplanes1D_both[0]->Draw("P");
        cEffLB_plane_both->cd(2);
        hEffLBplanes1D_both[1]->SetStats(0);
        hEffLBplanes1D_both[1]->Draw("P");
        cEffLB_plane_both->cd(3);
        hEffLBplanes1D_both[1]->SetStats(0);
        hEffLBplanes1D_both[2]->Draw("P");
        cEffLB_plane_both->cd(4);
        hEffLBplanes1D_both[1]->SetStats(0);
        hEffLBplanes1D_both[3]->Draw("P");

        TCanvas *cEffRPC_plane_both = new TCanvas();
        cEffRPC_plane_both->Divide(2,2);
        cEffRPC_plane_both->cd(1);
        hEffRPCplanes2D_both[0]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_both[0]->SetStats(0);
        hEffRPCplanes2D_both[0]->Draw("COLZ");
        cEffRPC_plane_both->cd(2);
        hEffRPCplanes2D_both[1]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_both[1]->SetStats(0);
        hEffRPCplanes2D_both[1]->Draw("COLZ");
        cEffRPC_plane_both->cd(3);
        hEffRPCplanes2D_both[2]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_both[2]->SetStats(0);
        hEffRPCplanes2D_both[2]->Draw("COLZ");
        cEffRPC_plane_both->cd(4);
        hEffRPCplanes2D_both[3]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_both[3]->SetStats(0);
        hEffRPCplanes2D_both[3]->Draw("COLZ");

        //Analyze the 4 periods of LHC23 pass4 skimmed QC1
        /*string periods[4] = {"za","zj","zs","zt"};

        TH1F *v_hEffLB_both[4];
        TH1F *v_hEffLB_BP[4];
        TH1F *v_hEffLB_NBP[4];

        THStack *hEffPeriodLB = new THStack();

        for (int period = 0; period < 4; period++) {
        
            string fileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC23"+periods[period]+".root";

            TFile *fInPeriod = new TFile(fileName.c_str(),"READ");
            fInPeriod->cd();

            TDirectoryFile *dPeriod = (TDirectoryFile*)fInPeriod->Get("mid-efficiency");
            dPeriod->cd();

            hFiredBothPlanesLB = (TH1F*)dPeriod->Get("nFiredBothperBoard");
            hFiredBPLB = (TH1F*)dPeriod->Get("nFiredBPperBoard");
            hFiredNBPLB = (TH1F*)dPeriod->Get("nFiredNBPperBoard");
            hTotLB = (TH1F*)dPeriod->Get("nTotperBoard");

            v_hEffLB_both[period] = new TH1F("effLB_both","effLB_both",nBinsBoard,-0.5,935.5);
            v_hEffLB_BP[period] = new TH1F("effLB_BP","effLB_BP",nBinsBoard,-0.5,935.5);
            v_hEffLB_NBP[period] = new TH1F("effLB_NBP","effLB_NBP",nBinsBoard,-0.5,935.5);

            for (int i = 1; i <= nBinsBoard; i++) {
                cout << "LB Both planes " <<  i << "\t" << hFiredBothPlanesLB->GetBinContent(i) << "\t" << hTotLB->GetBinContent(i) << endl;

                if (hTotLB->GetBinContent(i) != 0) {

                    effBothLB = (hFiredBothPlanesLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                    effBPLB = (hFiredBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;
                    effNBPLB = (hFiredNBPLB->GetBinContent(i)/hTotLB->GetBinContent(i))*100;

                    errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLB->GetBinContent(i));
                    errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLB->GetBinContent(i));
                    errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLB->GetBinContent(i));

                    v_hEffLB_both[period]->Fill(i,effBothLB);
                    v_hEffLB_both[period]->SetBinError(i,errEffBothLB);
                    v_hEffLB_BP[period]->Fill(i,effBPLB);
                    v_hEffLB_BP[period]->SetBinError(i,errEffBPLB);
                    v_hEffLB_NBP[period]->Fill(i,effNBPLB);
                    v_hEffLB_NBP[period]->SetBinError(i,errEffNBPLB);

                    v_hEffLB_both[period]->SetLineColor(period+1);
                    hEffPeriodLB->Add(v_hEffLB_both[period]);
                }
            }

            //fInPeriod->Close();
            //delete fInPeriod;

        }

        new TCanvas();
        hEffPeriodLB->Draw("nostack pfc");*/

        bool uploadToCCDB = false; //Only upload to CCDB if this is true
        if (uploadToCCDB) {
            //Variables used in the loop
            int LB936 = 0; //LB 1-> 936
            int plane = 0; //plane

            std::vector<o2::mid::ChEffCounter> counterVector; //Vector of efficiency counter structs
            struct o2::mid::ChEffCounter entry; //Struct for efficiency counters
            o2::mid::ChamberEfficiency effMap; //Chamber efficiency map

            //vector of struct, push back the struct populated with counts
            for (int ide = 0; ide < o2::mid::detparams::NDetectionElements; ++ide) {
                for (int icol = mapping.getFirstColumn(ide); icol < 7; ++icol) {
                    for (int iline = mapping.getFirstBoardBP(icol, ide); iline <= mapping.getLastBoardBP(icol, ide); ++iline) {                   
                        //Debug printout
                        //cout << "det ID " << ide << " col " << icol << " line " << iline << " LB " << mapping.getBoardId(iline,icol,ide) << endl;
                        //cout << "LB: " << mapping.getBoardId(iline,icol,ide) << "\t unique FEEID:" << o2::mid::detparams::makeUniqueFEEId(ide, icol, iline) << "\t";
                        
                        plane = o2::mid::detparams::getChamber(ide); //Get detection plane

                        LB936 = mapping.getBoardId(iline,icol,ide) + 234*plane; //LB translated to 1->936 from 1->234

                        entry.deId =  uint8_t(ide);
                        entry.columnId = uint8_t(icol);
                        entry.lineId = uint8_t(iline);

                        entry.counts[0] = hFiredBPLB->GetBinContent(LB936); //BP
                        entry.counts[1] = hFiredNBPLB->GetBinContent(LB936); //NBP
                        entry.counts[2] = hFiredBothPlanesLB->GetBinContent(LB936); //Both
                        entry.counts[3] = hTotLB->GetBinContent(LB936); //Total

                        counterVector.push_back(entry);
                        
                        //Debug printout
                        //cout << LB936 << "\t both: " << hFiredBothPlanesLB->GetBinContent(LB936) << "\t tot: " << hTotLB->GetBinContent(LB936) << endl;
                    }
                }
            }

            //Debug printouts to test vector filling
            //cout << "Size of counter vector: " << counterVector.size() << endl;
            //auto& polen = counterVector[21];
            //cout << "deId " << std::to_string(polen.deId) << " column " << std::to_string(polen.columnId) << " line " << std::to_string(polen.lineId) << " counts " << polen.counts[2] << endl;

            //Fill effMap from vectors
            effMap.setFromCounters(counterVector);
            
            //Debug, test of getEfficiency function
            //cout << "Test get efficiency: " << effMap.getEfficiency(68,5,1,o2::mid::ChamberEfficiency::EffType::BothPlanes);
            //det ID 68 col 5 line 1 LB 220        
            //getEfficiency(int deId, int columnId, int lineId, EffType type) const

            //Save in CCDB
            auto data = effMap.getCountersAsVector();
            o2::ccdb::CcdbApi api; //CCDB API
            api.init("http://ccdb-test.cern.ch:8080"); //Open connection to test CCDB
            std::map<std::string, std::string> md; //Metada map
            api.storeAsTFileAny(&counterVector, "MID/Calib/ChamberEfficiency", md, 1, o2::ccdb::CcdbObjectInfo::INFINITE_TIMESTAMP); //Upload in CCDB
        } //end if uploadToCCDB == true

    //} //If there is a loop on multiple runs

}