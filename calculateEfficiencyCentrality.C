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
#include "THnSparse.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include "TKey.h"
#include "TLine.h"

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameters
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter

using namespace std;

bool debug = false;

void calculateEfficiencyCentrality() {

    o2::mid::Mapping mapping;

    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    string period  = "DQ_LHC23PbPb_pass4";

    //string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC23_pass4_skimmed_QC1.root"; //pp
    //string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC23zzh_pass4_test1_QC1_small.root"; //Pb-Pb
    //string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResultsLHC23zzn_apass3_all_I-A11_some.root";
    //string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC23_PbPb_pass3_I-A11.root";
    
    string inFileName = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_"+period+".root";

    TFile *fIn = new TFile(inFileName.c_str(),"READ");
    
    //string outFileName = "efficiency.root";
    //TFile *fOut = new TFile(outFileName.c_str(),"RECREATE");

    int nBinsPlane = 4;
    int nBinsRPC = 72;
    int nBinsBoard = 936;

    string planeName[4] = {"MT11","MT12","MT21","MT22"};

    //Declare histo
    TH1F *hEffLBplanes1D_both[4]; //4 1D histograms for LB (one per plane eff on both planes)
    TH1F *hEffLBplanes1D_BP[4]; //4 1D histograms for LB (one per plane eff on BP)
    TH1F *hEffLBplanes1D_NBP[4]; //4 1D histograms for LB (one per plane eff on NBP)
    TH2F *hEffRPCplanes2D_both[4]; //4 2D histograms for RPC efficiency both planes, one per plane
    TH2F *hEffRPCplanes2D_BP[4]; //4 2D histograms for RPC efficiency BP, one per plane
    TH2F *hEffRPCplanes2D_NBP[4]; //4 2D histograms for RPC efficiency NBP, one per plane

    //Initialize them
    for (int i = 0; i < 4; i++) {
        hEffRPCplanes2D_both[i] = new TH2F(("RPC 2D efficiency both "+planeName[i]).c_str(),("RPC 2D efficiency both "+planeName[i]).c_str(),2,-1,1,9,0.5,9.5);
        hEffRPCplanes2D_BP[i] = new TH2F(("RPC 2D efficiency BP "+planeName[i]).c_str(),("RPC 2D efficiency BP "+planeName[i]).c_str(),2,-1,1,9,0.5,9.5);
        hEffRPCplanes2D_NBP[i] = new TH2F(("RPC 2D efficiency NBP "+planeName[i]).c_str(),("RPC 2D efficiency NBP "+planeName[i]).c_str(),2,-1,1,9,0.5,9.5);

        hEffLBplanes1D_both[i] = new TH1F(("LB efficiency both "+planeName[i]).c_str(),("LB efficiency both "+planeName[i]).c_str(),234,0.5,234.5);
        hEffLBplanes1D_BP[i] = new TH1F(("LB efficiency BP "+planeName[i]).c_str(),("LB efficiency BP "+planeName[i]).c_str(),234,0.5,234.5);
        hEffLBplanes1D_NBP[i] = new TH1F(("LB efficiency NBP "+planeName[i]).c_str(),("LB efficiency NBP "+planeName[i]).c_str(),234,0.5,234.5);
    }

    float effBothPlane, effBPPlane, effNBPPlane;
    float errEffBothPlane, errEffBPPlane, errEffNBPPlane;
    
    float effBothRPC, effBPRPC, effNBPRPC;
    float errEffBothRPC, errEffBPRPC, errEffNBPRPC;

    float effBothLB, effBPLB, effNBPLB;
    float errEffBothLB, errEffBPLB, errEffNBPLB;
 
    //vector<double> min = {0,10,20,30};
    //vector<double> max = {10,20,30,40};
    vector<double> min = {1,1,2,3,4};
    vector<double> max = {9,2,3,4,5};
    vector<int> color = {1,2,3,4,6};

    //Open the directory with the output of the task 
    TDirectoryFile *d = (TDirectoryFile*)fIn->Get("mid-efficiency");
    d->cd();

    THStack *hEffBPPlane_centr = new THStack();

    bool centrality = false;

    //Analyze Planes
    TH1F *hEffPlane_both = new TH1F("effPlane_both","effPlane_both",nBinsPlane,-0.5,3.5);
    TH1F *hEffPlane_BP = new TH1F("effPlane_BP","effPlane_BP",nBinsPlane,-0.5,3.5);
    TH1F *hEffPlane_NBP = new TH1F("effPlane_NBP","effPlane_NBP",nBinsPlane,-0.5,3.5);

    if (centrality) {
        //Get plane THnSparse(s)
        THnSparse *hSparseCentFiredTotPerPlane = (THnSparse*)d->Get("hSparseCentFiredTotperPlane");
        THnSparse *hSparseCentFiredBPPerPlane = (THnSparse*)d->Get("hSparseCentFiredBPperPlane");
        
        for (unsigned int j = 0; j < min.size(); j++) {

            //Analyze Planes
            TH1F *hEffPlane_both = new TH1F("effPlane_both","effPlane_both",nBinsPlane,-0.5,3.5);
            TH1F *hEffPlane_BP = new TH1F("effPlane_BP","effPlane_BP",nBinsPlane,-0.5,3.5);
            TH1F *hEffPlane_NBP = new TH1F("effPlane_NBP","effPlane_NBP",nBinsPlane,-0.5,3.5);

            hSparseCentFiredTotPerPlane->GetAxis(1)->SetRange(min.at(j),max.at(j));
            hSparseCentFiredBPPerPlane->GetAxis(1)->SetRange(min.at(j),max.at(j));

            TH1D* totPlaneCountsProj = hSparseCentFiredTotPerPlane->Projection(0); 
            TH1D* BPPlaneCountsProj = hSparseCentFiredBPPerPlane->Projection(0); 

            for (int i = 1; i <= nBinsPlane; i++) {
                //effBothPlane = (bothPlaneCountsProj->GetBinContent(i)/totPlaneCountsProj->GetBinContent(i))*100;
                effBPPlane = (BPPlaneCountsProj->GetBinContent(i)/totPlaneCountsProj->GetBinContent(i))*100;
                //effNBPPlane = (NBPBPPlaneCountsProj->GetBinContent(i)/totPlaneCountsProj->GetBinContent(i))*100;

                //errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/totPlaneCountsProj->GetBinContent(i));
                errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/totPlaneCountsProj->GetBinContent(i));
                //errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/totPlaneCountsProj->GetBinContent(i));

                //hEffPlane_both->SetBinContent(i,effBothPlane);
                //hEffPlane_both->SetBinError(i,errEffBothPlane);
                hEffPlane_BP->SetBinContent(i,effBPPlane);
                hEffPlane_BP->SetBinError(i,errEffBPPlane);
                hEffPlane_BP->SetLineColor(color.at(j));
                //hEffPlane_NBP->SetBinContent(i,effNBPPlane);
                //hEffPlane_NBP->SetBinError(i,errEffNBPPlane);
            }
            hEffBPPlane_centr->Add(hEffPlane_BP);
        }
        
    }
    
    else {
        TH1F *hFiredBothPlanesPlane = (TH1F*)d->Get("nFiredBothperPlane");
        TH1F *hFiredBPPlane= (TH1F*)d->Get("nFiredBPperPlane");
        TH1F *hFiredNBPPlane = (TH1F*)d->Get("nFiredNBPperPlane");
        TH1F *hTotPlane = (TH1F*)d->Get("nTotperPlane");

        for (int i = 1; i <= nBinsPlane; i++) {
            //effBothPlane = (hFiredBothPlanesPlane->GetBinContent(i)/hTotPlane->GetBinContent(i))*100;
            effBPPlane = (hFiredBPPlane->GetBinContent(i)/hTotPlane->GetBinContent(i))*100;
            //effNBPPlane = (hFiredNBPPlane->GetBinContent(i)/hTotPlane->GetBinContent(i))*100;

            //errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/hTotPlane->GetBinContent(i));
            errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/hTotPlane->GetBinContent(i));
            //errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/hTotPlane->GetBinContent(i));

            //hEffPlane_both->SetBinContent(i,effBothPlane);
            //hEffPlane_both->SetBinError(i,errEffBothPlane);
            hEffPlane_BP->SetBinContent(i,effBPPlane);
            hEffPlane_BP->SetBinError(i,errEffBPPlane);
            //hEffPlane_NBP->SetBinContent(i,effNBPPlane);
            //hEffPlane_NBP->SetBinError(i,errEffNBPPlane);
        }

    }

    //THnSparseD* my_thn = (THnSparseD*) _file0->Get("thn_name"); 
    //root [2] my_thn->GetAxis(6)->SetRange(1,1); 
    //root [3] TH1D* my_thn_projE = my_thn->Projection(0); 
    //root [4] my_thn_projE->Draw();

    //Get Plane counts
    /*TH1F *hFiredBothPlanesPlane = (TH1F*)d->Get("nFiredBothperPlane");
    TH1F *hFiredBPPlane= (TH1F*)d->Get("nFiredBPperPlane");
    TH1F *hFiredNBPPlane = (TH1F*)d->Get("nFiredNBPperPlane");
    TH1F *hTotPlane = (TH1F*)d->Get("nTotperPlane");
    
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
    TH1F *hEffLB_NBP = new TH1F("effBoard_NBP","effBoard_NBP",nBinsBoard,0.5,936.5);*/

    //fIn->Close();

    /*for (int i = 1; i <= nBinsPlane; i++) {
        //cout << i << "\t" << hFiredBothPlanesRPC->GetBinContent(i) << "\t" << hTotRPC->GetBinContent(i) << endl;

        if (hTotRPC->GetBinContent(i) != 0) {
            effBothPlane = (hFiredBothPlanesPlane->GetBinContent(i)/hTotPlane->GetBinContent(i))*100;
            effBPPlane = (hFiredBPPlane->GetBinContent(i)/hTotPlane->GetBinContent(i))*100;
            effNBPPlane = (hFiredNBPPlane->GetBinContent(i)/hTotPlane->GetBinContent(i))*100;

            errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/hTotPlane->GetBinContent(i));
            errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/hTotPlane->GetBinContent(i));
            errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/hTotPlane->GetBinContent(i));

            hEffPlane_both->SetBinContent(i,effBothPlane);
            hEffPlane_both->SetBinError(i,errEffBothPlane);
            hEffPlane_BP->SetBinContent(i,effBPPlane);
            hEffPlane_BP->SetBinError(i,errEffBPPlane);
            hEffPlane_NBP->SetBinContent(i,effNBPPlane);
            hEffPlane_NBP->SetBinError(i,errEffNBPPlane);
        }   
    }

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
            hEffRPCplanes2D_BP[0]->SetBinContent(2,i,effBPRPC);
            hEffRPCplanes2D_NBP[0]->SetBinContent(2,i,effNBPRPC);
        }
        else if (i >= 37 && i <= 45) {
            hEffRPCplanes2D_both[0]->SetBinContent(1,i-36,effBothRPC);
            hEffRPCplanes2D_BP[0]->SetBinContent(1,i-36,effBPRPC);
            hEffRPCplanes2D_NBP[0]->SetBinContent(1,i-36,effNBPRPC);
        }

        //MT12
        else if (i >= 10 && i <= 18) {
            hEffRPCplanes2D_both[1]->SetBinContent(2,i-9,effBothRPC);
            hEffRPCplanes2D_BP[1]->SetBinContent(2,i-9,effBPRPC);
            hEffRPCplanes2D_NBP[1]->SetBinContent(2,i-9,effNBPRPC);
        }
        else if (i >= 46 && i <= 54) {
            hEffRPCplanes2D_both[1]->SetBinContent(1,i-45,effBothRPC);
            hEffRPCplanes2D_BP[1]->SetBinContent(1,i-45,effBPRPC);
            hEffRPCplanes2D_NBP[1]->SetBinContent(1,i-45,effNBPRPC);
        }

        //MT21
        else if (i >= 19 && i <= 27) {
            hEffRPCplanes2D_both[2]->SetBinContent(2,i-18,effBothRPC);
            hEffRPCplanes2D_BP[2]->SetBinContent(2,i-18,effBPRPC);
            hEffRPCplanes2D_NBP[2]->SetBinContent(2,i-18,effNBPRPC);
        }
        else if (i >= 55 && i <= 63) {
            hEffRPCplanes2D_both[2]->SetBinContent(1,i-54,effBothRPC);
            hEffRPCplanes2D_BP[2]->SetBinContent(1,i-54,effBPRPC);
            hEffRPCplanes2D_NBP[2]->SetBinContent(1,i-54,effNBPRPC);
        }

        //MT22
        else if (i >= 28 && i <= 36) {
            hEffRPCplanes2D_both[3]->SetBinContent(2,i-27,effBothRPC);
            hEffRPCplanes2D_BP[3]->SetBinContent(2,i-27,effBPRPC);
            hEffRPCplanes2D_NBP[3]->SetBinContent(2,i-27,effNBPRPC);
        }
        else if (i >= 64 && i <= 72) {
            hEffRPCplanes2D_both[3]->SetBinContent(1,i-63,effBothRPC);
            hEffRPCplanes2D_BP[3]->SetBinContent(1,i-63,effBPRPC);
            hEffRPCplanes2D_NBP[3]->SetBinContent(1,i-63,effNBPRPC);
        }
    }

    for (int i = 1; i <= nBinsBoard; i++) {
        //cout << "LB Both planes " <<  i << "\t" << hFiredBothPlanesLB->GetBinContent(i) << "\t" << hTotLB->GetBinContent(i) << endl;
        
        if (hTotLB->GetBinContent(i) == 0) {
            cout << "LB " <<  i << "\t" << hTotLB->GetBinContent(i) << endl;
        }
            
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
                hEffLBplanes1D_BP[0]->SetBinContent(i,effBPLB);
                hEffLBplanes1D_BP[0]->SetBinError(i,errEffBPLB);
                hEffLBplanes1D_NBP[0]->SetBinContent(i,effNBPLB);
                hEffLBplanes1D_NBP[0]->SetBinError(i,errEffNBPLB);
            }
            else if (i >= 235 && i <= 468) {
                hEffLBplanes1D_both[1]->SetBinContent(i-234,effBothLB);
                hEffLBplanes1D_both[1]->SetBinError(i-234,errEffBothLB);
                hEffLBplanes1D_BP[1]->SetBinContent(i-234,effBPLB);
                hEffLBplanes1D_BP[1]->SetBinError(i-234,errEffBPLB);
                hEffLBplanes1D_NBP[1]->SetBinContent(i-234,effNBPLB);
                hEffLBplanes1D_NBP[1]->SetBinError(i-234,errEffNBPLB);
            }
            else if (i>= 469 && i <= 702) {
                hEffLBplanes1D_both[2]->SetBinContent(i-468,effBothLB);
                hEffLBplanes1D_both[2]->SetBinError(i-468,errEffBothLB);
                hEffLBplanes1D_BP[2]->SetBinContent(i-468,effBPLB);
                hEffLBplanes1D_BP[2]->SetBinError(i-468,errEffBPLB);
                hEffLBplanes1D_NBP[2]->SetBinContent(i-468,effNBPLB);
                hEffLBplanes1D_NBP[2]->SetBinError(i-468,errEffNBPLB);
            }

            else {
                hEffLBplanes1D_both[3]->SetBinContent(i-702,effBothLB);
                hEffLBplanes1D_both[3]->SetBinError(i-702,errEffBothLB);
                hEffLBplanes1D_BP[3]->SetBinContent(i-702,effBPLB);
                hEffLBplanes1D_BP[3]->SetBinError(i-702,errEffBPLB);
                hEffLBplanes1D_NBP[3]->SetBinContent(i-702,effNBPLB);
                hEffLBplanes1D_NBP[3]->SetBinError(i-702,errEffNBPLB);
            }
        }
    }*/

    //Eff per plane
    /*TCanvas *cEffBotPlanesPlane = new TCanvas(); //Both
    cEffBotPlanesPlane->cd();
    hEffPlane_both->SetStats(0);
    hEffPlane_both->GetXaxis()->SetTitle("Plane");
    hEffPlane_both->GetYaxis()->SetTitle("Efficiency [%]");
    hEffPlane_both->GetYaxis()->SetRangeUser(0,100);
    hEffPlane_both->Draw("P");*/

    TCanvas *cEffBPPlane = new TCanvas(); //BP
    cEffBPPlane->cd();
    if (centrality) {
        hEffBPPlane_centr->Draw("nostack pfc");
    }
    else {
        hEffPlane_BP->SetStats(0);
        hEffPlane_BP->GetXaxis()->SetTitle("Plane");
        hEffPlane_BP->GetYaxis()->SetTitle("Efficiency [%]");
        hEffPlane_BP->GetYaxis()->SetRangeUser(0,100);
        hEffPlane_BP->Draw("HISTO");
    }
    

    /*TCanvas *cEffNBPPlane = new TCanvas(); //NBP
    cEffNBPPlane->cd();
    hEffPlane_NBP->SetStats(0);
    hEffPlane_NBP->GetXaxis()->SetTitle("Plane");
    hEffPlane_NBP->GetYaxis()->SetTitle("Efficiency [%]");
    hEffPlane_NBP->GetYaxis()->SetRangeUser(0,100);
    hEffPlane_NBP->Draw("P");*/

    //Eff per RPC
    /*TCanvas *cEffBotPlanesRPC = new TCanvas(); //Both
    cEffBotPlanesRPC->cd();
    hEffRPC_both->SetStats(0);
    hEffRPC_both->GetXaxis()->SetTitle("RPC");
    hEffRPC_both->GetYaxis()->SetTitle("Efficiency [%]");
    hEffRPC_both->GetYaxis()->SetRangeUser(0,100);
    hEffRPC_both->Draw("P");

    TCanvas *cEffBPRPC = new TCanvas(); //BP
    cEffBPRPC->cd();
    hEffRPC_BP->SetStats(0);
    hEffRPC_BP->GetXaxis()->SetTitle("RPC");
    hEffRPC_BP->GetYaxis()->SetTitle("Efficiency [%]");
    hEffRPC_BP->GetYaxis()->SetRangeUser(0,100);
    hEffRPC_BP->Draw("P");

    TCanvas *cEffNBPRPC = new TCanvas(); //NBP
    cEffNBPRPC->cd();
    hEffRPC_NBP->SetStats(0);
    hEffRPC_NBP->GetXaxis()->SetTitle("RPC");
    hEffRPC_NBP->GetYaxis()->SetTitle("Efficiency [%]");
    hEffRPC_NBP->GetYaxis()->SetRangeUser(0,100);
    hEffRPC_NBP->Draw("P");*/

    //Eff per LB (1-936 index)
    /*TCanvas *cEffBotPlanesBoard = new TCanvas(); //Both
    cEffBotPlanesBoard->cd();
    hEffLB_both->SetStats(0);
    hEffLB_both->GetXaxis()->SetTitle("Local Board");
    hEffLB_both->GetYaxis()->SetTitle("Efficiency [%]");
    hEffLB_both->Draw("HISTO");

    TCanvas *cEffBPBoard = new TCanvas(); //BP
    cEffBPBoard->cd();
    hEffLB_BP->SetStats(0);
    hEffLB_BP->GetXaxis()->SetTitle("Local Board");
    hEffLB_BP->GetYaxis()->SetTitle("Efficiency [%]");
    hEffLB_BP->Draw("HISTO");

    TCanvas *cEffNBPBoard = new TCanvas(); //NBP
    cEffNBPBoard->cd();
    hEffLB_NBP->SetStats(0);
    hEffLB_NBP->GetXaxis()->SetTitle("Local Board");
    hEffLB_NBP->GetYaxis()->SetTitle("Efficiency [%]");
    hEffLB_NBP->Draw("HISTO");*/

    //Canvas for LB efficiency per plane
    /*TCanvas *cEffLB_plane_both[4];
    TCanvas *cEffLB_plane_BP[4];
    TCanvas *cEffLB_plane_NBP[4];

    for (int i = 0; i < 4; i++) {
        cEffLB_plane_both[i] = new TCanvas();
        cEffLB_plane_BP[i] = new TCanvas();
        cEffLB_plane_NBP[i] = new TCanvas();
        //both
        cEffLB_plane_both[i]->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        cEffLB_plane_both[i]->SetCanvasSize(1200,1200);
        hEffLBplanes1D_both[i]->SetTitle((planeName[i]+" both").c_str());
        hEffLBplanes1D_both[i]->SetStats(0);
        //hEffLBplanes1D_both[i]->SetLineColor(color[period]);
        //hEffLBplanes1D_both[i]->SetMarkerColor(color[period]);
        hEffLBplanes1D_both[i]->SetMarkerStyle(8);
        hEffLBplanes1D_both[i]->SetMarkerSize(.8);
        hEffLBplanes1D_both[i]->GetXaxis()->SetTitle("Local board");
        hEffLBplanes1D_both[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEffLBplanes1D_both[i]->GetXaxis()->CenterTitle(true);
        hEffLBplanes1D_both[i]->GetYaxis()->CenterTitle(true);
        hEffLBplanes1D_both[i]->GetYaxis()->SetRangeUser(0,105);
        hEffLBplanes1D_both[i]->GetXaxis()->SetTitleFont(62);
        hEffLBplanes1D_both[i]->GetYaxis()->SetTitleFont(62);
        hEffLBplanes1D_both[i]->GetXaxis()->SetLabelFont(62);
        hEffLBplanes1D_both[i]->GetYaxis()->SetLabelFont(62);
        hEffLBplanes1D_both[i]->GetYaxis()->SetTitleOffset(1.1);
        //if (i == 0) {
        //    lTot->AddEntry(hEffLBplanes1D_both[0],"LHC23 pass4 skimmed QC1","p");
        //}
        hEffLBplanes1D_both[i]->Draw("P");
        lTot->Draw("SAME");
        cEffLB_plane_both[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/LB_bothPlanes_total_"+planeName[i]+".pdf").c_str());
        
        //BP
        cEffLB_plane_BP[i]->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        cEffLB_plane_BP[i]->SetCanvasSize(1200,1200);
        hEffLBplanes1D_BP[i]->SetTitle((planeName[i]+" BP").c_str());
        hEffLBplanes1D_BP[i]->SetStats(0);
        //hEffLBplanes1D_BP[i]->SetLineColor(color[period]);
        //hEffLBplanes1D_BP[i]->SetMarkerColor(color[period]);
        hEffLBplanes1D_BP[i]->SetMarkerStyle(8);
        hEffLBplanes1D_BP[i]->SetMarkerSize(.8);
        hEffLBplanes1D_BP[i]->GetXaxis()->SetTitle("Local board");
        hEffLBplanes1D_BP[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEffLBplanes1D_BP[i]->GetXaxis()->CenterTitle(true);
        hEffLBplanes1D_BP[i]->GetYaxis()->CenterTitle(true);
        hEffLBplanes1D_BP[i]->GetYaxis()->SetRangeUser(0,105);
        hEffLBplanes1D_BP[i]->GetXaxis()->SetTitleFont(62);
        hEffLBplanes1D_BP[i]->GetYaxis()->SetTitleFont(62);
        hEffLBplanes1D_BP[i]->GetXaxis()->SetLabelFont(62);
        hEffLBplanes1D_BP[i]->GetYaxis()->SetLabelFont(62);
        hEffLBplanes1D_BP[i]->GetYaxis()->SetTitleOffset(1.1);
        hEffLBplanes1D_BP[i]->Draw("P");
        lTot->Draw("SAME");
        cEffLB_plane_BP[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/LB_BP_total_"+planeName[i]+".pdf").c_str());

        //BP
        cEffLB_plane_NBP[i]->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        cEffLB_plane_NBP[i]->SetCanvasSize(1200,1200);
        hEffLBplanes1D_NBP[i]->SetTitle((planeName[i]+" BP").c_str());
        hEffLBplanes1D_NBP[i]->SetStats(0);
        //hEffLBplanes1D_NBP[i]->SetLineColor(color[period]);
        //hEffLBplanes1D_NBP[i]->SetMarkerColor(color[period]);
        hEffLBplanes1D_NBP[i]->SetMarkerStyle(8);
        hEffLBplanes1D_NBP[i]->SetMarkerSize(.8);
        hEffLBplanes1D_NBP[i]->GetXaxis()->SetTitle("Local board");
        hEffLBplanes1D_NBP[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEffLBplanes1D_NBP[i]->GetXaxis()->CenterTitle(true);
        hEffLBplanes1D_NBP[i]->GetYaxis()->CenterTitle(true);
        hEffLBplanes1D_NBP[i]->GetYaxis()->SetRangeUser(0,105);
        hEffLBplanes1D_NBP[i]->GetXaxis()->SetTitleFont(62);
        hEffLBplanes1D_NBP[i]->GetYaxis()->SetTitleFont(62);
        hEffLBplanes1D_NBP[i]->GetXaxis()->SetLabelFont(62);
        hEffLBplanes1D_NBP[i]->GetYaxis()->SetLabelFont(62);
        hEffLBplanes1D_NBP[i]->GetYaxis()->SetTitleOffset(1.1);
        hEffLBplanes1D_NBP[i]->Draw("P");
        lTot->Draw("SAME");
        cEffLB_plane_NBP[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/LB_NBP_total_"+planeName[i]+".pdf").c_str());   
    }

    //2D eff map per RPC both planes
    TCanvas *cEffRPC_plane_both = new TCanvas();
    cEffRPC_plane_both->SetCanvasSize(1200,1200);
    cEffRPC_plane_both->Divide(2,2);
    
    for (int plane = 0; plane < nBinsPlane; plane++) {
        cEffRPC_plane_both->cd(plane+1);
        hEffRPCplanes2D_both[plane]->SetTitle((planeName[plane]+" both").c_str());
        hEffRPCplanes2D_both[plane]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_both[plane]->SetStats(0);
        hEffRPCplanes2D_both[plane]->GetXaxis()->SetNdivisions(503);
        hEffRPCplanes2D_both[plane]->Draw("COLZ");
    } 
    cEffRPC_plane_both->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/2D_both_"+period+".png").c_str());

    //2D eff map per RPC BP
    TCanvas *cEffRPC_plane_BP = new TCanvas();
    cEffRPC_plane_BP->SetCanvasSize(1200,1200);
    cEffRPC_plane_BP->Divide(2,2);
    
    for (int plane = 0; plane < nBinsPlane; plane++) {
        cEffRPC_plane_BP->cd(plane+1);
        hEffRPCplanes2D_BP[plane]->SetTitle((planeName[plane]+" BP").c_str());
        hEffRPCplanes2D_BP[plane]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_BP[plane]->SetStats(0);
        hEffRPCplanes2D_BP[plane]->GetXaxis()->SetNdivisions(503);
        hEffRPCplanes2D_BP[plane]->Draw("COLZ");
    }
    cEffRPC_plane_BP->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/2D_BP_"+period+".png").c_str());

    //2D eff map per RPC NBP
    TCanvas *cEffRPC_plane_NBP = new TCanvas();
    cEffRPC_plane_NBP->SetCanvasSize(1200,1200);
    cEffRPC_plane_NBP->Divide(2,2);
    
    for (int plane = 0; plane < nBinsPlane; plane++) {
        cEffRPC_plane_NBP->cd(plane+1);
        hEffRPCplanes2D_NBP[plane]->SetTitle((planeName[plane]+" NBP").c_str());
        hEffRPCplanes2D_NBP[plane]->GetZaxis()->SetRangeUser(0,100);
        hEffRPCplanes2D_NBP[plane]->SetStats(0);
        hEffRPCplanes2D_NBP[plane]->GetXaxis()->SetNdivisions(503);
        hEffRPCplanes2D_NBP[plane]->Draw("COLZ");
    }
    cEffRPC_plane_NBP->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/2D_NBP_"+period+".png").c_str());

    //Analyze the 4 periods of LHC23 pass4 skimmed QC1
    string periods[4] = {"za","zj","zs","zt"};

    TH1F *v_hEffLB_both[4];
    TH1F *v_hEffLB_BP[4];
    TH1F *v_hEffLB_NBP[4];

    TH1F *hFiredBothPlanesLBPeriod[4];
    TH1F *hFiredBPLBPeriod[4];
    TH1F *hFiredNBPLBPeriod[4];
    TH1F *hTotLBPeriod[4];

    string fileName[4];

    //Array of TFiles to open period-by-period .root files
    TFile *fInPeriod[4];

    //THStack for all planes in all periods (LB 1->936)
    THStack *hEffPeriodLB = new THStack();
    //Array of THStacks for the four planes (LB 1->234)
    THStack *hEffPeriodLB_planes_both[4];
    THStack *hEffPeriodLB_planes_BP[4];
    THStack *hEffPeriodLB_planes_NBP[4];
    //Matrix of histos: 4 elements because of 4 periods and 4 planes per period
    TH1F *hEffLB_period_planes_both[4][4]; //both
    TH1F *hEffLB_period_planes_BP[4][4]; //BP
    TH1F *hEffLB_period_planes_NBP[4][4]; //NBP

    //Legend of the different periods
    TLegend *lPeriods = new TLegend(0.720,0.195,0.885,0.357,"","rNDC");
    lPeriods->SetFillStyle(0); //Transparent background
    lPeriods->SetTextFont(62); //Bold legend

    //Array of colors for histos in the stack
    int color[4] = {2,3,4,6};

    for (int i = 0; i < 4; i++) {
        hEffPeriodLB_planes_both[i] = new THStack();
        hEffPeriodLB_planes_BP[i] = new THStack();
        hEffPeriodLB_planes_NBP[i] = new THStack();
    }

    for (int period = 0; period < 4; period++) {
        
        for (int j = 0; j < 4; j++) {
            hEffLB_period_planes_both[period][j] = new TH1F(("hEffLB_period_plane_both "+to_string(j)+periods[period]).c_str(),("hEffLB_period_plane_both "+to_string(j)+periods[period]).c_str(),234,-0.5,235.5);
            hEffLB_period_planes_BP[period][j] = new TH1F(("hEffLB_period_plane_BP "+to_string(j)+periods[period]).c_str(),("hEffLB_period_plane_BP "+to_string(j)+periods[period]).c_str(),234,-0.5,235.5);
            hEffLB_period_planes_NBP[period][j] = new TH1F(("hEffLB_period_plane_NBP "+to_string(j)+periods[period]).c_str(),("hEffLB_period_plane_NBP "+to_string(j)+periods[period]).c_str(),234,-0.5,235.5);
        }
    
        fileName[period] = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/AnalysisResults_LHC23"+periods[period]+".root";

        fInPeriod[period] = new TFile(fileName[period].c_str(),"READ");
        fInPeriod[period]->cd();
        
        //cout << fileName[period] << endl;

        TDirectoryFile *dPeriod = (TDirectoryFile*)fInPeriod[period]->Get("mid-efficiency");
        dPeriod->cd();

        hFiredBothPlanesLBPeriod[period] = new TH1F();
        hFiredBPLBPeriod[period] = new TH1F();
        hFiredNBPLBPeriod[period] = new TH1F();
        hTotLBPeriod[period] = new TH1F();

        hFiredBothPlanesLBPeriod[period] = (TH1F*)dPeriod->Get("nFiredBothperBoard");
        hFiredBPLBPeriod[period] = (TH1F*)dPeriod->Get("nFiredBPperBoard");
        hFiredNBPLBPeriod[period] = (TH1F*)dPeriod->Get("nFiredNBPperBoard");
        hTotLBPeriod[period] = (TH1F*)dPeriod->Get("nTotperBoard");

        v_hEffLB_both[period] = new TH1F(("effLB_both"+periods[period]).c_str(),("effLB_both"+periods[period]).c_str(),nBinsBoard,-0.5,935.5);
        v_hEffLB_BP[period] = new TH1F(("effLB_BP"+periods[period]).c_str(),("effLB_BP"+periods[period]).c_str(),nBinsBoard,-0.5,935.5);
        v_hEffLB_NBP[period] = new TH1F(("effLB_NBP"+periods[period]).c_str(),("effLB_NBP"+periods[period]).c_str(),nBinsBoard,-0.5,935.5);

        for (int i = 1; i <= nBinsBoard; i++) {
            if (hTotLBPeriod[period]->GetBinContent(i) != 0) {

                effBothLB = (hFiredBothPlanesLBPeriod[period]->GetBinContent(i)/hTotLBPeriod[period]->GetBinContent(i))*100;
                effBPLB = (hFiredBPLBPeriod[period]->GetBinContent(i)/hTotLBPeriod[period]->GetBinContent(i))*100;
                effNBPLB = (hFiredNBPLBPeriod[period]->GetBinContent(i)/hTotLBPeriod[period]->GetBinContent(i))*100;

                errEffBothLB = TMath::Sqrt(effBothLB*(100-effBothLB)/hTotLBPeriod[period]->GetBinContent(i));
                errEffBPLB = TMath::Sqrt(effBPLB*(100-effBPLB)/hTotLBPeriod[period]->GetBinContent(i));
                errEffNBPLB = TMath::Sqrt(effNBPLB*(100-effNBPLB)/hTotLBPeriod[period]->GetBinContent(i));

                v_hEffLB_both[period]->Fill(i,effBothLB);
                v_hEffLB_both[period]->SetBinError(i,errEffBothLB);
                v_hEffLB_BP[period]->Fill(i,effBPLB);
                v_hEffLB_BP[period]->SetBinError(i,errEffBPLB);
                v_hEffLB_NBP[period]->Fill(i,effNBPLB);
                v_hEffLB_NBP[period]->SetBinError(i,errEffNBPLB);

                //cout <<periods[period] << "\t" << effBothLB << "\t" << errEffBothLB << "\t" << hTotLBPeriod[period]->GetBinContent(i) << endl;

                //Efficiency per plane
                if (i <= 234) {
                    hEffLB_period_planes_both[period][0]->SetBinContent(i,effBothLB);
                    hEffLB_period_planes_both[period][0]->SetBinError(i,errEffBothLB);
                    hEffLB_period_planes_BP[period][0]->SetBinContent(i,effBPLB);
                    hEffLB_period_planes_BP[period][0]->SetBinError(i,errEffBPLB);
                    hEffLB_period_planes_NBP[period][0]->SetBinContent(i,effBPLB);
                    hEffLB_period_planes_NBP[period][0]->SetBinError(i,errEffBPLB);
                }
                else if (i >= 235 && i <= 468) {
                    hEffLB_period_planes_both[period][1]->SetBinContent(i-234,effBothLB);
                    hEffLB_period_planes_both[period][1]->SetBinError(i-234,errEffBothLB);
                    hEffLB_period_planes_BP[period][1]->SetBinContent(i-234,effBPLB);
                    hEffLB_period_planes_BP[period][1]->SetBinError(i-234,errEffBPLB);
                    hEffLB_period_planes_NBP[period][1]->SetBinContent(i-234,effNBPLB);
                    hEffLB_period_planes_NBP[period][1]->SetBinError(i-234,errEffNBPLB);
                }
                else if (i>= 469 && i <= 702) {
                    hEffLB_period_planes_both[period][2]->SetBinContent(i-468,effBothLB);
                    hEffLB_period_planes_both[period][2]->SetBinError(i-468,errEffBothLB);
                    hEffLB_period_planes_BP[period][2]->SetBinContent(i-468,effBPLB);
                    hEffLB_period_planes_BP[period][2]->SetBinError(i-468,errEffBPLB);
                    hEffLB_period_planes_NBP[period][2]->SetBinContent(i-468,effNBPLB);
                    hEffLB_period_planes_NBP[period][2]->SetBinError(i-468,errEffNBPLB);
                }

                else {
                    hEffLB_period_planes_both[period][3]->SetBinContent(i-702,effBothLB);
                    hEffLB_period_planes_both[period][3]->SetBinError(i-702,errEffBothLB);
                    hEffLB_period_planes_BP[period][3]->SetBinContent(i-702,effBPLB);
                    hEffLB_period_planes_BP[period][3]->SetBinError(i-702,errEffBPLB);
                    hEffLB_period_planes_NBP[period][3]->SetBinContent(i-702,effNBPLB);
                    hEffLB_period_planes_NBP[period][3]->SetBinError(i-702,errEffNBPLB);
                }
            }
        }

        //Add histo to THStack and set marker color/size/style
        for (int j = 0; j < 4; j++) {
            //Both planes
            hEffLB_period_planes_both[period][j]->SetLineColor(color[period]);
            hEffLB_period_planes_both[period][j]->SetMarkerColor(color[period]);
            hEffLB_period_planes_both[period][j]->SetMarkerStyle(8);
            hEffLB_period_planes_both[period][j]->SetMarkerSize(.8);
            hEffPeriodLB_planes_both[j]->Add(hEffLB_period_planes_both[period][j]);
            //BP
            hEffLB_period_planes_BP[period][j]->SetLineColor(color[period]);
            hEffLB_period_planes_BP[period][j]->SetMarkerColor(color[period]);
            hEffLB_period_planes_BP[period][j]->SetMarkerStyle(8);
            hEffLB_period_planes_BP[period][j]->SetMarkerSize(.8);
            hEffPeriodLB_planes_BP[j]->Add(hEffLB_period_planes_BP[period][j]);
            //NBP
            hEffLB_period_planes_NBP[period][j]->SetLineColor(color[period]);
            hEffLB_period_planes_NBP[period][j]->SetMarkerColor(color[period]);
            hEffLB_period_planes_NBP[period][j]->SetMarkerStyle(8);
            hEffLB_period_planes_NBP[period][j]->SetMarkerSize(.8);
            hEffPeriodLB_planes_NBP[j]->Add(hEffLB_period_planes_NBP[period][j]);
        }

        //Add Legend entry taking the MT11 historgam (taking any other would work, this is just randomly chosed)
        lPeriods->AddEntry(hEffLB_period_planes_both[period][0],("LHC23_"+periods[period]).c_str(),"p");

        v_hEffLB_both[period]->SetLineColor(period+1);
        hEffPeriodLB->Add(v_hEffLB_both[period]);
        
    }
    
    TCanvas *cTot = new TCanvas();
    cTot->cd();
    hEffPeriodLB->Draw("nostack");

    //4 canvases, one per plane
    TCanvas *cEffLBPeriod_both[4];
    TCanvas *cEffLBPeriod_BP[4];
    TCanvas *cEffLBPeriod_NBP[4];

    for (int i = 0; i < 4; i++) {
        cEffLBPeriod_both[i] =  new TCanvas();
        cEffLBPeriod_BP[i] =  new TCanvas();
        cEffLBPeriod_NBP[i] =  new TCanvas();
        
        //Both
        cEffLBPeriod_both[i]->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        cEffLBPeriod_both[i]->SetCanvasSize(1200,1200);
        hEffPeriodLB_planes_both[i]->SetTitle((planeName[i]+" both").c_str());
        hEffPeriodLB_planes_both[i]->Draw("nostack");
        lPeriods->Draw("SAME");
        hEffPeriodLB_planes_both[i]->GetXaxis()->SetTitle("Local board");
        hEffPeriodLB_planes_both[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEffPeriodLB_planes_both[i]->GetXaxis()->CenterTitle(true);
        hEffPeriodLB_planes_both[i]->GetYaxis()->CenterTitle(true);
        hEffPeriodLB_planes_both[i]->GetXaxis()->SetTitleFont(62);
        hEffPeriodLB_planes_both[i]->GetYaxis()->SetTitleFont(62);
        hEffPeriodLB_planes_both[i]->GetXaxis()->SetLabelFont(62);
        hEffPeriodLB_planes_both[i]->GetYaxis()->SetLabelFont(62);
        hEffPeriodLB_planes_both[i]->GetYaxis()->SetTitleOffset(1.1);
        cEffLBPeriod_both[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/LB_bothPlanes_periods_"+planeName[i]+".pdf").c_str());
    
        //BP
        cEffLBPeriod_BP[i]->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        cEffLBPeriod_BP[i]->SetCanvasSize(1200,1200);
        hEffPeriodLB_planes_BP[i]->SetTitle((planeName[i]+" BP").c_str());
        hEffPeriodLB_planes_BP[i]->Draw("nostack");
        lPeriods->Draw("SAME");
        hEffPeriodLB_planes_BP[i]->GetXaxis()->SetTitle("Local board");
        hEffPeriodLB_planes_BP[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEffPeriodLB_planes_BP[i]->GetXaxis()->CenterTitle(true);
        hEffPeriodLB_planes_BP[i]->GetYaxis()->CenterTitle(true);
        hEffPeriodLB_planes_BP[i]->GetXaxis()->SetTitleFont(62);
        hEffPeriodLB_planes_BP[i]->GetYaxis()->SetTitleFont(62);
        hEffPeriodLB_planes_BP[i]->GetXaxis()->SetLabelFont(62);
        hEffPeriodLB_planes_BP[i]->GetYaxis()->SetLabelFont(62);
        hEffPeriodLB_planes_BP[i]->GetYaxis()->SetTitleOffset(1.1);
        cEffLBPeriod_BP[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/LB_BP_periods_"+planeName[i]+".pdf").c_str());

        //NBP
        cEffLBPeriod_NBP[i]->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        cEffLBPeriod_NBP[i]->SetCanvasSize(1200,1200);
        hEffPeriodLB_planes_NBP[i]->SetTitle((planeName[i]+" NBP").c_str());
        hEffPeriodLB_planes_NBP[i]->Draw("nostack");
        lPeriods->Draw("SAME");
        hEffPeriodLB_planes_NBP[i]->GetXaxis()->SetTitle("Local board");
        hEffPeriodLB_planes_NBP[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEffPeriodLB_planes_NBP[i]->GetXaxis()->CenterTitle(true);
        hEffPeriodLB_planes_NBP[i]->GetYaxis()->CenterTitle(true);
        hEffPeriodLB_planes_NBP[i]->GetXaxis()->SetTitleFont(62);
        hEffPeriodLB_planes_NBP[i]->GetYaxis()->SetTitleFont(62);
        hEffPeriodLB_planes_NBP[i]->GetXaxis()->SetLabelFont(62);
        hEffPeriodLB_planes_NBP[i]->GetYaxis()->SetLabelFont(62);
        hEffPeriodLB_planes_NBP[i]->GetYaxis()->SetTitleOffset(1.1);
        cEffLBPeriod_NBP[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/LB_NBP_periods_"+planeName[i]+".pdf").c_str());
    }*/

}