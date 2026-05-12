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
#include <unordered_map>
#include <chrono>
#include <thread>

#include "MIDEfficiency/Efficiency.h" //MID efficiency
#include "MIDBase/DetectorParameters.h" //Detector parameters
#include "MIDBase/Mapping.h" //MID mapping
#include "DataFormatsMID/Track.h" //MID track from O2
#include "DataFormatsMID/ChEffCounter.h" //Chamber efficiency counter


//# ev vs phi 

using namespace std;

bool debug = false;

void calculateEfficiencyCentrality() {

    //MID RPC mapping
    o2::mid::Mapping mapping;

    //General path to add flexibility to the code + period name
    //string period = "LHC23_pass4_skimmed_QC1"; //pp skimmed QC data of 2023 pass 4
    //string period = "LHC23_PbPb_pass3_I-A11"; //Pb-Pb dataset - one of the two used for the analyses of Nazar
    //string period = "LHC23_PbPb_pass3_fullTPC"; //Pb-Pb dataset - other used for the analyses of Nazar
    //string period = "LHC22o_pass7_minBias";
    //string period = "LHC22_pass7_skimmed";
    //string period = "LHC23_pass4_skimmed";
    string period = "LHC23_PbPb_pass4";
    //string period = "LHC24_pass1_skimmed";
    //string period = "LHC24_ppref_pass1"; 
    //string period = "LHC25ad_pass2"; //pO
    //string period = "LHC25ae_pass2"; //O-O
    //string period = "LHC25af_pass2"; //Ne-Ne

    bool isPbPb;

    if (period == "LHC23_PbPb_pass4" || period == "LHC24_PbPb_pass1" || period == "LHC25ae_pass2" || period == "LHC25af_pass2") {
        isPbPb = true;
    }
    else {
        isPbPb = false;
    }

    //string inFileName = "/media/luca/Elements/MIDefficiency/AnalysisResults_"+period+".root";
    string inFileName = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/AnalysisResults.root";

    cout << "Analyzing file: " << inFileName << " isPbPb: " << isPbPb << endl;

    TFile *fIn = new TFile(inFileName.c_str(),"READ");
    
    //string outFileName = "efficiency.root";
    //TFile *fOut = new TFile(outFileName.c_str(),"RECREATE");

    float maxCentrality, maxEta, maxPhi, maxPt;

    float minCentrality, minEta, minPhi, minPt;

    int binsCentrality = 9;
    int binsPt = 150;
    int binsEta = 100;
    int binsPhi = 100;

    int nBinsPlane = 4;
    int nBinsRPC = 72;
    int nBinsBoard = 936;

    if (isPbPb) {
        maxCentrality = 90.;
        maxEta = 0.;
        maxPhi = 4.;
        maxPt = 20.;

        minCentrality = 0.;
        minEta = -5;
        minPhi = -4;
        minPt = 0.; 
    }

    else {
        maxCentrality = 90.;
        maxEta = 5.;
        maxPhi = 6.28319;
        maxPt = 20.;

        minCentrality = 0.;
        minEta = -5;
        minPhi = -6.28319;
        minPt = 0.;
    }

    float ptStep = (maxPt-minPt)/binsPt;
    float etaStep = (maxEta-minEta)/binsEta;
    float phiStep = (maxPhi-minPhi)/binsPhi;
    float centStep = (maxCentrality-minCentrality)/binsCentrality;

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
 
    //1-2 = 0-10%
    //2-3 = 10-20%
    //3-4 = 20-30% 
    //4-5 = 30-40%
    //5-6 = 40-50%
    //6-7 = 60-70% 
    //7-8 = 70-80%
    //8-9 = 80-90%
    //Map between bin number and corresponding centrality
    //centr[1] -> prints 0
    unordered_map<int, int> centrLow = {{0,0},{1,10},{2,20},{3,30},{4,40},{5,50},{6,60},{7,70},{8,80}};
    unordered_map<int, int> centrHigh = {{0,10},{1,20},{2,30},{3,40},{4,50},{5,60},{6,70},{7,80},{8,90}};
    
    //vector<double> min = {1,1,2,3,8};
    //vector<double> max = {9,2,3,4,9};
    vector<double> min = {1,1,2,3,4};
    vector<double> max = {9,1,2,3,4};
    vector<int> color = {1,2,3,4,6};

    //Open the directory with the output of the task 
    TDirectoryFile *d = (TDirectoryFile*)fIn->Get("mid-efficiency");
    d->cd();

    THStack *hEffBPPlane_centr = new THStack();

    //Analyze Planes
    TH1F *hEffPlane_both = new TH1F("effPlane_both","effPlane_both",nBinsPlane,-0.5,3.5);
    TH1F *hEffPlane_BP = new TH1F("effPlane_BP","effPlane_BP",nBinsPlane,-0.5,3.5);
    TH1F *hEffPlane_NBP = new TH1F("effPlane_NBP","effPlane_NBP",nBinsPlane,-0.5,3.5);

    //Get plane THNsparses once and for all to be used later on
    THnSparse *hSparseCentFiredTotPerPlane = (THnSparse*)d->Get("hSparseCentFiredTotperPlane");
    THnSparse *hSparseCentFiredBothPerPlane = (THnSparse*)d->Get("hSparseCentFiredBothperPlane");
    THnSparse *hSparseCentFiredBPPerPlane = (THnSparse*)d->Get("hSparseCentFiredBPperPlane");
    THnSparse *hSparseCentFiredNBPPerPlane = (THnSparse*)d->Get("hSparseCentFiredNBPperPlane");

    THnSparse *hSparseCentFiredTotPerRPC = (THnSparse*)d->Get("hSparseCentFiredTotperRPC");
    THnSparse *hSparseCentFiredBothPerRPC = (THnSparse*)d->Get("hSparseCentFiredBothperRPC");
    THnSparse *hSparseCentFiredBPPerRPC = (THnSparse*)d->Get("hSparseCentFiredBPperRPC");
    THnSparse *hSparseCentFiredNBPPerRPC = (THnSparse*)d->Get("hSparseCentFiredNBPperRPC");

    //Centrality on x axis, one histo per plane
    TH1F *hEff_bothPlanes_cent[4];
    TH1F *hEff_BPPlanes_cent[4];
    TH1F *hEff_NBPPlanes_cent[4];
    //pt
    TH1F *hEff_bothPlanes_pt[4];
    TH1F *hEff_BPPlanes_pt[4];
    TH1F *hEff_NBPPlanes_pt[4];
    //eta
    TH1F *hEff_bothPlanes_eta[4];
    TH1F *hEff_BPPlanes_eta[4];
    TH1F *hEff_NBPPlanes_eta[4];
    //phi
    TH1F *hEff_bothPlanes_phi[4];
    TH1F *hEff_BPPlanes_phi[4];
    TH1F *hEff_NBPPlanes_phi[4];

    //Centrality on x axis, one histo per RPC
    TH1F *hEff_bothRPC_cent[72];
    TH1F *hEff_BPRPC_cent[72];
    TH1F *hEff_NBPRPC_cent[72];
    //pt
    TH1F *hEff_bothRPC_pt[72];
    TH1F *hEff_BPRPC_pt[72];
    TH1F *hEff_NBPRPC_pt[72];
    //eta
    TH1F *hEff_bothRPC_eta[72];
    TH1F *hEff_BPRPC_eta[72];
    TH1F *hEff_NBPRPC_eta[72];
    //phi
    TH1F *hEff_bothRPC_phi[72];
    TH1F *hEff_BPRPC_phi[72];
    TH1F *hEff_NBPRPC_phi[72];


    for (int i = 0; i < nBinsPlane; i++) {
        hEff_bothPlanes_cent[i] = new TH1F(("hEff_"+planeName[i]+"_bothPlanes_cent").c_str(),("hEff_"+planeName[i]+"_bothPlanes_cent").c_str(),9,0.,90.);
        hEff_BPPlanes_cent[i] = new TH1F(("hEff_"+planeName[i]+"_BP_cent").c_str(),("hEff_"+planeName[i]+"_BP_cent").c_str(),9,0.,90.);
        hEff_NBPPlanes_cent[i] = new TH1F(("hEff_"+planeName[i]+"_NBP_cent").c_str(),("hEff_"+planeName[i]+"_NBP_cent").c_str(),9,0.,90.);
        //pt
        hEff_bothPlanes_pt[i] = new TH1F(("hEff_"+planeName[i]+"_bothPlanes_pt").c_str(),("hEff_"+planeName[i]+"_bothPlanes_pt").c_str(),150,0.,20.);
        hEff_BPPlanes_pt[i] = new TH1F(("hEff_"+planeName[i]+"_BP_pt").c_str(),("hEff_"+planeName[i]+"_BP_pt").c_str(),150,0.,20.);
        hEff_NBPPlanes_pt[i] = new TH1F(("hEff_"+planeName[i]+"_NBP_pt").c_str(),("hEff_"+planeName[i]+"_NBP_pt").c_str(),150,0.,20.);
        //eta
        hEff_bothPlanes_eta[i] = new TH1F(("hEff_"+planeName[i]+"_bothPlanes_eta").c_str(),("hEff_"+planeName[i]+"_bothPlanes_eta").c_str(),binsEta,minEta,maxEta);
        hEff_BPPlanes_eta[i] = new TH1F(("hEff_"+planeName[i]+"_BP_eta").c_str(),("hEff_"+planeName[i]+"_BP_eta").c_str(),binsEta,minEta,maxEta);
        hEff_NBPPlanes_eta[i] = new TH1F(("hEff_"+planeName[i]+"_NBP_eta").c_str(),("hEff_"+planeName[i]+"_NBP_eta").c_str(),binsEta,minEta,maxEta);
        //phi
        hEff_bothPlanes_phi[i] = new TH1F(("hEff_"+planeName[i]+"_bothPlanes_phi").c_str(),("hEff_"+planeName[i]+"_bothPlanes_phi").c_str(),binsPhi,minPhi,maxPhi);
        hEff_BPPlanes_phi[i] = new TH1F(("hEff_"+planeName[i]+"_BP_phi").c_str(),("hEff_"+planeName[i]+"_BP_phi").c_str(),binsPhi,minPhi,maxPhi);
        hEff_NBPPlanes_phi[i] = new TH1F(("hEff_"+planeName[i]+"_NBP_phi").c_str(),("hEff_"+planeName[i]+"_NBP_phi").c_str(),binsPhi,minPhi,maxPhi);
    }

    for (int i = 0; i < nBinsRPC; i++) {
        hEff_bothRPC_cent[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_bothRPC_cent").c_str(),("hEff_RPC_"+to_string(i)+"_bothRPC_cent").c_str(),9,0.,90.);
        hEff_BPRPC_cent[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_BP_RPC_cent").c_str(),("hEff_RPC_"+to_string(i)+"_BP_RPC_cent").c_str(),9,0.,90.);
        hEff_NBPRPC_cent[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_NBP_RPC_cent").c_str(),("hEff_RPC_"+to_string(i)+"_NBP_RPC_cent").c_str(),9,0.,90.);
        //pt
        hEff_bothRPC_pt[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_bothRPC_pt").c_str(),("hEff_RPC_"+to_string(i)+"_bothRPC_pt").c_str(),150,0.,20.);
        hEff_BPRPC_pt[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_BP_RPC_pt").c_str(),("hEff_RPC_"+to_string(i)+"_BP_RPC_pt").c_str(),150,0.,20.);
        hEff_NBPRPC_pt[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_NBP_RPC_pt").c_str(),("hEff_RPC_"+to_string(i)+"_NBP_RPC_pt").c_str(),150,0.,20.);
        //eta
        hEff_bothRPC_eta[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_bothRPC_eta").c_str(),("hEff_RPC_"+to_string(i)+"_bothRPC_eta").c_str(),binsEta,minEta,maxEta);
        hEff_BPRPC_eta[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_BP_RPC_eta").c_str(),("hEff_RPC_"+to_string(i)+"_BP_RPC_eta").c_str(),binsEta,minEta,maxEta);
        hEff_NBPRPC_eta[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_NBP_RPC_eta").c_str(),("hEff_RPC_"+to_string(i)+"_NBP_RPC_eta").c_str(),binsEta,minEta,maxEta);
        //phi
        hEff_bothRPC_phi[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_bothRPC_phi").c_str(),("hEff_RPC_"+to_string(i)+"_bothRPC_phi").c_str(),binsPhi,minPhi,maxPhi);
        hEff_BPRPC_phi[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_BP_RPC_phi").c_str(),("hEff_RPC_"+to_string(i)+"_BP_RPC_phi").c_str(),binsPhi,minPhi,maxPhi);
        hEff_NBPRPC_phi[i] = new TH1F(("hEff_RPC_"+to_string(i)+"_NBP_RPC_phi").c_str(),("hEff_RPC_"+to_string(i)+"_NBP_RPC_phi").c_str(),binsPhi,minPhi,maxPhi);
    }

    TLegend *lCentrality = new TLegend(0.129,0.626,0.206,0.881,"","rNDC");
    lCentrality->SetBorderSize(0);	
    lCentrality->SetFillStyle(0);	
    lCentrality->SetTextFont(62);	

    //Perform a centrality cut on the data
    //Only true for Pb-Pb or O-O   
    bool centralityCut = false;
    
    if (!isPbPb) {
        centralityCut = false;
    }

    if(centralityCut) {
        hSparseCentFiredTotPerPlane->GetAxis(1)->SetRange(50/centStep,90/centStep);
        hSparseCentFiredBothPerPlane->GetAxis(1)->SetRange(50/centStep,90/centStep);
        hSparseCentFiredBPPerPlane->GetAxis(1)->SetRange(50/centStep,90/centStep);
        hSparseCentFiredNBPPerPlane->GetAxis(1)->SetRange(50/centStep,90/centStep);

        hSparseCentFiredTotPerRPC->GetAxis(1)->SetRange(50/centStep,90/centStep);
        hSparseCentFiredBothPerRPC->GetAxis(1)->SetRange(50/centStep,90/centStep);
        hSparseCentFiredBPPerRPC->GetAxis(1)->SetRange(50/centStep,90/centStep);
        hSparseCentFiredNBPPerRPC->GetAxis(1)->SetRange(50/centStep,90/centStep);
    }

    //Pt cut, can be used on all datasets (pp, O-O, Pb-Pb, p-O, Ne-Ne)
    bool ptCut = true; //If true, a pt cut is also applied in the computation of the efficiency vs eta
    //Add a cut in pt in the efficiency computation

    //Change the value here with the desired pt cut x/ptStep where x is the desired pt value and the step is calculated earlier
    int ptMinCut = 1.8/ptStep;

    if (ptCut) {

        int ptAxis;

        //If the data is pp or Pb-Pb (O-O, Ne-Ne) the axes on the THNsparse are different
        if (isPbPb) {
            ptAxis = 2;
        }
        else {
            ptAxis = 1;
        }
        
        cout << "Cutting on pt axis (2 for PbPb, 1 for pp): " << ptAxis << " at pt = " << ptMinCut*ptStep << " GeV/c" << endl;
        
        hSparseCentFiredTotPerPlane->GetAxis(ptAxis)->SetRange(ptMinCut,150);
        hSparseCentFiredBothPerPlane->GetAxis(ptAxis)->SetRange(ptMinCut,150);
        hSparseCentFiredBPPerPlane->GetAxis(ptAxis)->SetRange(ptMinCut,150);
        hSparseCentFiredNBPPerPlane->GetAxis(ptAxis)->SetRange(ptMinCut,150);

        hSparseCentFiredTotPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
        hSparseCentFiredBothPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
        hSparseCentFiredBPPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
        hSparseCentFiredNBPPerRPC->GetAxis(ptAxis)->SetRange(ptMinCut,150);
    }

    //Plots vs centrality/pt/eta/phi
    //Depending on the cuts selected earlier the plots might be cut or not

    bool centrality = true;

    if (!isPbPb) {
        centrality = false;
    }

    if (centrality) {  

        TH2D* totPlaneCountsProj = hSparseCentFiredTotPerPlane->Projection(0,1); 
        TH2D* BothPlaneCountsProj = hSparseCentFiredBothPerPlane->Projection(0,1);
        TH2D* BPPlaneCountsProj = hSparseCentFiredBPPerPlane->Projection(0,1);
        TH2D* NBPPlaneCountsProj = hSparseCentFiredNBPPerPlane->Projection(0,1);
        
        TH2D* totRPCCountsProj = hSparseCentFiredTotPerRPC->Projection(0,1); 
        TH2D* BothRPCCountsProj = hSparseCentFiredBothPerRPC->Projection(0,1);
        TH2D* BPRPCCountsProj = hSparseCentFiredBPPerRPC->Projection(0,1);
        TH2D* NBPRPCCountsProj = hSparseCentFiredNBPPerRPC->Projection(0,1);

        for (int centBin = 1; centBin <= totPlaneCountsProj->GetNbinsX(); centBin++) {
 
            for (int i = 1; i <= nBinsPlane; i++) {
                effBothPlane = (BothPlaneCountsProj->GetBinContent(centBin,i)/totPlaneCountsProj->GetBinContent(centBin,i))*100;
                effBPPlane = (BPPlaneCountsProj->GetBinContent(centBin,i)/totPlaneCountsProj->GetBinContent(centBin,i))*100;
                effNBPPlane = (NBPPlaneCountsProj->GetBinContent(centBin,i)/totPlaneCountsProj->GetBinContent(centBin,i))*100;

                errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/totPlaneCountsProj->GetBinContent(centBin,i));
                errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/totPlaneCountsProj->GetBinContent(centBin,i));
                errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/totPlaneCountsProj->GetBinContent(centBin,i));

                cout << effBothPlane << "\t +- \t" << errEffBothPlane << endl;
                cout << effBPPlane << "\t +- \t" << errEffBPPlane << endl;
                cout << effNBPPlane << "\t +- \t" << errEffNBPPlane << endl;

                hEff_bothPlanes_cent[i-1]->SetBinContent(centBin,effBothPlane);
                hEff_bothPlanes_cent[i-1]->SetBinError(centBin,errEffBothPlane);
                hEff_BPPlanes_cent[i-1]->SetBinContent(centBin,effBPPlane);
                hEff_BPPlanes_cent[i-1]->SetBinError(centBin,errEffBPPlane);
                hEff_NBPPlanes_cent[i-1]->SetBinContent(centBin,effNBPPlane);
                hEff_NBPPlanes_cent[i-1]->SetBinError(centBin,errEffNBPPlane); 
            } 
            for (int i = 1; i <= nBinsRPC; i++) {
                effBothRPC = (BothRPCCountsProj->GetBinContent(centBin,i)/totRPCCountsProj->GetBinContent(centBin,i))*100;
                effBPRPC = (BPRPCCountsProj->GetBinContent(centBin,i)/totRPCCountsProj->GetBinContent(centBin,i))*100;
                effNBPRPC = (NBPRPCCountsProj->GetBinContent(centBin,i)/totRPCCountsProj->GetBinContent(centBin,i))*100;

                errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/totRPCCountsProj->GetBinContent(centBin,i));
                errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/totRPCCountsProj->GetBinContent(centBin,i));
                errEffNBPRPC = TMath::Sqrt(effNBPRPC*(100-effNBPRPC)/totRPCCountsProj->GetBinContent(centBin,i));

                cout << effBothRPC << "\t +- \t" << errEffBothRPC << endl;
                cout << effBPRPC << "\t +- \t" << errEffBPRPC << endl;
                cout << effNBPRPC << "\t +- \t" << errEffNBPRPC << endl;

                hEff_bothRPC_cent[i-1]->SetBinContent(centBin,effBothRPC);
                hEff_bothRPC_cent[i-1]->SetBinError(centBin,errEffBothRPC);
                hEff_BPRPC_cent[i-1]->SetBinContent(centBin,effBPRPC);
                hEff_BPRPC_cent[i-1]->SetBinError(centBin,errEffBPRPC);
                hEff_NBPRPC_cent[i-1]->SetBinContent(centBin,effNBPRPC);
                hEff_NBPRPC_cent[i-1]->SetBinError(centBin,errEffNBPRPC); 
            }   
        }   
        delete totPlaneCountsProj;
        delete BothPlaneCountsProj;
        delete BPPlaneCountsProj;
        delete NBPPlaneCountsProj;
        //-----//
        delete totRPCCountsProj;
        delete BothRPCCountsProj;
        delete BPRPCCountsProj;
        delete NBPRPCCountsProj;            
    } //End of if(centrality)

    int ptRange = 1;
    //float ptMinCut = (int)2/ptStep;
    cout << "Pt step (GeV/c): " << ptStep << " phi step (rad): " << phiStep << " eta step (a.u.): " << etaStep << " cent step (%): " << centStep << endl;
    cout << "ptMinCut " << ptMinCut << " maxPt/ptStep " << maxPt/ptStep << endl;
    cout << "min eta: " << TMath::Abs(-2/etaStep) << " max eta: " << TMath::Abs(-4/etaStep) << endl;
    cout << "cent 50%: " << 50/centStep << " cent 90%: " << 90/centStep << endl;
    cout << "min phi: " << (int)(3.3/phiStep) << "max phi: " << (int)(6.6/phiStep) <<endl;

    //Analyze eff vs pt
    bool pt = true;
    if (pt) {
        
        int ptAxis;

        if (isPbPb) {
            ptAxis = 2;
        }
        else {
            ptAxis = 1;
        }
        
        cout << "Analyzing pt with axis (2 for PbPb, 1 for pp): " << ptAxis << endl;

        TH2D* totPlaneCountsProj = hSparseCentFiredTotPerPlane->Projection(0,ptAxis); 
        TH2D* BothPlaneCountsProj = hSparseCentFiredBothPerPlane->Projection(0,ptAxis);
        TH2D* BPPlaneCountsProj = hSparseCentFiredBPPerPlane->Projection(0,ptAxis);
        TH2D* NBPPlaneCountsProj = hSparseCentFiredNBPPerPlane->Projection(0,ptAxis);
            
        TH2D* totRPCCountsProj = hSparseCentFiredTotPerRPC->Projection(0,ptAxis); 
        TH2D* BothRPCCountsProj = hSparseCentFiredBothPerRPC->Projection(0,ptAxis);
        TH2D* BPRPCCountsProj = hSparseCentFiredBPPerRPC->Projection(0,ptAxis);
        TH2D* NBPRPCCountsProj = hSparseCentFiredNBPPerRPC->Projection(0,ptAxis);
    
        cout << "Conteggi pt: " << totPlaneCountsProj->GetNbinsX() << endl << endl;
       
        for (int ptBin = 1; ptBin <= totPlaneCountsProj->GetNbinsX(); ptBin++) {
            cout << ptBin << endl;

            for (int i = 1; i <= nBinsPlane; i++) {

                effBothPlane = (BothPlaneCountsProj->GetBinContent(ptBin,i)/totPlaneCountsProj->GetBinContent(ptBin,i))*100;
                effBPPlane = (BPPlaneCountsProj->GetBinContent(ptBin,i)/totPlaneCountsProj->GetBinContent(ptBin,i))*100;
                effNBPPlane = (NBPPlaneCountsProj->GetBinContent(ptBin,i)/totPlaneCountsProj->GetBinContent(ptBin,i))*100;

                if (std::isnan(effBothPlane)) {
                    continue;
                }

                errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/totPlaneCountsProj->GetBinContent(ptBin,i));
                errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/totPlaneCountsProj->GetBinContent(ptBin,i));
                errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/totPlaneCountsProj->GetBinContent(ptBin,i));

                cout << "pT analysis, plane: " << i-1 << endl;
                //cout << effBothPlane << "\t +- \t" << errEffBothPlane << endl;
                //cout << effBPPlane << "\t +- \t" << errEffBPPlane << endl;
                //cout << effNBPPlane << "\t +- \t" << errEffNBPPlane << endl;

                if (ptCut) {
                    hEff_bothPlanes_pt[i-1]->SetBinContent(ptBin+14,effBothPlane);
                    hEff_bothPlanes_pt[i-1]->SetBinError(ptBin+14,errEffBothPlane);
                    hEff_BPPlanes_pt[i-1]->SetBinContent(ptBin+14,effBPPlane);
                    hEff_BPPlanes_pt[i-1]->SetBinError(ptBin+14,errEffBPPlane);
                    hEff_NBPPlanes_pt[i-1]->SetBinContent(ptBin+14,effNBPPlane);
                    hEff_NBPPlanes_pt[i-1]->SetBinError(ptBin+14,errEffNBPPlane);
                }
                else {
                    hEff_bothPlanes_pt[i-1]->SetBinContent(ptBin,effBothPlane);
                    hEff_bothPlanes_pt[i-1]->SetBinError(ptBin,errEffBothPlane);
                    hEff_BPPlanes_pt[i-1]->SetBinContent(ptBin,effBPPlane);
                    hEff_BPPlanes_pt[i-1]->SetBinError(ptBin,errEffBPPlane);
                    hEff_NBPPlanes_pt[i-1]->SetBinContent(ptBin,effNBPPlane);
                    hEff_NBPPlanes_pt[i-1]->SetBinError(ptBin,errEffNBPPlane);
                }
                 
            }

            for (int i = 1; i <= nBinsRPC; i++) {

                effBothRPC = (BothRPCCountsProj->GetBinContent(ptBin,i)/totRPCCountsProj->GetBinContent(ptBin,i))*100;
                effBPRPC = (BPRPCCountsProj->GetBinContent(ptBin,i)/totRPCCountsProj->GetBinContent(ptBin,i))*100;
                effNBPRPC = (NBPRPCCountsProj->GetBinContent(ptBin,i)/totRPCCountsProj->GetBinContent(ptBin,i))*100;

                if (std::isnan(effBothRPC)) {
                    continue;
                }

                errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/totRPCCountsProj->GetBinContent(ptBin,i));
                errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/totRPCCountsProj->GetBinContent(ptBin,i));
                errEffNBPRPC = TMath::Sqrt(effNBPRPC*(100-effNBPRPC)/totRPCCountsProj->GetBinContent(ptBin,i));

                cout << "pT analysis, RPC: " << i-1 << endl;
                //cout << effBothRPC << "\t +- \t" << errEffBothRPC << endl;
                //cout << effBPRPC << "\t +- \t" << errEffBPRPC << endl;
                //cout << effNBPRPC << "\t +- \t" << errEffNBPRPC << endl;

                hEff_bothRPC_pt[i-1]->SetBinContent(ptBin+14,effBothRPC);
                hEff_bothRPC_pt[i-1]->SetBinError(ptBin+14,errEffBothRPC);
                hEff_BPRPC_pt[i-1]->SetBinContent(ptBin+14,effBPRPC);
                hEff_BPRPC_pt[i-1]->SetBinError(ptBin+14,errEffBPRPC);
                hEff_NBPRPC_pt[i-1]->SetBinContent(ptBin+14,effNBPRPC);
                hEff_NBPRPC_pt[i-1]->SetBinError(ptBin+14,errEffNBPRPC); 

            }
        }
        delete totPlaneCountsProj;
        delete BothPlaneCountsProj;
        delete BPPlaneCountsProj;
        delete NBPPlaneCountsProj;
        //-----//
        delete totRPCCountsProj;
        delete BothRPCCountsProj;
        delete BPRPCCountsProj;
        delete NBPRPCCountsProj;
    } //end of if(pt)

    int etaRange = 1;

    bool eta = false;
    if (eta) {

        int etaAxis;

        if (isPbPb) {
            etaAxis = 3;
        }
        else {
            etaAxis = 2;
        }
        
        cout << "Analyzing eta with axis (2 for PbPb, 1 for pp): " << etaAxis << endl;

        TH2D* totPlaneCountsProj = hSparseCentFiredTotPerPlane->Projection(0,etaAxis); 
        TH2D* BothPlaneCountsProj = hSparseCentFiredBothPerPlane->Projection(0,etaAxis);
        TH2D* BPPlaneCountsProj = hSparseCentFiredBPPerPlane->Projection(0,etaAxis);
        TH2D* NBPPlaneCountsProj = hSparseCentFiredNBPPerPlane->Projection(0,etaAxis);
        
        TH2D* totRPCCountsProj = hSparseCentFiredTotPerRPC->Projection(0,etaAxis); 
        TH2D* BothRPCCountsProj = hSparseCentFiredBothPerRPC->Projection(0,etaAxis);
        TH2D* BPRPCCountsProj = hSparseCentFiredBPPerRPC->Projection(0,etaAxis);
        TH2D* NBPRPCCountsProj = hSparseCentFiredNBPPerRPC->Projection(0,etaAxis);

        //Loop through the eta bins
        for (int etaBin = 1; etaBin <= totPlaneCountsProj->GetNbinsX(); etaBin++) {          
            cout << "eta analysis, bin: " << etaBin << endl;

            for (int i = 1; i <= nBinsPlane; i++) {
                effBothPlane = (BothPlaneCountsProj->GetBinContent(etaBin,i)/totPlaneCountsProj->GetBinContent(etaBin,i))*100;
                effBPPlane = (BPPlaneCountsProj->GetBinContent(etaBin,i)/totPlaneCountsProj->GetBinContent(etaBin,i))*100;
                effNBPPlane = (NBPPlaneCountsProj->GetBinContent(etaBin,i)/totPlaneCountsProj->GetBinContent(etaBin,i))*100;

                if (std::isnan(effBothPlane)) {
                    continue;
                }

                errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/totPlaneCountsProj->GetBinContent(etaBin,i));
                errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/totPlaneCountsProj->GetBinContent(etaBin,i));
                errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/totPlaneCountsProj->GetBinContent(etaBin,i));

                cout << effBothPlane << "\t +- \t" << errEffBothPlane << endl;
                cout << effBPPlane << "\t +- \t" << errEffBPPlane << endl;
                cout << effNBPPlane << "\t +- \t" << errEffNBPPlane << endl;

                hEff_bothPlanes_eta[i-1]->SetBinContent(etaBin,effBothPlane);
                hEff_bothPlanes_eta[i-1]->SetBinError(etaBin,errEffBothPlane);
                hEff_BPPlanes_eta[i-1]->SetBinContent(etaBin,effBPPlane);
                hEff_BPPlanes_eta[i-1]->SetBinError(etaBin,errEffBPPlane);
                hEff_NBPPlanes_eta[i-1]->SetBinContent(etaBin,effNBPPlane);
                hEff_NBPPlanes_eta[i-1]->SetBinError(etaBin,errEffNBPPlane); 
            }

            for (int i = 1; i <= nBinsRPC; i++) {

                effBothRPC = (BothRPCCountsProj->GetBinContent(etaBin,i)/totRPCCountsProj->GetBinContent(etaBin,i))*100;
                effBPRPC = (BPRPCCountsProj->GetBinContent(etaBin,i)/totRPCCountsProj->GetBinContent(etaBin,i))*100;
                effNBPRPC = (NBPRPCCountsProj->GetBinContent(etaBin,i)/totRPCCountsProj->GetBinContent(etaBin,i))*100;

                if (std::isnan(effBothRPC)) {
                    continue;
                }

                errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/totRPCCountsProj->GetBinContent(etaBin,i));
                errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/totRPCCountsProj->GetBinContent(etaBin,i));
                errEffNBPRPC = TMath::Sqrt(effNBPRPC*(100-effNBPRPC)/totRPCCountsProj->GetBinContent(etaBin,i));

                cout << "eta analysis, RPC" << endl;
                cout << effBothRPC << "\t +- \t" << errEffBothRPC << endl;
                cout << effBPRPC << "\t +- \t" << errEffBPRPC << endl;
                cout << effNBPRPC << "\t +- \t" << errEffNBPRPC << endl;

                hEff_bothRPC_eta[i-1]->SetBinContent(etaBin,effBothRPC);
                hEff_bothRPC_eta[i-1]->SetBinError(etaBin,errEffBothRPC);
                hEff_BPRPC_eta[i-1]->SetBinContent(etaBin,effBPRPC);
                hEff_BPRPC_eta[i-1]->SetBinError(etaBin,errEffBPRPC);
                hEff_NBPRPC_eta[i-1]->SetBinContent(etaBin,effNBPRPC);
                hEff_NBPRPC_eta[i-1]->SetBinError(etaBin,errEffNBPRPC); 

            }
        }
        delete totPlaneCountsProj;
        delete BothPlaneCountsProj;
        delete BPPlaneCountsProj;
        delete NBPPlaneCountsProj;
        //-----//
        delete totRPCCountsProj;
        delete BothRPCCountsProj;
        delete BPRPCCountsProj;
        delete NBPRPCCountsProj;
    } //end of if(eta)

    int phiRange = 1;

    bool phi = false;
    if (phi) {

        int phiAxis;

        if (isPbPb) {
            phiAxis = 4;
        }
        else {
            phiAxis = 3;
        }
        
        cout << "Analyzing phi with axis (2 for PbPb, 1 for pp): " << phiAxis << endl;

        TH2D* totPlaneCountsProj = hSparseCentFiredTotPerPlane->Projection(0,phiAxis); 
        TH2D* BothPlaneCountsProj = hSparseCentFiredBothPerPlane->Projection(0,phiAxis);
        TH2D* BPPlaneCountsProj = hSparseCentFiredBPPerPlane->Projection(0,phiAxis);
        TH2D* NBPPlaneCountsProj = hSparseCentFiredNBPPerPlane->Projection(0,phiAxis);
        
        TH2D* totRPCCountsProj = hSparseCentFiredTotPerRPC->Projection(0,phiAxis); 
        TH2D* BothRPCCountsProj = hSparseCentFiredBothPerRPC->Projection(0,phiAxis);
        TH2D* BPRPCCountsProj = hSparseCentFiredBPPerRPC->Projection(0,phiAxis);
        TH2D* NBPRPCCountsProj = hSparseCentFiredNBPPerRPC->Projection(0,phiAxis);
    
        //Loop through all the phi bins
        for (int phiBin = 1; phiBin <= totPlaneCountsProj->GetNbinsX(); phiBin++) {
            cout << "Analyzing phi, bin: " << phiBin << endl;

            for (int i = 1; i <= nBinsPlane; i++) {
                effBothPlane = (BothPlaneCountsProj->GetBinContent(phiBin,i)/totPlaneCountsProj->GetBinContent(phiBin,i))*100;
                effBPPlane = (BPPlaneCountsProj->GetBinContent(phiBin,i)/totPlaneCountsProj->GetBinContent(phiBin,i))*100;
                effNBPPlane = (NBPPlaneCountsProj->GetBinContent(phiBin,i)/totPlaneCountsProj->GetBinContent(phiBin,i))*100;

                if (std::isnan(effBothPlane)) {
                    continue;
                }

                errEffBothPlane = TMath::Sqrt(effBothPlane*(100-effBothPlane)/totPlaneCountsProj->GetBinContent(phiBin,i));
                errEffBPPlane = TMath::Sqrt(effBPPlane*(100-effBPPlane)/totPlaneCountsProj->GetBinContent(phiBin,i));
                errEffNBPPlane = TMath::Sqrt(effNBPPlane*(100-effNBPPlane)/totPlaneCountsProj->GetBinContent(phiBin,i));

                cout << effBothPlane << "\t +- \t" << errEffBothPlane << endl;
                cout << effBPPlane << "\t +- \t" << errEffBPPlane << endl;
                cout << effNBPPlane << "\t +- \t" << errEffNBPPlane << endl;

                hEff_bothPlanes_phi[i-1]->SetBinContent(phiBin,effBothPlane);
                hEff_bothPlanes_phi[i-1]->SetBinError(phiBin,errEffBothPlane);
                hEff_BPPlanes_phi[i-1]->SetBinContent(phiBin,effBPPlane);
                hEff_BPPlanes_phi[i-1]->SetBinError(phiBin,errEffBPPlane);
                hEff_NBPPlanes_phi[i-1]->SetBinContent(phiBin,effNBPPlane);
                hEff_NBPPlanes_phi[i-1]->SetBinError(phiBin,errEffNBPPlane); 
            }
            for (int i = 1; i <= nBinsRPC; i++) {

                effBothRPC = (BothRPCCountsProj->GetBinContent(phiBin,i)/totRPCCountsProj->GetBinContent(phiBin,i))*100;
                effBPRPC = (BPRPCCountsProj->GetBinContent(phiBin,i)/totRPCCountsProj->GetBinContent(phiBin,i))*100;
                effNBPRPC = (NBPRPCCountsProj->GetBinContent(phiBin,i)/totRPCCountsProj->GetBinContent(phiBin,i))*100;

                if (std::isnan(effBothRPC)) {
                    continue;
                }

                errEffBothRPC = TMath::Sqrt(effBothRPC*(100-effBothRPC)/totRPCCountsProj->GetBinContent(phiBin,i));
                errEffBPRPC = TMath::Sqrt(effBPRPC*(100-effBPRPC)/totRPCCountsProj->GetBinContent(phiBin,i));
                errEffNBPRPC = TMath::Sqrt(effNBPRPC*(100-effNBPRPC)/totRPCCountsProj->GetBinContent(phiBin,i));

                cout << "phi analysis, RPC" << endl;
                cout << effBothRPC << "\t +- \t" << errEffBothRPC << endl;
                cout << effBPRPC << "\t +- \t" << errEffBPRPC << endl;
                cout << effNBPRPC << "\t +- \t" << errEffNBPRPC << endl;

                hEff_bothRPC_phi[i-1]->SetBinContent(phiBin,effBothRPC);
                hEff_bothRPC_phi[i-1]->SetBinError(phiBin,errEffBothRPC);
                hEff_BPRPC_phi[i-1]->SetBinContent(phiBin,effBPRPC);
                hEff_BPRPC_phi[i-1]->SetBinError(phiBin,errEffBPRPC);
                hEff_NBPRPC_phi[i-1]->SetBinContent(phiBin,effNBPRPC);
                hEff_NBPRPC_phi[i-1]->SetBinError(phiBin,errEffNBPRPC); 

            }
        }
        delete totPlaneCountsProj;
        delete BothPlaneCountsProj;
        delete BPPlaneCountsProj;
        delete NBPPlaneCountsProj;
        //-----//
        delete totRPCCountsProj;
        delete BothRPCCountsProj;
        delete BPRPCCountsProj;
        delete NBPRPCCountsProj;
    } //end of if(phi)

    TCanvas *cEffBPPlane = new TCanvas(); //BP
    cEffBPPlane->cd();
    if (centrality) {
        hEffBPPlane_centr->Draw("nostack pfc");
        lCentrality->Draw("SAME");
    }
    else {
        hEffPlane_BP->SetStats(0);
        hEffPlane_BP->GetXaxis()->SetTitle("Plane");
        hEffPlane_BP->GetYaxis()->SetTitle("Efficiency [%]");
        hEffPlane_BP->GetYaxis()->SetRangeUser(0,100);
        hEffPlane_BP->Draw("HISTO");
    }

    //Eff per plane at different centralities
    TCanvas *cEffBothPlane_cent[4];
    TCanvas *cEffBPPlane_cent[4];
    TCanvas *cEffNBPPlane_cent[4];

    TCanvas *cEffBothPlane_pt[4];
    TCanvas *cEffBPPlane_pt[4];
    TCanvas *cEffNBPPlane_pt[4];

    TCanvas *cEffBothPlane_eta[4];
    TCanvas *cEffBPPlane_eta[4];
    TCanvas *cEffNBPPlane_eta[4];

    TCanvas *cEffBothPlane_phi[4];
    TCanvas *cEffBPPlane_phi[4];
    TCanvas *cEffNBPPlane_phi[4];

    //Eff per RPC at different centralities
    TCanvas *cEffBothRPC_cent[72];
    TCanvas *cEffBPRPC_cent[72];
    TCanvas *cEffNBPRPC_cent[72];

    TCanvas *cEffBothRPC_pt[72];
    TCanvas *cEffBPRPC_pt[72];
    TCanvas *cEffNBPRPC_pt[72];

    TCanvas *cEffBothRPC_eta[72];
    TCanvas *cEffBPRPC_eta[72];
    TCanvas *cEffNBPRPC_eta[72];

    TCanvas *cEffBothRPC_phi[72];
    TCanvas *cEffBPRPC_phi[72];
    TCanvas *cEffNBPRPC_phi[72];

    //Per plane
    //Eff vs centrality
    for (int i = 0; i < nBinsPlane; i++) {
        cEffBothPlane_cent[i] = new TCanvas();
        cEffBPPlane_cent[i] = new TCanvas();
        cEffNBPPlane_cent[i] = new TCanvas();

        cEffBothPlane_cent[i]->cd();
        cEffBothPlane_cent[i]->SetGridx();
        cEffBothPlane_cent[i]->SetGridy();
        hEff_bothPlanes_cent[i]->SetStats(0);
        hEff_bothPlanes_cent[i]->GetXaxis()->SetTitle("FT0C [%]");
        hEff_bothPlanes_cent[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothPlanes_cent[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothPlanes_cent[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothPlanes_cent[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothPlanes_cent[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothPlanes_cent[i]->GetXaxis()->CenterTitle(true);
        hEff_bothPlanes_cent[i]->GetYaxis()->CenterTitle(true);
        hEff_bothPlanes_cent[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothPlanes_cent[i]->SetLineColor(kBlack);
        hEff_bothPlanes_cent[i]->SetLineWidth(2);
        hEff_bothPlanes_cent[i]->Draw("E");

        cEffBPPlane_cent[i]->cd();
        cEffBPPlane_cent[i]->SetGridx();
        cEffBPPlane_cent[i]->SetGridy();
        hEff_BPPlanes_cent[i]->SetStats(0);
        hEff_BPPlanes_cent[i]->GetXaxis()->SetTitle("FT0C [%]");
        hEff_BPPlanes_cent[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPPlanes_cent[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPPlanes_cent[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPPlanes_cent[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPPlanes_cent[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPPlanes_cent[i]->GetXaxis()->CenterTitle(true);
        hEff_BPPlanes_cent[i]->GetYaxis()->CenterTitle(true);
        hEff_BPPlanes_cent[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPPlanes_cent[i]->SetLineColor(kRed);
        hEff_BPPlanes_cent[i]->SetLineWidth(2);
        hEff_BPPlanes_cent[i]->Draw("E");

        cEffNBPPlane_cent[i]->cd();
        cEffNBPPlane_cent[i]->SetGridx();
        cEffNBPPlane_cent[i]->SetGridy();
        hEff_NBPPlanes_cent[i]->SetStats(0);
        hEff_NBPPlanes_cent[i]->GetXaxis()->SetTitle("FT0C [%]");
        hEff_NBPPlanes_cent[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPPlanes_cent[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPPlanes_cent[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPPlanes_cent[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPPlanes_cent[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPPlanes_cent[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPPlanes_cent[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPPlanes_cent[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPPlanes_cent[i]->SetLineColor(kGreen+3);
        hEff_NBPPlanes_cent[i]->SetLineWidth(2);
        hEff_NBPPlanes_cent[i]->Draw("E");
    }

    //Eff vs pt
    for (int i = 0; i < nBinsPlane; i++) {
        cEffBothPlane_pt[i] = new TCanvas();
        cEffBPPlane_pt[i] = new TCanvas();
        cEffNBPPlane_pt[i] = new TCanvas();

        cEffBothPlane_pt[i]->cd();
        cEffBothPlane_pt[i]->SetGridx();
        cEffBothPlane_pt[i]->SetGridy();
        hEff_bothPlanes_pt[i]->SetStats(0);
        hEff_bothPlanes_pt[i]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hEff_bothPlanes_pt[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothPlanes_pt[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothPlanes_pt[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothPlanes_pt[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothPlanes_pt[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothPlanes_pt[i]->GetXaxis()->CenterTitle(true);
        hEff_bothPlanes_pt[i]->GetYaxis()->CenterTitle(true);
        hEff_bothPlanes_pt[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothPlanes_pt[i]->SetLineColor(kBlack);
        hEff_bothPlanes_pt[i]->SetLineWidth(2);
        hEff_bothPlanes_pt[i]->Draw("E");

        cEffBPPlane_pt[i]->cd();
        cEffBPPlane_pt[i]->SetGridx();
        cEffBPPlane_pt[i]->SetGridy();
        hEff_BPPlanes_pt[i]->SetStats(0);
        hEff_BPPlanes_pt[i]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hEff_BPPlanes_pt[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPPlanes_pt[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPPlanes_pt[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPPlanes_pt[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPPlanes_pt[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPPlanes_pt[i]->GetXaxis()->CenterTitle(true);
        hEff_BPPlanes_pt[i]->GetYaxis()->CenterTitle(true);
        hEff_BPPlanes_pt[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPPlanes_pt[i]->SetLineColor(kRed);
        hEff_BPPlanes_pt[i]->SetLineWidth(2);
        hEff_BPPlanes_pt[i]->Draw("E");

        cEffNBPPlane_pt[i]->cd();
        cEffNBPPlane_pt[i]->SetGridx();
        cEffNBPPlane_pt[i]->SetGridy();
        hEff_NBPPlanes_pt[i]->SetStats(0);
        hEff_NBPPlanes_pt[i]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hEff_NBPPlanes_pt[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPPlanes_pt[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPPlanes_pt[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPPlanes_pt[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPPlanes_pt[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPPlanes_pt[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPPlanes_pt[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPPlanes_pt[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPPlanes_pt[i]->SetLineColor(kGreen+3);
        hEff_NBPPlanes_pt[i]->SetLineWidth(2);
        hEff_NBPPlanes_pt[i]->Draw("E");
    }

    //Eff vs eta
    for (int i = 0; i < nBinsPlane; i++) {
        cEffBothPlane_eta[i] = new TCanvas();
        cEffBPPlane_eta[i] = new TCanvas();
        cEffNBPPlane_eta[i] = new TCanvas();

        cEffBothPlane_eta[i]->cd();
        cEffBothPlane_eta[i]->SetGridx();
        cEffBothPlane_eta[i]->SetGridy();
        hEff_bothPlanes_eta[i]->SetStats(0);
        hEff_bothPlanes_eta[i]->GetXaxis()->SetTitle("#eta");
        hEff_bothPlanes_eta[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothPlanes_eta[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothPlanes_eta[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothPlanes_eta[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothPlanes_eta[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothPlanes_eta[i]->GetXaxis()->CenterTitle(true);
        hEff_bothPlanes_eta[i]->GetYaxis()->CenterTitle(true);
        hEff_bothPlanes_eta[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothPlanes_eta[i]->SetLineColor(kBlack);
        hEff_bothPlanes_eta[i]->SetLineWidth(2);
        hEff_bothPlanes_eta[i]->Draw("E");

        cEffBPPlane_eta[i]->cd();
        cEffBPPlane_eta[i]->SetGridx();
        cEffBPPlane_eta[i]->SetGridy();
        hEff_BPPlanes_eta[i]->SetStats(0);
        hEff_BPPlanes_eta[i]->GetXaxis()->SetTitle("#eta");
        hEff_BPPlanes_eta[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPPlanes_eta[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPPlanes_eta[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPPlanes_eta[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPPlanes_eta[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPPlanes_eta[i]->GetXaxis()->CenterTitle(true);
        hEff_BPPlanes_eta[i]->GetYaxis()->CenterTitle(true);
        hEff_BPPlanes_eta[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPPlanes_eta[i]->SetLineColor(kRed);
        hEff_BPPlanes_eta[i]->SetLineWidth(2);
        hEff_BPPlanes_eta[i]->Draw("E");

        cEffNBPPlane_eta[i]->cd();
        cEffNBPPlane_eta[i]->SetGridx();
        cEffNBPPlane_eta[i]->SetGridy();
        hEff_NBPPlanes_eta[i]->SetStats(0);
        hEff_NBPPlanes_eta[i]->GetXaxis()->SetTitle("#eta");
        hEff_NBPPlanes_eta[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPPlanes_eta[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPPlanes_eta[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPPlanes_eta[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPPlanes_eta[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPPlanes_eta[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPPlanes_eta[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPPlanes_eta[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPPlanes_eta[i]->SetLineColor(kGreen+3);
        hEff_NBPPlanes_eta[i]->SetLineWidth(2);
        hEff_NBPPlanes_eta[i]->Draw("E");
    }

    //Eff vs phi
    for (int i = 0; i < nBinsPlane; i++) {
        cEffBothPlane_phi[i] = new TCanvas();
        cEffBPPlane_phi[i] = new TCanvas();
        cEffNBPPlane_phi[i] = new TCanvas();

        cEffBothPlane_phi[i]->cd();
        cEffBothPlane_phi[i]->SetGridx();
        cEffBothPlane_phi[i]->SetGridy();
        hEff_bothPlanes_phi[i]->SetStats(0);
        hEff_bothPlanes_phi[i]->GetXaxis()->SetTitle("#phi [rad]");
        hEff_bothPlanes_phi[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothPlanes_phi[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothPlanes_phi[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothPlanes_phi[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothPlanes_phi[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothPlanes_phi[i]->GetXaxis()->CenterTitle(true);
        hEff_bothPlanes_phi[i]->GetYaxis()->CenterTitle(true);
        hEff_bothPlanes_phi[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothPlanes_phi[i]->SetLineColor(kBlack);
        hEff_bothPlanes_phi[i]->SetLineWidth(2);
        hEff_bothPlanes_phi[i]->Draw("E");

        cEffBPPlane_phi[i]->cd();
        cEffBPPlane_phi[i]->SetGridx();
        cEffBPPlane_phi[i]->SetGridy();
        hEff_BPPlanes_phi[i]->SetStats(0);
        hEff_BPPlanes_phi[i]->GetXaxis()->SetTitle("phi [rad]");
        hEff_BPPlanes_phi[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPPlanes_phi[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPPlanes_phi[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPPlanes_phi[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPPlanes_phi[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPPlanes_phi[i]->GetXaxis()->CenterTitle(true);
        hEff_BPPlanes_phi[i]->GetYaxis()->CenterTitle(true);
        hEff_BPPlanes_phi[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPPlanes_phi[i]->SetLineColor(kRed);
        hEff_BPPlanes_phi[i]->SetLineWidth(2);
        hEff_BPPlanes_phi[i]->Draw("E");

        cEffNBPPlane_phi[i]->cd();
        cEffNBPPlane_phi[i]->SetGridx();
        cEffNBPPlane_phi[i]->SetGridy();
        hEff_NBPPlanes_phi[i]->SetStats(0);
        hEff_NBPPlanes_phi[i]->GetXaxis()->SetTitle("#phi [rad]");
        hEff_NBPPlanes_phi[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPPlanes_phi[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPPlanes_phi[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPPlanes_phi[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPPlanes_phi[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPPlanes_phi[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPPlanes_phi[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPPlanes_phi[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPPlanes_phi[i]->SetLineColor(kGreen+3);
        hEff_NBPPlanes_phi[i]->SetLineWidth(2);
        hEff_NBPPlanes_phi[i]->Draw("E");
    }

    //Per RPC
    //Eff vs centrality
    for (int i = 0; i < nBinsRPC; i++) {
        cEffBothRPC_cent[i] = new TCanvas();
        cEffBPRPC_cent[i] = new TCanvas();
        cEffNBPRPC_cent[i] = new TCanvas();

        cEffBothRPC_cent[i]->cd();
        cEffBothRPC_cent[i]->SetGridx();
        cEffBothRPC_cent[i]->SetGridy();
        hEff_bothRPC_cent[i]->SetStats(0);
        hEff_bothRPC_cent[i]->GetXaxis()->SetTitle("FT0C [%]");
        hEff_bothRPC_cent[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothRPC_cent[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothRPC_cent[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothRPC_cent[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothRPC_cent[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothRPC_cent[i]->GetXaxis()->CenterTitle(true);
        hEff_bothRPC_cent[i]->GetYaxis()->CenterTitle(true);
        hEff_bothRPC_cent[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothRPC_cent[i]->SetLineColor(kBlack);
        hEff_bothRPC_cent[i]->SetLineWidth(2);
        hEff_bothRPC_cent[i]->Draw("E");

        cEffBPRPC_cent[i]->cd();
        cEffBPRPC_cent[i]->SetGridx();
        cEffBPRPC_cent[i]->SetGridy();
        hEff_BPRPC_cent[i]->SetStats(0);
        hEff_BPRPC_cent[i]->GetXaxis()->SetTitle("FT0C [%]");
        hEff_BPRPC_cent[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPRPC_cent[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPRPC_cent[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPRPC_cent[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPRPC_cent[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPRPC_cent[i]->GetXaxis()->CenterTitle(true);
        hEff_BPRPC_cent[i]->GetYaxis()->CenterTitle(true);
        hEff_BPRPC_cent[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPRPC_cent[i]->SetLineColor(kRed);
        hEff_BPRPC_cent[i]->SetLineWidth(2);
        hEff_BPRPC_cent[i]->Draw("E");

        cEffNBPRPC_cent[i]->cd();
        cEffNBPRPC_cent[i]->SetGridx();
        cEffNBPRPC_cent[i]->SetGridy();
        hEff_NBPRPC_cent[i]->SetStats(0);
        hEff_NBPRPC_cent[i]->GetXaxis()->SetTitle("FT0C [%]");
        hEff_NBPRPC_cent[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPRPC_cent[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPRPC_cent[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPRPC_cent[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPRPC_cent[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPRPC_cent[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPRPC_cent[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPRPC_cent[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPRPC_cent[i]->SetLineColor(kGreen+3);
        hEff_NBPRPC_cent[i]->SetLineWidth(2);
        hEff_NBPRPC_cent[i]->Draw("E");
    }

    //Eff vs pt
    for (int i = 0; i < nBinsRPC; i++) {
        cEffBothRPC_pt[i] = new TCanvas();
        cEffBPRPC_pt[i] = new TCanvas();
        cEffNBPRPC_pt[i] = new TCanvas();

        cEffBothRPC_pt[i]->cd();
        cEffBothRPC_pt[i]->SetGridx();
        cEffBothRPC_pt[i]->SetGridy();
        hEff_bothRPC_pt[i]->SetStats(0);
        hEff_bothRPC_pt[i]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hEff_bothRPC_pt[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothRPC_pt[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothRPC_pt[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothRPC_pt[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothRPC_pt[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothRPC_pt[i]->GetXaxis()->CenterTitle(true);
        hEff_bothRPC_pt[i]->GetYaxis()->CenterTitle(true);
        hEff_bothRPC_pt[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothRPC_pt[i]->SetLineColor(kBlack);
        hEff_bothRPC_pt[i]->SetLineWidth(2);
        hEff_bothRPC_pt[i]->Draw("E");

        cEffBPRPC_pt[i]->cd();
        cEffBPRPC_pt[i]->SetGridx();
        cEffBPRPC_pt[i]->SetGridy();
        hEff_BPRPC_pt[i]->SetStats(0);
        hEff_BPRPC_pt[i]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hEff_BPRPC_pt[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPRPC_pt[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPRPC_pt[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPRPC_pt[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPRPC_pt[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPRPC_pt[i]->GetXaxis()->CenterTitle(true);
        hEff_BPRPC_pt[i]->GetYaxis()->CenterTitle(true);
        hEff_BPRPC_pt[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPRPC_pt[i]->SetLineColor(kRed);
        hEff_BPRPC_pt[i]->SetLineWidth(2);
        hEff_BPRPC_pt[i]->Draw("E");

        cEffNBPRPC_pt[i]->cd();
        cEffNBPRPC_pt[i]->SetGridx();
        cEffNBPRPC_pt[i]->SetGridy();
        hEff_NBPRPC_pt[i]->SetStats(0);
        hEff_NBPRPC_pt[i]->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]");
        hEff_NBPRPC_pt[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPRPC_pt[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPRPC_pt[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPRPC_pt[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPRPC_pt[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPRPC_pt[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPRPC_pt[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPRPC_pt[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPRPC_pt[i]->SetLineColor(kGreen+3);
        hEff_NBPRPC_pt[i]->SetLineWidth(2);
        hEff_NBPRPC_pt[i]->Draw("E");
    }

    //Eff vs eta
    for (int i = 0; i < nBinsRPC; i++) {
        cEffBothRPC_eta[i] = new TCanvas();
        cEffBPRPC_eta[i] = new TCanvas();
        cEffNBPRPC_eta[i] = new TCanvas();

        cEffBothRPC_eta[i]->cd();
        cEffBothRPC_eta[i]->SetGridx();
        cEffBothRPC_eta[i]->SetGridy();
        hEff_bothRPC_eta[i]->SetStats(0);
        hEff_bothRPC_eta[i]->GetXaxis()->SetTitle("#eta");
        hEff_bothRPC_eta[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothRPC_eta[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothRPC_eta[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothRPC_eta[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothRPC_eta[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothRPC_eta[i]->GetXaxis()->CenterTitle(true);
        hEff_bothRPC_eta[i]->GetYaxis()->CenterTitle(true);
        hEff_bothRPC_eta[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothRPC_eta[i]->SetLineColor(kBlack);
        hEff_bothRPC_eta[i]->SetLineWidth(2);
        hEff_bothRPC_eta[i]->Draw("E");

        cEffBPRPC_eta[i]->cd();
        cEffBPRPC_eta[i]->SetGridx();
        cEffBPRPC_eta[i]->SetGridy();
        hEff_BPRPC_eta[i]->SetStats(0);
        hEff_BPRPC_eta[i]->GetXaxis()->SetTitle("#eta");
        hEff_BPRPC_eta[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPRPC_eta[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPRPC_eta[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPRPC_eta[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPRPC_eta[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPRPC_eta[i]->GetXaxis()->CenterTitle(true);
        hEff_BPRPC_eta[i]->GetYaxis()->CenterTitle(true);
        hEff_BPRPC_eta[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPRPC_eta[i]->SetLineColor(kRed);
        hEff_BPRPC_eta[i]->SetLineWidth(2);
        hEff_BPRPC_eta[i]->Draw("E");

        cEffNBPRPC_eta[i]->cd();
        cEffNBPRPC_eta[i]->SetGridx();
        cEffNBPRPC_eta[i]->SetGridy();
        hEff_NBPRPC_eta[i]->SetStats(0);
        hEff_NBPRPC_eta[i]->GetXaxis()->SetTitle("#eta");
        hEff_NBPRPC_eta[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPRPC_eta[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPRPC_eta[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPRPC_eta[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPRPC_eta[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPRPC_eta[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPRPC_eta[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPRPC_eta[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPRPC_eta[i]->SetLineColor(kGreen+3);
        hEff_NBPRPC_eta[i]->SetLineWidth(2);
        hEff_NBPRPC_eta[i]->Draw("E");
    }

    //Eff vs phi
    for (int i = 0; i < nBinsRPC; i++) {
        cEffBothRPC_phi[i] = new TCanvas();
        cEffBPRPC_phi[i] = new TCanvas();
        cEffNBPRPC_phi[i] = new TCanvas();

        cEffBothRPC_phi[i]->cd();
        cEffBothRPC_phi[i]->SetGridx();
        cEffBothRPC_phi[i]->SetGridy();
        hEff_bothRPC_phi[i]->SetStats(0);
        hEff_bothRPC_phi[i]->GetXaxis()->SetTitle("#phi [rad]");
        hEff_bothRPC_phi[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_bothRPC_phi[i]->GetXaxis()->SetTitleFont(62);
        hEff_bothRPC_phi[i]->GetXaxis()->SetLabelFont(62);
        hEff_bothRPC_phi[i]->GetYaxis()->SetTitleFont(62);
        hEff_bothRPC_phi[i]->GetYaxis()->SetLabelFont(62);
        hEff_bothRPC_phi[i]->GetXaxis()->CenterTitle(true);
        hEff_bothRPC_phi[i]->GetYaxis()->CenterTitle(true);
        hEff_bothRPC_phi[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_bothRPC_phi[i]->SetLineColor(kBlack);
        hEff_bothRPC_phi[i]->SetLineWidth(2);
        hEff_bothRPC_phi[i]->Draw("E");

        cEffBPRPC_phi[i]->cd();
        cEffBPRPC_phi[i]->SetGridx();
        cEffBPRPC_phi[i]->SetGridy();
        hEff_BPRPC_phi[i]->SetStats(0);
        hEff_BPRPC_phi[i]->GetXaxis()->SetTitle("phi [rad]");
        hEff_BPRPC_phi[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_BPRPC_phi[i]->GetXaxis()->SetTitleFont(62);
        hEff_BPRPC_phi[i]->GetXaxis()->SetLabelFont(62);
        hEff_BPRPC_phi[i]->GetYaxis()->SetTitleFont(62);
        hEff_BPRPC_phi[i]->GetYaxis()->SetLabelFont(62);
        hEff_BPRPC_phi[i]->GetXaxis()->CenterTitle(true);
        hEff_BPRPC_phi[i]->GetYaxis()->CenterTitle(true);
        hEff_BPRPC_phi[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_BPRPC_phi[i]->SetLineColor(kRed);
        hEff_BPRPC_phi[i]->SetLineWidth(2);
        hEff_BPRPC_phi[i]->Draw("E");

        cEffNBPRPC_phi[i]->cd();
        cEffNBPRPC_phi[i]->SetGridx();
        cEffNBPRPC_phi[i]->SetGridy();
        hEff_NBPRPC_phi[i]->SetStats(0);
        hEff_NBPRPC_phi[i]->GetXaxis()->SetTitle("#phi [rad]");
        hEff_NBPRPC_phi[i]->GetYaxis()->SetTitle("Efficiency [%]");
        hEff_NBPRPC_phi[i]->GetXaxis()->SetTitleFont(62);
        hEff_NBPRPC_phi[i]->GetXaxis()->SetLabelFont(62);
        hEff_NBPRPC_phi[i]->GetYaxis()->SetTitleFont(62);
        hEff_NBPRPC_phi[i]->GetYaxis()->SetLabelFont(62);
        hEff_NBPRPC_phi[i]->GetXaxis()->CenterTitle(true);
        hEff_NBPRPC_phi[i]->GetYaxis()->CenterTitle(true);
        hEff_NBPRPC_phi[i]->GetYaxis()->SetRangeUser(80,100);
        hEff_NBPRPC_phi[i]->SetLineColor(kGreen+3);
        hEff_NBPRPC_phi[i]->SetLineWidth(2);
        hEff_NBPRPC_phi[i]->Draw("E");
    }

    bool write = true;
    if(write) {
        
        TFile *fOutTot;
        
        if (centralityCut == true && ptCut == false) {
            fOutTot = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/outEffAll_" + period + "_centCut.root").c_str(),"RECREATE");
        }
        else if (centralityCut == false && ptCut == true) {
            fOutTot = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/outEffAll_" + period + "_ptCut.root").c_str(),"RECREATE");
        }
        else if (centralityCut == true && ptCut == true) {
            fOutTot = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/outEffAll_" + period + "_pt_cent_cut.root").c_str(),"RECREATE");
        }
        else if (centralityCut == false && ptCut == false) {
            fOutTot = new TFile(("/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/outEffAll_" + period + "_noCuts.root").c_str(),"RECREATE");
        }
        
        
        fOutTot->cd();

        for (int i = 0; i < nBinsPlane; i++) {

            cEffBothPlane_cent[i] ->Write(("eff_bothPlanes_Plane" + to_string(i) + "_vs_cent").c_str());
            cEffBPPlane_cent[i] ->Write(("eff_BP_Plane" + to_string(i) + "_vs_cent").c_str());
            cEffNBPPlane_cent[i] ->Write(("eff_NBP_Plane" + to_string(i) + "_vs_cent").c_str());

            cEffBothPlane_pt[i] ->Write(("eff_bothPlanes_Plane" + to_string(i) + "_vs_pt").c_str());
            cEffBPPlane_pt[i] ->Write(("eff_BP_Plane" + to_string(i) + "_vs_pt").c_str());
            cEffNBPPlane_pt[i] ->Write(("eff_NBP_Plane" + to_string(i) + "_vs_pt").c_str());

            cEffBothPlane_eta[i] ->Write(("eff_bothPlanes_Plane" + to_string(i) + "_vs_eta").c_str());
            cEffBPPlane_eta[i] ->Write(("eff_BP_Plane" + to_string(i) + "_vs_eta").c_str());
            cEffNBPPlane_eta[i] ->Write(("eff_NBP_Plane" + to_string(i) + "_vs_eta").c_str());

            cEffBothPlane_phi[i] ->Write(("eff_bothPlanes_Plane" + to_string(i) + "_vs_phi").c_str());
            cEffBPPlane_phi[i] ->Write(("eff_BP_Plane" + to_string(i) + "_vs_phi").c_str());
            cEffNBPPlane_phi[i] ->Write(("eff_NBP_Plane" + to_string(i) + "_vs_phi").c_str());
        }

        for (int i = 0; i < nBinsRPC; i++) {

            cEffBothRPC_cent[i] ->Write(("eff_bothPlanes_RPC" + to_string(i) + "_vs_cent").c_str());
            cEffBPRPC_cent[i] ->Write(("eff_BP_RPC" + to_string(i) + "_vs_cent").c_str());
            cEffNBPRPC_cent[i] ->Write(("eff_NBP_RPC" + to_string(i) + "_vs_cent").c_str());

            cEffBothRPC_pt[i] ->Write(("eff_bothPlanes_RPC" + to_string(i) + "_vs_pt").c_str());
            cEffBPRPC_pt[i] ->Write(("eff_BP_RPC" + to_string(i) + "_vs_pt").c_str());
            cEffNBPRPC_pt[i] ->Write(("eff_NBP_RPC" + to_string(i) + "_vs_pt").c_str());

            cEffBothRPC_eta[i] ->Write(("eff_bothPlanes_RPC" + to_string(i) + "_vs_eta").c_str());
            cEffBPRPC_eta[i] ->Write(("eff_BP_RPC" + to_string(i) + "_vs_eta").c_str());
            cEffNBPRPC_eta[i] ->Write(("eff_NBP_RPC" + to_string(i) + "_vs_eta").c_str());

            cEffBothRPC_phi[i] ->Write(("eff_bothPlanes_RPC" + to_string(i) + "_vs_phi").c_str());
            cEffBPRPC_phi[i] ->Write(("eff_BP_RPC" + to_string(i) + "_vs_phi").c_str());
            cEffNBPRPC_phi[i] ->Write(("eff_NBP_RPC" + to_string(i) + "_vs_phi").c_str());
        }

        fOutTot->Close();
    }

    bool save = true;
    
    if (save) {
        for (int i = 0; i < nBinsPlane; i++) {

            if (ptCut) {
                //centrality
                if (centrality) {
                    cEffBothPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_cent_pTcut.png").c_str());
                    cEffBPPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_cent_pTcut.png").c_str());
                    cEffNBPPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_cent_pTcut.png").c_str());
                }
                //pt
                if (pt) {
                    cEffBothPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_pt_pTcut.png").c_str());
                    cEffBPPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_pt_pTcut.png").c_str());
                    cEffNBPPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_pt_pTcut.png").c_str());
                }
                //eta
                if (eta) {
                    cEffBothPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_eta_pTcut.png").c_str());
                    cEffBPPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_eta_pTcut.png").c_str());
                    cEffNBPPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_eta_pTcut.png").c_str());
                }
                //phi
                if (phi) {
                    cEffBothPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_phi_pTcut.png").c_str());
                    cEffBPPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_phi_pTcut.png").c_str());
                    cEffNBPPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_phi_pTcut.png").c_str());
                }
            }

            if (ptCut && centralityCut) {
                //centrality
                if (centrality) {
                    cEffBothPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_cent_pT_cent_cut.png").c_str());
                    cEffBPPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_cent_pT_cent_cut.png").c_str());
                    cEffNBPPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_cent_pT_cent_cut.png").c_str());
                }
                //pt
                if (pt) {
                    cEffBothPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_pt_pT_cent_cut.png").c_str());
                    cEffBPPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_pt_pT_cent_cut.png").c_str());
                    cEffNBPPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_pt_pT_cent_cut.png").c_str());
                }
                //eta
                if (eta) {
                    cEffBothPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_eta_pT_cent_cut.png").c_str());
                    cEffBPPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_eta_pT_cent_cut.png").c_str());
                    cEffNBPPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_eta_pT_cent_cut.png").c_str());
                }
                //phi
                if (phi) {
                    cEffBothPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_phi_pT_cent_cut.png").c_str());
                    cEffBPPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_phi_pT_cent_cut.png").c_str());
                    cEffNBPPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_phi_pT_cent_cut.png").c_str());
                }
            }
        
            else {
                //centrality
                if (centrality) {
                    cEffBothPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_cent.png").c_str());
                    cEffBPPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_cent.png").c_str());
                    cEffNBPPlane_cent[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_cent.png").c_str());
                }
                //pt
                if (pt) {
                    cEffBothPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_pt.png").c_str());
                    cEffBPPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_pt.png").c_str());
                    cEffNBPPlane_pt[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_pt.png").c_str());
                }
                //eta
                if (eta) {
                    cEffBothPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_eta.png").c_str());
                    cEffBPPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_eta.png").c_str());
                    cEffNBPPlane_eta[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_eta.png").c_str());
                }
                //phi
                if (phi) {
                    cEffBothPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_bothPlanes_phi.png").c_str());
                    cEffBPPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_BP_planes_phi.png").c_str());
                    cEffNBPPlane_phi[i]->SaveAs(("/home/luca/cernbox/assegnoTorino/MIDefficiency/presentations/images/"+planeName[i]+"_NBP_planes_phi.png").c_str());
                }
            } //end of else on ptcut
        }
    }
    
    
    //---------------------//

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