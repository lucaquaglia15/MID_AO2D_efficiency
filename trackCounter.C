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

void trackCounter() {

    int totEntries = 0; //total number of entries
    int totFolders = 0; //total number of folders

    TFile *f = new TFile("AO2D.root","READ"); //Open AO2D.root file

    TIter keyList(f->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)keyList())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TDirectoryFile")) 
            continue;
        TDirectoryFile *d = (TDirectoryFile*)key->ReadObj();
        d->cd();
        TTree *t = (TTree*)d->Get("O2midefftracks");
        int nEntries = t->GetEntries();
        totEntries += nEntries; //Sum up all the entries
        totFolders++; //Increase by one the number of folders
        cout << "Folder name:" << key->GetTitle() << " entires in the tree: " << nEntries << endl;
    }
    f->Close();
    cout << "Total entries = " << totEntries << " total folders = " << totFolders << " done!" << endl;
}