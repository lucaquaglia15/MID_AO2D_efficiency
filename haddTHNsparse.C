/// \file
/// \ingroup tutorial_io
/// \notebook -nodraw
/// Macro to add histogram files
/// This macro is kept for didactical purposes only: use instead the executable $ROOTSYS/bin/hadd !
/// 
/// This macro will add histograms from a list of root files and write them
/// to a target root file. The target file is newly created and must not be
/// identical to one of the source files.
/// This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
/// which had a problem with directories more than one level deep.
/// The macro from Sven has been enhanced by Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
/// to automatically add Trees (via a chain of trees).
///
/// \macro_code
///
/// \author Sven A. Schmidt, sven.schmidt@cern.ch, 13.2.2001

/*
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TError.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );


void haddTHNsparse(const string mergedfile, const string filename_list) {

    //Only show fatal errors and not other warnings
    gErrorIgnoreLevel = kFatal;

    TH1::AddDirectory(kFALSE);  // histograms not kept in memory by ROOT

    // Prepare the files to me merged
    if(gSystem->AccessPathName("hsimple1.root")) {
    gSystem->CopyFile("hsimple.root", "hsimple1.root");
    gSystem->CopyFile("hsimple.root", "hsimple2.root");
    }

    // in an interactive ROOT session, edit the file names
    // Target and FileList, then
    // root > .L hadd.C
    // root > hadd()


    //const char* mergedfile = "Scan000775_HV1_CAEN.root"; 

    string line;
    ifstream myfile;
    myfile.open(filename_list.c_str());

    if(!myfile.is_open()) {
    	perror("Error open");
    	exit(EXIT_FAILURE);
    }

    Target = TFile::Open( mergedfile.c_str() , "RECREATE" );
    FileList = new TList();

    //cout << "File to merge: \n";
    while(getline(myfile, line)) {
    	//cout << "\t" <<line << endl;
    	FileList->Add( TFile::Open(line.c_str()) );
    }

//    Target = TFile::Open( "result.root", "RECREATE" );
//    FileList = new TList();
    // FileList->Add( TFile::Open("3DCZT_prova.root") );
    // FileList->Add( TFile::Open("GunAction_4_60mm_4Cryst_Phan_Tissue_n0_5_2021_06_20_13_46_02.root") );

    MergeRootfile( Target, FileList );

    // --- memory cleanup ---
    TIter nextfile(FileList);
    TFile* f;
    while ((f = (TFile*)nextfile())) {
        f->Close();
        delete f;
    }
}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {

    //  cout << "Target path: " << target->GetPath() << endl;
    TString path( (char*)strstr( target->GetPath(), ":" ) );
    path.Remove( 0, 2 );

    TFile *first_source = (TFile*)sourcelist->First();
    first_source->cd( path );
    TDirectory *current_sourcedir = gDirectory;
    //gain time, do not add the objects in the list in memory
    Bool_t status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // loop over all keys in this directory
    TChain *globChain = 0;
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key, *oldkey=0;
    while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
        // descendant of TH1 -> merge it

        //      cout << "Merging histogram " << obj->GetName() << endl;
        TH1 *h1 = (TH1*)obj;

        // loop over all source files and add the content of the
        // correspondant histogram to the one pointed to by "h1"
        TFile *nextsource = (TFile*)sourcelist->After( first_source );
        while ( nextsource ) {

            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
                TH1 *h2 = (TH1*)key2->ReadObj();
                h1->Add( h2 );
                delete h2;
            }

            nextsource = (TFile*)sourcelist->After( nextsource );
        }
    }
    else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

        // loop over all source files create a chain of Trees "globChain"
        const char* obj_name= obj->GetName();

        globChain = new TChain(obj_name);
        globChain->Add(first_source->GetName());
        TFile *nextsource = (TFile*)sourcelist->After( first_source );
        //      const char* file_name = nextsource->GetName();
        // cout << "file name  " << file_name << endl;
        while ( nextsource ) {

            globChain->Add(nextsource->GetName());
            nextsource = (TFile*)sourcelist->After( nextsource );
        }

    } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
        // it's a subdirectory

        //cout << "Found subdirectory " << obj->GetName() << endl;

        // create a new subdir of same name and title in the target file
        target->cd();
        TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

        // newdir is now the starting point of another round of merging
        // newdir still knows its depth within the target file via
        // GetPath(), so we can still figure out where we are in the recursion
        MergeRootfile( newdir, sourcelist );

    } else if ( obj->IsA()->InheritsFrom(THnSparse::Class()) ) {
        THnSparse *h1 = (THnSparse*)obj;

        TFile *nextsource = (TFile*)sourcelist->After( first_source );
        while ( nextsource ) {
            nextsource->cd(path);
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
                THnSparse *h2 = (THnSparse*)key2->ReadObj();
                h1->Add(h2);  // THnSparse has Add() just like TH1
                delete h2;
            }
            nextsource = (TFile*)sourcelist->After( nextsource );
        }
    } else {

        // object is of no type that we know or can handle
        //cout << "Unknown object type, name: "
        //<< obj->GetName() << " title: " << obj->GetTitle() << endl;
    }

    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
        target->cd();

        //!!if the object is a tree, it is stored in globChain...
        if(obj->IsA()->InheritsFrom( TTree::Class() ))
            globChain->Merge(target->GetFile(),0,"keep");
        else
            obj->Write( key->GetName() );
    }

    } // while ( ( TKey *key = (TKey*)nextkey() ) )

    // save modifications to target file
    target->SaveSelf(kTRUE);
    TH1::AddDirectory(status);
}*/

//chatgpt versione

#include <string>
#include <vector>
#include <fstream>
#include "TFile.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TChain.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TError.h"

using namespace std;

// Forward declaration
void MergeRootfile(TDirectory *target, const vector<string> &filenames, const TString &path="");

void haddTHNsparse(const string &mergedfile, const string &filename_list) {
    // Only show fatal errors
    //gErrorIgnoreLevel = kFatal;
    
    TH1::AddDirectory(kFALSE);

    // Read the filenames
    vector<string> filenames;
    ifstream myfile(filename_list);
    if (!myfile.is_open()) {
        perror("Error opening filename list");
        exit(EXIT_FAILURE);
    }
    string line;
    while (getline(myfile, line)) {
        if (!line.empty()) filenames.push_back(line);
    }
    myfile.close();

    // Open target file
    TFile* Target = TFile::Open(mergedfile.c_str(), "RECREATE");
    if (!Target || Target->IsZombie()) {
        cerr << "Error creating target file " << mergedfile << endl;
        return;
    }

    // Start recursive merge from root directory
    MergeRootfile(Target, filenames);

    // Save and close
    Target->Write();
    Target->Close();
    delete Target;

    // Force ROOT cleanup
    //gROOT->GetListOfFiles()->Delete();
    //gROOT->GetListOfCleanups()->Delete();
}

// ---------------- MergeRootfile ----------------
void MergeRootfile(TDirectory *target, const vector<string> &filenames, const TString &path) {

    // Open the first file to get the keys in this directory
    TFile* first_source = TFile::Open(filenames[0].c_str());
    if (!first_source || first_source->IsZombie()) return;

    first_source->cd(path);
    TDirectory* current_sourcedir = gDirectory;

    TIter nextkey(current_sourcedir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;

        TString name = obj->GetName();
        cout << "Name is: " << name << endl;

        // ---------------- TH1 ----------------
        if (obj->InheritsFrom(TH1::Class())) {
            TH1* h1 = (TH1*)obj;
            cout << "I am in TH1" << endl;

            // Loop over remaining files and add
            for (size_t i = 1; i < filenames.size(); ++i) {
                TFile* f = TFile::Open(filenames[i].c_str());
                if (!f || f->IsZombie()) continue;
                cout << "path in TH1 " << path << endl;
                f->cd(path);
                TKey* k2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(name);
                if (k2) {
                    TH1* h2 = (TH1*)k2->ReadObj();
                    h1->Add(h2);
                    delete h2;
                }
                f->Close();
                delete f;
            }

            cout << "I am at the end of TH1" << endl;
            target->cd();
            h1->Write(name);
            delete h1;  // free accumulator
            cout << "Deleted TH1 in TH1" << endl;

        }
        // ---------------- THnSparse ----------------
        else if (obj->InheritsFrom(THnSparse::Class())) {
            THnSparse* h1 = (THnSparse*)obj;
            cout << "I am in THNsparse" << endl;

            for (size_t i = 1; i < filenames.size(); ++i) {
                TFile* f = TFile::Open(filenames[i].c_str());
                if (!f || f->IsZombie()) continue;
                f->cd(path);
                cout << "path in THNSparse " << path << endl;
                TKey* k2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(name);
                if (k2) {
                    THnSparse* h2 = (THnSparse*)k2->ReadObj();
                    h1->Add(h2);
                    delete h2;
                }
                f->Close();
                delete f;
            }

            target->cd();
            h1->Write(name);
            delete h1;

        }
        // ---------------- TTree ----------------
        else if (obj->InheritsFrom(TTree::Class())) {
            const char* tname = obj->GetName();
            TChain* chain = new TChain(tname);
            for (const auto &fname : filenames) chain->Add(fname.c_str());

            chain->Merge(target->GetFile(), 0, "keep");
            delete chain;
        }
        // ---------------- TDirectory ----------------
        else if (obj->InheritsFrom(TDirectory::Class())) {
            TString subdir_name = obj->GetName();
            target->cd();
            TDirectory* newdir = target->mkdir(subdir_name, obj->GetTitle());
            // Recursive call
            MergeRootfile(newdir, filenames, path + "/" + subdir_name);
        }

        cout << "End of loop" << endl;
    }

    first_source->Close();
    delete first_source;
}

