#include <stdio.h>

void make_fake_data() {
  // The simulations and normalizations used to make the fake data
  static const int nFiles = 4;
  TString simDir = "simulations/";
  TString filenames [nFiles] = {"10mK-co60_cnaf113.root",
				"10mK-k40_cnaf121.root",
				"10mKFlan-th232_cnaf103.root",
				"10mKFlan-u238_cnaf104.root"};
  double normalization [nFiles] = {0.00005, 0.001, 0.05, 0.1};

  TString treeName = "outTree";
  TString tempFilePrefix = "temp_tree_less_elts_";

  
  TFile * files [nFiles];
  TFile * tempfiles [nFiles];
  TTree * oldTrees [nFiles];
  TTree * newTrees [nFiles];
  int counts [nFiles];
  for(int i=0; i<nFiles; i++){
    files[i] = new TFile(simDir+filenames[i]);
    oldTrees[i] = (TTree*) files[i]->Get("outTree");
    counts[i] = oldTrees[i]->GetEntries()* normalization[i];
    cout << "File " << filenames[i]
	 << ", nEntries=" << oldTrees[i]->GetEntries()
	 << ", Norm=" << normalization[i]
	 << ", Counts=" << counts[i]
	 << std::endl;
    // ROOT somehow does not let you put trees into a chain, so we need
    // to save the trees to file, then open those files...
    TString tempfilename(tempFilePrefix);
    tempfilename = tempfilename + Form("%i",i) + ".root";
    tempfiles[i] = new TFile(tempfilename, "recreate");
    newTrees[i] = oldTrees[i]->CloneTree(counts[i]);
    tempfiles[i]->Write();
    files[i]->Close();
    tempfiles[i]->Close();
  }

  TFile *newfile = new TFile("fake_data.root","recreate");
  TChain *ch = new TChain(treeName);
  
  for(int i=0; i<nFiles; i++){
    if(counts[i]>0){
      TString tempfilename(tempFilePrefix);
      tempfilename = tempfilename + Form("%i",i) + ".root";
      ch->Add(tempfilename);
    }
  }

  ch->CloneTree();
  newfile->Write();
  
  for(int i=0; i<nFiles; i++){
    TString tempfilename(tempFilePrefix);
    tempfilename = tempfilename + Form("%i",i) + ".root";
    remove(tempfilename);
  }

  delete newfile;
}
