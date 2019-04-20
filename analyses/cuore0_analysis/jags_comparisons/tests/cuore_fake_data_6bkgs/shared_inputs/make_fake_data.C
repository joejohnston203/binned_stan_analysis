#include <stdio.h>

void make_fake_data() {
  // The simulations and normalizations used to make the fake data
  static const int nFiles = 6;
  //TString simDir = "../simulations_cuore/reduced_sims/";
  TString simDir = "/nfs/cuore1/data/simulation/CUORE/2017/ntp/";
  TString filenames [nFiles] = {"2nuFIJK-100M.root",
                                "CuNOSV-k40.root",
                                "CuNOSV-co60.root",
                                "TopPb5cm_210Bi_cnaf136.root",
                                "TeO2Sx-pb210-.001_cnaf14.root",
                                "CuNOSVTowerSx-pb210-.01.root"};
  double normalization [nFiles] = {1.86e-3,
                                   9.52e-3,
                                   5.13e-4,
                                   7.46e-3,
                                   1.23e-2,
                                   1.18e-2};

  TString outFileName = "fake_data_cuore_6_sim.root";
  TString treeName = "outTree";
  TString tempFilePrefix = "temp_file_delete_me_";

  
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
    newTrees[i]->SetAlias("Energy", "Ener2");
    newTrees[i]->SetAlias("TotalEnergy", "ESum2");
    newTrees[i]->SetAlias("PSA", "0*Ener2+1");
    newTrees[i]->SetAlias("Included", "0*Ener2+1");
    newTrees[i]->SetAlias("Multiplicity", "Multiplicity");
    newTrees[i]->SetAlias("MultipletIndex", "MultipletIndex");
    tempfiles[i]->Write();
    files[i]->Close();
    tempfiles[i]->Close();
  }

  TFile *newfile = new TFile(outFileName,"recreate");
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
