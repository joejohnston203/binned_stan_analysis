void gen_cns_time(){
  TString filename = "cns_time.root";
  TFile * myfile = TFile::Open(filename, "RECREATE");
  TTree mytree("cns_time", "cns shape tree");

  Float_t cns_time;
  mytree.Branch("cns_time", &cns_time, "cns_time/F");

  TRandom3 * rd = new TRandom3();
  
  int nhigh = 10000;
  int nlow = 5000;
  int years = 1;

  for(int y=0; y<years; y++){
    double start = y*365.0;
    double mid = start+365.0/2.0;
    double end = (y+1)*365.0;
    for(int i=0; i<nlow; i++){
      cns_time = start+(mid-start)*rd->Rndm();
      mytree.Fill();
    }
    for(int j=0; j<nhigh; j++){
      cns_time = mid+(end-mid)*rd->Rndm();
      mytree.Fill();
    }
  }
  mytree.Write();
  myfile->Close();
}
