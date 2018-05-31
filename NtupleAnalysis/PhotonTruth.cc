/**
   This program produces energy response plots from Monte-Carlo simulations
*/

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TGraphAsymmErrors.h>

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];
    
  dummyv[0] = strdup("main");

  //Histogram Binning
  const int nbinseta = 10;
  Double_t etabins[nbinseta+1] = {};
  double etamin = -0.9;
  double etamax = 0.9;
  double etastep = (etamax-etamin)/nbinseta;
  for(int i=0; i<nbinseta+1; i++){
    etabins[i] = etamin + i*etastep;
  }
  
  const int nbinsphi = 40;
  Double_t phibins[nbinsphi+1] = {};
  double phimin = -1.0*TMath::Pi();
  double phimax = 1.0*TMath::Pi();
  double phistep = (phimax-phimin)/nbinsphi;
  for(int i=0; i<nbinsphi+1; i++){
    phibins[i] = phimin + i*phistep;
  }
  
  const int nbinstrack = 20;
  Double_t trackbins[nbinstrack+1] = {};
  double ptmin = 1;
  double ptmax = 21;
  double ptstep = (ptmax-ptmin)/nbinstrack;
  for(int i = 0; i < nbinstrack+1; i++){
    trackbins[i] = ptmin + i*ptstep;
  }
  
  TH1D h_Den("h_Den", "truth photons", 20, 10, 30);    
  TH1D h_Num("h_Num", "reco photons filled with truthpt reco", 20, 10, 30);    
  TH1D h_Num2("h_Num2", "reco photons filled with pt reco", 20, 10, 30);    
  TH2D h_Correlation("h_Correlation", "", 20, 10, 30, 20, 10, 30);    

  h_Den.Sumw2();
  h_Num.Sumw2();
  h_Num2.Sumw2(); 

  h_Den.SetTitle("; p_{T}^{truth} [GeV]; entries");
  h_Num.SetTitle("; p_{T}^{Reco,truthpt} [GeV]; entries");
  h_Num2.SetTitle("; p_{T}^{Reco,pt} [GeV]; entries");
  h_Correlation.SetTitle("; True p_{T} [GeV]; Reconstructed p_{T} [GeV]");

  TH1F hDen("hDen", "", nbinstrack, trackbins);
  TH1F hNum("hNum","", nbinstrack, trackbins);
  TH1F hNum2("hNum2","",nbinstrack, trackbins);
  TH1F hReco("hReco","", nbinstrack,trackbins);
  TH1F hFake("hFake", "", nbinstrack, trackbins);
  TH1F hTrackCut("hTrackCut", "", 10, -0.5, 9.5);
  TH2F hCorrelation("hCorrelation", "", nbinstrack, trackbins, nbinstrack, trackbins);
  TH2F hRes_Pt("hRes_Pt", "", 200, 0, 10.0, 80, -50, 50);  
  
  hDen.Sumw2();
  hNum.Sumw2();
  hNum2.Sumw2();
  hFake.Sumw2();
  hReco.Sumw2();

  hDen.SetTitle("; p_{T}^{true} [GeV/c]; entries");
  hNum.SetTitle("; p_{T}^{Reco,Embed} [GeV/c]; entries");
  hNum2.SetTitle("; p_{T}^{Reco} [GeV/c]; entries");
  hFake.SetTitle("; p_{T}^{reco} [GeV/c]; Fake Rate");
  hCorrelation.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
 
  int TrackBit = 16;
  int numEvents_tracks = 0;
  int numEvents_clusters = 0;
  TApplication application("", &dummyc, dummyv);
  TCanvas* canvas = new TCanvas();
  
  //Looping over ntuples
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile *file = TFile::Open((TString)argv[iarg]);
        
    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
        
    // Get all the TTree variables from the file to open, I guess
    TTree *_tree_event = NULL;
    _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
    if (_tree_event == NULL) {
      std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
      _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " fail " << std::endl;
	exit(EXIT_FAILURE);
      }
    } 
	//TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
        
    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }
            
    //you define variables
    Double_t primary_vertex[3];
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
    Float_t track_dca_xy[NTRACK_MAX];
    Float_t track_dca_z[NTRACK_MAX];
    Float_t track_its_chi_square[NTRACK_MAX];
    UChar_t track_its_ncluster[NTRACK_MAX];
    unsigned short track_mc_truth_index[NTRACK_MAX];
        
    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX];
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];   
 
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];
    Float_t cell_e[17664];
        
    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
        
    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    Float_t eg_cross_section;
    Int_t   eg_ntrial;
    
    
        
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);
    _tree_event->SetBranchAddress("track_its_ncluster", track_its_ncluster);
    _tree_event->SetBranchAddress("track_dca_xy", track_dca_xy);
    _tree_event->SetBranchAddress("track_dca_z", track_dca_z);
    _tree_event->SetBranchAddress("track_its_chi_square", track_its_chi_square);
    _tree_event->SetBranchAddress("track_mc_truth_index", track_mc_truth_index);
        
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);        
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
        
    _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
    _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);        
    _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",mc_truth_first_parent_pdg_code);
    _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
    

    const double maxEta = 0.8;
    Long64_t numEntries = 1000000;//_tree_event->GetEntries();
    // Loop over events
    for(Long64_t ievent = 0; ievent < numEntries ; ievent++){
    // for(Long64_t ievent = 0; ievent < 10000 ; ievent++){
      _tree_event->GetEntry(ievent);
      
      bool eventChange = true;
      double weight = (double)eg_cross_section/(double)eg_ntrial;

      //loop over tracks
      for (int n = 0;  n< ntrack; n++){
	//Track Cuts
	/*
	  track quality cut
	  eta cut
	  pt cut
	  trigger selection cut
	  its cluster chi^2 cut
	  its cluster cut
	  DCA_r = DCA_xy cut
	  DCA_z cut
	 */
	hTrackCut.Fill(0);
	if((track_quality[n]&TrackBit)==0) continue; hTrackCut.Fill(1);//track quality cut
	if(TMath::Abs(track_eta[n])> maxEta) continue; hTrackCut.Fill(2);//eta cut
	if(track_pt[n] < 0.15) continue; hTrackCut.Fill(3);//pt cut
	if(track_its_chi_square[n]>36.0) continue; hTrackCut.Fill(5);//its cluster chi^2 cut
	if(TrackBit == 16)
	  {
	    if(track_its_ncluster[n] < 5) continue; 
	    hTrackCut.Fill(6);//its cluster cut
	  }
	if(TMath::Abs(track_dca_xy[n]) > 2.4) continue; hTrackCut.Fill(7);//distance of closest approach cut
	if(TMath::Abs(track_dca_z[n]) > 3.2) continue; hTrackCut.Fill(8);//distance of closest approach cut
	
	hReco.Fill(track_pt[n],weight);

	unsigned short index = track_mc_truth_index[n];
	//particles not associated with MC particle (i.e, secondaries or fakes)
	if(index>65534){ 
	  hFake.Fill(track_pt[n],weight);
	}//end if noMCParticle
	
	//particles associated with MC particle
	if(index<65534){ 
	  if((TMath::Abs(mc_truth_pdg_code[index])!= 211)  && 
	     (TMath::Abs(mc_truth_pdg_code[index])!=321) && 
	     (TMath::Abs(mc_truth_pdg_code[index])!=2212)) continue;//*/
	  if (eventChange) {numEvents_tracks++; eventChange = false;}
	  
	  hNum.Fill(mc_truth_pt[index],weight);
	  hNum2.Fill(track_pt[n],weight);
	  hCorrelation.Fill(mc_truth_pt[index], track_pt[n], weight);
	  hRes_Pt.Fill(mc_truth_pt[index], 100*(track_pt[n]-mc_truth_pt[index])/(mc_truth_pt[index]),weight);
	}//end if hasMCParticle
      }//end track loop
      


      //Loop over MC particles (all are primaries), pick charged ones with |eta|<0.8
      for (int nTru = 0;  nTru< nmc_truth; nTru++){
        int pdgcode = mc_truth_pdg_code[nTru];
        if((TMath::Abs(mc_truth_pdg_code[nTru])!= 211) && 
	   (TMath::Abs(mc_truth_pdg_code[nTru])!=321) && 
	   (TMath::Abs(mc_truth_pdg_code[nTru])!=2212)) continue;//*/
	//if(int(mc_truth_charge[n])==0) continue;
	//if(TMath::Abs(mc_truth_pdg_code[index])!= 211)  continue;
        if(TMath::Abs(mc_truth_eta[nTru])> maxEta) continue; //skip particles with |eta|>0.8
        
	
	hDen.Fill(mc_truth_pt[nTru],weight);
      }//end loop over MC particles



      eventChange = true;
      //loop over clusters
      for (ULong64_t n = 0; n < ncluster; n++) {
        //Photon Selection
        //if( not(cluster_pt[n]>8)) {continue;} //select pt of photons
        if( not(cluster_ncell[n]>2)) continue;   //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>0.05)) continue; //removes "spiky" clusters
        if( not(cluster_nlocal_maxima[n]<= 2)) continue; //require to have at most 2 local maxima.
	if( not(cluster_distance_to_bad_channel[n]>=2.0)) continue;

        //Isolation and shower shape selection:
        if( not(cluster_iso_its_04[n] < 1.0)) continue;
        if( not(cluster_lambda_square[n][0]<0.27)) continue; //single-photon selection (as opposed to merged photon).
	// Access the corresonding mc_truth particle; skip if index is 65535, which is invalid, 
        //or the truth particle pT is less than 10, or the mc_truth_pdg_code is not 22 (it's not a photon)
	//std::cout << " cluster pt " << cluster_pt[n] << " phi" << cluster_phi[n] << " eta " << cluster_eta[n] << std::endl;
	Bool_t isTruePhoton = false;
        Float_t truth_pt = -999.0;
	for(int counter = 0 ; counter<32; counter++){
	  unsigned short index = cluster_mc_truth_index[n][counter];                   

          if(isTruePhoton) break;
          if(index==65535) continue;
	  //std::cout<<"truth, pt: " << mc_truth_pt[index] << "phi " << mc_truth_phi[index] << " eta " << mc_truth_eta[index] << 
	  //  " code: " << mc_truth_pdg_code[index] << " status " << int(mc_truth_status[index]) << " parentpdg " << mc_truth_first_parent_pdg_code[index] << std::endl; 
          if(mc_truth_pdg_code[index]!=22) continue;
          if(mc_truth_first_parent_pdg_code[index]!=22) continue;
          if( not (mc_truth_status[index] >0)) continue;        
          isTruePhoton = true;
          truth_pt     = mc_truth_pt[index];
	}//end loop over indices
	
	if(isTruePhoton){
	  //fill in this histogram only photons that can be traced to a generated non-decay photon.	
            h_Num.Fill(truth_pt);
            h_Num2.Fill(cluster_pt[n]); 
	    h_Correlation.Fill(truth_pt, cluster_pt[n]);
	    if (eventChange) {numEvents_clusters++; eventChange = false;}
	} 
      }//end loop on clusters
       

      


      //loop over truth particles
      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
        //if(mc_truth_pt[nmc]<5.0) continue;
	if(mc_truth_pdg_code[nmc]==22 && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==22){
	  //std::cout << mc_truth_pt[nmc] << "phi " << mc_truth_phi[nmc] << " eta " << mc_truth_eta[nmc] << 
	  //  " code: " << mc_truth_pdg_code[nmc] << " status " << int(mc_truth_status[nmc]) << " parentpdg " << mc_truth_first_parent_pdg_code[nmc] << std::endl;    
	  h_Den.Fill(mc_truth_pt[nmc]);
	}
      } //end loop over mc truth particles
        
      // std::cout<<" ----------- "<< std::endl;
     // Create the file label, to be used within the filenames, to represent the source file
      std::string opened_files = "";
      for (int iarg = 1; iarg < argc; iarg++) {
      std::string filepath = argv[iarg];
            
      opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
      }

      h_Den.SetLineColor(2);
      THStack* hs = new THStack("hs","stack histo for plotting");
      hs->Add(&h_Num2);
      hs->Add(&h_Den);

      if (ievent % 10000 == 0) {
	hs->Draw("e1x0nostack");
	std::cout << ievent << " " << _tree_event->GetEntries() << std::endl;
        canvas->Update();
      } 

    }//end over events

  }//end loop over ntuples
  cout << numEvents_tracks << endl;
  cout << numEvents_clusters << endl;
  /*const double tot_eta = 1.6;  
  for(int i = 1; i < hNum2->GetNbinsX()+1; i++)
    {
      double dpt = hNum->GetBinWidth(i);
      double content = hNum->GetBinContent(i);
      double temp = content/((double)numEvents*dpt*tot_eta);
      hNum->SetBinContent(i, temp);
      double error = hNum->GetBinError(i);
      double tempErr = error/((double)numEvents*dpt*tot_eta);
      hNum->SetBinError(i, tempErr);

      dpt = hNum2->GetBinWidth(i);
      content = hNum2->GetBinContent(i);
      temp = content/((double)numEvents*dpt*tot_eta);
      hNum2->SetBinContent(i, temp);
      error = hNum2->GetBinError(i);
      tempErr = error/((double)numEvents*dpt*tot_eta);
      hNum2->SetBinError(i, tempErr);

      dpt = hDen->GetBinWidth(i);
      content = hDen->GetBinContent(i);
      temp = content/((double)numEvents*dpt*tot_eta);
      hDen->SetBinContent(i, temp);
      error = hDen->GetBinError(i);
      tempErr = error/((double)numEvents*dpt*tot_eta);
      hDen->SetBinError(i, tempErr);
      }*/

  TFile* fout_cluster = new TFile("PhotonEfficiency_1Meve.root","RECREATE");

  TGraphAsymmErrors* eff_cluster = new TGraphAsymmErrors(&h_Num, &h_Den);
  eff_cluster->Write("Efficiency");
  h_Den.Write("hTruth");
  h_Num.Write("hRecoEmbed");
  h_Num2.Write("hReco");
  h_Correlation.Write();
  fout_cluster->Close();
 

  TFile* fout_track = new TFile("TrackEfficiency_20GeV_start1GeV_1Meve.root","RECREATE");

  TGraphAsymmErrors* eff = new TGraphAsymmErrors(&hNum, &hDen);
  eff->SetTitle("; p_{T}^{true} ; #epsilon");
  eff->Write("Efficiency");
  hDen.Write("hTruth");
  hNum.Write("hRecoEmbed");
  hNum2.Write("hReco");
  hCorrelation.Write();
  hFake.Write("hFake");
  hFake.Divide(&hReco);
  hFake.Write("FakeRate");
  fout_track->Close();
  
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}//end main
