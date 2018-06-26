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
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>


#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

using namespace std;

double SetPthatWeights(TString filename, double Xsection, double ntrial)
{

  TString MC = filename(filename.Last('/')+1,2);
  TString MCname = "";
  if(MC == "16")
    MCname = filename(filename.Last('/')+1,12);
  else
      MCname = filename(filename.Last('/')+1,13);
  
  double pthatWeight = 1.0;
  
  //17g6a3 weights
  if(MCname == "17g6a3_pthat1")
    return 4.47e-11;
  if(MCname == "17g6a3_pthat2")
    return 9.83e-11;
  if(MCname == "17g6a3_pthat3")
    return 1.04e-10;
  if(MCname == "17g6a3_pthat4")
    return 1.01e-11;
  
  //16c3c weights
  if(MCname == "16c3c_pthat1")
    return 3.941701e-03;
  if(MCname == "16c3c_pthat2")
    return 2.001984e-03;
  if(MCname == "16c3c_pthat3")
    return 9.862765e-04 ;
  if(MCname == "16c3c_pthat4")
    return 9.862765e-04 ;

  //18b10ab
  if(MC == "18")
    return Xsection/ntrial;
 
  return 0.0;
  
}

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
  
  const int nbinsphi = 80;
  Double_t phibins[nbinsphi+1] = {};
  double phimin = -1.0*TMath::Pi();
  double phimax = 1.0*TMath::Pi();
  double phistep = (phimax-phimin)/nbinsphi;
  for(int i=0; i<nbinsphi+1; i++){
    phibins[i] = phimin + i*phistep;
  }
  
  const int nbinstrack = 15;
  Double_t trackbins[nbinstrack+1] = {};
  //Double_t trackbins[nbinstrack+1] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};//pPB
  //Double_t trackbins[nbinstrack+1] = {0.5, 1.0,2.0,3.0,4.0,5.0,6.0, 8.0, 10.0, 13.0, 20.0};//pp bining
  double ptmin = 1;
  double ptmax = 16;
  double ptstep = (ptmax-ptmin)/nbinstrack;
  for(int i = 0; i < nbinstrack+1; i++){
    trackbins[i] = ptmin + i*ptstep;
  }
  
  TH1D h_Den("h_Den", "", 20, 10, 30);    
  TH1D h_Num("h_Num", "", 20, 10, 30);    
  TH1D h_Reco("h_Reco", "", 20, 10, 30);    
  TH2D h_Correlation("h_Correlation", "", 20, 10, 30, 20, 10, 30);
  TH2D h_Num2D("h_Num2D","", nbinsphi, phibins, nbinseta, etabins);
  TH2D h_Den2D("h_Den2D","", nbinsphi, phibins, nbinseta, etabins);

  h_Den.Sumw2();
  h_Num.Sumw2();
  h_Reco.Sumw2();

  h_Den.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  h_Num.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  h_Reco.SetTitle("reco photons filled with pt reco; p_{T}^{Reco,pt} [GeV]; entries");
  h_Correlation.SetTitle("; True p_{T} [GeV]; Reconstructed p_{T} [GeV]");
  h_Num2D.SetTitle(";#phi_{true}; #eta_{true}");
  h_Den2D.SetTitle(";#phi_{true}; #eta_{true}");

  TH1F hTruth("hTruth", "", nbinstrack, trackbins);
  TH1F hRecoTruth("hRecoTruth","", nbinstrack, trackbins);
  TH1F hRecof("hRecof","",nbinstrack, trackbins);
  TH1F hReco("hReco","", nbinstrack,trackbins);
  TH1F hFake("hFake", "", nbinstrack, trackbins);
  TH1F hTrackCut("hTrackCut", "", 10, -0.5, 9.5);
  TH1F hTrackQuality("hTrackQuality", "", 20, -0.5, 19.5);
  
  TH1F hTruth_eta("hTruth_eta","", nbinseta, etabins);
  TH1F hRecoTruth_eta("hRecoTruth_eta","", nbinseta, etabins);
  TH1F hReco_eta("hReco_eta","", nbinseta, etabins);
  TH1F hTruth_phi("hTruth_phi","", nbinsphi, phibins);
  TH1F hRecoTruth_phi("hRecoTruth_phi","", nbinsphi, phibins);
  TH1F hReco_phi("hReco_phi","", nbinsphi, phibins);
  
  TH2D hRecoTruth2D("hRecoTruth2D","", nbinsphi, phibins, nbinseta, etabins);
  TH2D hReco2D("hReco2D","", nbinsphi, phibins, nbinseta, etabins);
  TH2D hTruth2D("hTruth2D","", nbinsphi, phibins, nbinseta, etabins);
  TH2F hCorrelation("hCorrelation", "", nbinstrack, trackbins, nbinstrack, trackbins);
  TH2F hCorrelation_cor("hCorrelation_cor", "", nbinstrack, trackbins, nbinstrack, trackbins);
  TH2F hRes_Pt("hRes_Pt", "", 200, 0, 10.0, 80, -50, 50);
  
  TH1F hZvertex("hZvertex","",60, -30, 30);

  hTruth.Sumw2();
  hRecoTruth.Sumw2();
  hReco.Sumw2();
  hFake.Sumw2();
  hRecof.Sumw2();
  hTruth_eta.Sumw2();
  hRecoTruth_eta.Sumw2();
  hReco_eta.Sumw2();
  hTruth_phi.Sumw2();
  hRecoTruth_phi.Sumw2();
  hReco_phi.Sumw2();
  hCorrelation.Sumw2();
  hCorrelation_cor.Sumw2();
  hZvertex.Sumw2();


  hTruth.SetTitle("; p_{T}^{true} [GeV/c]; entries");
  hRecoTruth.SetTitle("; p_{T}^{Reco,Embed} [GeV/c]; entries");
  hReco.SetTitle("; p_{T}^{Reco} [GeV/c]; entries");
  hFake.SetTitle("; p_{T}^{reco} [GeV/c]; Fake Rate");
  hCorrelation.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
  hCorrelation_cor.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
  hRecoTruth2D.SetTitle(";#phi;#eta");
  hReco2D.SetTitle(";#phi;#eta");
  hTruth2D.SetTitle(";#phi;#eta");
  hZvertex.SetTitle(";Z_{v} [cm]; counts");
  
  const int TrackBit = 16;//3 for TPC+ITS, 16 for ITS only 
  int numEvents, numEvents_tracks, numEvents_clusters;
  numEvents = numEvents_tracks = numEvents_clusters = 0;
  double sumWeight = 0.0;
  
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
    Long64_t totEvents = _tree_event->GetEntries();
    Long64_t restrictEvents = 100000;
    Long64_t numEntries = TMath::Min(totEvents,restrictEvents);
    cout << numEntries << endl;
    // Loop over events
    for(Long64_t ievent = 0; ievent < numEntries ; ievent++){
    // for(Long64_t ievent = 0; ievent < 10000 ; ievent++){
      _tree_event->GetEntry(ievent);
      
      bool eventChange = true;
      
      //Selecting pthat weights
      double weight = (SetPthatWeights((TString)argv[iarg], (double)eg_cross_section, (double)eg_ntrial))/totEvents;
      if(ievent%10000 == 0)
	cout << weight << endl;
      sumWeight += weight;

      //event selection
      //if(not( TMath::Abs(primary_vertex[2])<10.0)) continue; //vertex z position
      hZvertex.Fill(primary_vertex[2]);
      numEvents++;

      //loop over tracks
      for (int n = 0;  n< ntrack; n++){

	//Track Cuts
	hTrackCut.Fill(0);
	//cout << (int)track_quality[n] << endl;
	hTrackQuality.Fill((int)track_quality[n]);
	if((track_quality[n]&TrackBit)==0) continue; hTrackCut.Fill(1);//track quality cut
	if(TMath::Abs(track_eta[n])> maxEta) continue; hTrackCut.Fill(2);//eta cut
	if(track_pt[n] < 0.15) continue; hTrackCut.Fill(3);//pt cut
	if(track_its_chi_square[n]>36.0) continue; hTrackCut.Fill(4);//its cluster chi^2 cut
	if(TrackBit == 16)
	  {
	    
	    if(track_its_ncluster[n] < 4) continue; 
	    hTrackCut.Fill(6);//its cluster cut
	  }
	if(TMath::Abs(track_dca_xy[n]) > 2.4) continue; hTrackCut.Fill(7);//distance of closest approach cut
	if(TMath::Abs(track_dca_z[n]) > 3.2) continue; hTrackCut.Fill(8);//distance of closest approach cut
	
	hRecof.Fill(track_pt[n],weight);
	
	if(track_pt[n] > 1){
	  hReco_eta.Fill(track_eta[n]);
	  hReco_phi.Fill(track_phi[n]);
	  hReco2D.Fill(track_phi[n], track_eta[n]);
	}

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
	  
	  hRecoTruth.Fill(mc_truth_pt[index],weight);
	  hReco.Fill(track_pt[n],weight);
	  hCorrelation.Fill(mc_truth_pt[index], track_pt[n], weight);
	  hCorrelation_cor.Fill(mc_truth_pt[index], track_pt[n], weight);
	  hRes_Pt.Fill(mc_truth_pt[index], 100*(track_pt[n]-mc_truth_pt[index])/(mc_truth_pt[index]),weight);
	  hRecoTruth_eta.Fill(mc_truth_eta[index]);
	  hRecoTruth_phi.Fill(mc_truth_phi[index]);
	  hRecoTruth2D.Fill(mc_truth_phi[index], mc_truth_eta[index]);
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
        
	
	hTruth.Fill(mc_truth_pt[nTru],weight);
	hTruth_eta.Fill(mc_truth_eta[nTru]);
	hTruth_phi.Fill(mc_truth_phi[nTru]);
	hTruth2D.Fill(mc_truth_phi[nTru], mc_truth_eta[nTru]);
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
        Float_t truth_eta = -999.0;
        Float_t truth_phi = -999.0;
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
	  truth_phi    = mc_truth_phi[index];
	  truth_eta    = mc_truth_eta[index];
	}//end loop over indices
	
	if(isTruePhoton){
	  //fill in this histogram only photons that can be traced to a generated non-decay photon.	
	  h_Num.Fill(truth_pt,weight);
	  h_Reco.Fill(cluster_pt[n],weight); 
	  h_Correlation.Fill(truth_pt, cluster_pt[n],weight);
	  h_Num2D.Fill(truth_phi, truth_eta);
	  if (eventChange) {numEvents_clusters++; eventChange = false;}
	} 
      }//end loop on clusters
       

      


      //loop over truth particles
      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
        //if(mc_truth_pt[nmc]<5.0) continue;
	if(mc_truth_pdg_code[nmc]==22 && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==22){
	  //std::cout << mc_truth_pt[nmc] << "phi " << mc_truth_phi[nmc] << " eta " << mc_truth_eta[nmc] << 
	  //  " code: " << mc_truth_pdg_code[nmc] << " status " << int(mc_truth_status[nmc]) << " parentpdg " << mc_truth_first_parent_pdg_code[nmc] << std::endl;    
	  h_Den.Fill(mc_truth_pt[nmc],weight);
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
      hs->Add(&h_Reco);
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
  cout << sumWeight << endl;

  TH1D eventSelection("eventSelection","", 8, -0.5, 7.5);
  

  const double tot_eta = 1.6;  
  for(int i = 1; i < hReco.GetNbinsX()+1; i++)
    {
      double dpt, content, temp, error, tempErr;
      dpt = content = temp = error = tempErr = 0.0;
      
      dpt = hRecoTruth.GetBinWidth(i);
      content = hRecoTruth.GetBinContent(i);
      temp = content/dpt;
      hRecoTruth.SetBinContent(i, temp);
      error = hRecoTruth.GetBinError(i);
      tempErr = error/dpt;
      hRecoTruth.SetBinError(i, tempErr);

      dpt = hReco.GetBinWidth(i);
      content = hReco.GetBinContent(i);
      temp = content/dpt;
      hReco.SetBinContent(i, temp);
      error = hReco.GetBinError(i);
      tempErr = error/dpt;
      hReco.SetBinError(i, tempErr);

      dpt = hTruth.GetBinWidth(i);
      content = hTruth.GetBinContent(i);
      temp = content/dpt;
      hTruth.SetBinContent(i, temp);
      error = hTruth.GetBinError(i);
      tempErr = error/dpt;
      hTruth.SetBinError(i, tempErr);

    }

  for(int x = 1; x < hCorrelation_cor.GetNbinsX()+1; x++)
    {
      for(int y = 1; y < hCorrelation_cor.GetNbinsY()+1; y++)
	{
	  double dx , dy, content, temp, contentErr, tempErr;
	  dx = dy = content = temp = contentErr = tempErr = 0.0;
	  dx = hCorrelation_cor.GetXaxis()->GetBinWidth(x);
	  dy = hCorrelation_cor.GetYaxis()->GetBinWidth(y);
	  content = hCorrelation_cor.GetBinContent(x, y);
	  contentErr = hCorrelation_cor.GetBinError(x, y);

	  temp = content/(dx*dy);
	  hCorrelation_cor.SetBinContent(x, y, temp);
	  tempErr = tempErr/(dx*dy);
	  hCorrelation_cor.SetBinError(x, y, tempErr);

	}
    }
  hCorrelation_cor.Scale(1/sumWeight);

  hTrackCut.GetXaxis()->SetBinLabel(1,"All Tracks");
  hTrackCut.GetXaxis()->SetBinLabel(2,"Track quality cut");
  hTrackCut.GetXaxis()->SetBinLabel(3,"Track #eta cut");
  hTrackCut.GetXaxis()->SetBinLabel(4,"pt cut");
  hTrackCut.GetXaxis()->SetBinLabel(5,"Trigger cut");
  hTrackCut.GetXaxis()->SetBinLabel(6,"ITS nCluster cut");
  hTrackCut.GetXaxis()->SetBinLabel(7,"ITS #chi^{2} cut");
  hTrackCut.GetXaxis()->SetBinLabel(8,"DCAr cut");
  hTrackCut.GetXaxis()->SetBinLabel(9,"DCAz cut");

  //TF1* gausfit = new TF1("gaus","gaus", -25,25);
  //gausfit->SetLineColor(kRed);
  //TGraphErrors* g_mean = new TGraphErrors();
  //TGraphErrors* g_sigma = new TGraphErrors();
  
  //auto c1 = new TCanvas();
   
  /*//Study of the ITS-only track pT resolution
  const Double_t bins[10] = {  0.1,          0.16681005,   0.27825594,   0.46415888 ,  0.77426368, 1.29154967,   2.15443469,   3.59381366,   5.9948425,   10.0};

  int nbins = 9;
  for(int i=0; i<nbins; i++){
    double minpt = bins[i];
    double maxpt = bins[i+1];
    double binwidth = maxpt-minpt;
    int minbin =  hRes_Pt.GetXaxis()->FindBin(minpt);
    int maxbin =  hRes_Pt.GetXaxis()->FindBin(maxpt);
    TH1D* h1 = hRes_Pt.ProjectionY("h", minbin, maxbin);
    
    h1->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    h1->Draw();  
    h1->GetYaxis()->SetNdivisions(5);
    h1->GetXaxis()->SetNdivisions(5);
    h1->GetYaxis()->SetTitle("counts");
    h1->Fit(gausfit,"R");
    h1->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    gausfit->Draw("same");
    //myText(0.18, 0.8, kBlack, Form("%2.1f < p_{T}^{truth} < %2.1f GeV", minpt, maxpt));
    //myText(0.18, 0.74, kRed, Form("#mu = %2.1f [%]", gausfit->GetParameter(1)));
    //myText(0.18, 0.68, kRed, Form("#sigma = %2.1f [%]", gausfit->GetParameter(2))); 
    g_sigma->SetPoint(g_sigma->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(2));     
    g_sigma->SetPointError(g_sigma->GetN()-1, binwidth/2.0, gausfit->GetParError(2));
    g_mean->SetPoint(g_mean->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(1));
    g_mean->SetPointError(g_mean->GetN()-1, binwidth/2.0, gausfit->GetParError(1));
    
    }*/
  
  //g_sigma->SetTitle("Relative resolution vs p_{T} ; p_{T}^{true} [GeV]; #sigma(p_{T})/p_{T} [%]"); 
  //g_mean->SetTitle("; p_{T}^{true} [GeV]; Relative bias [%]");
  
  
  bool makeClusterFile = false;
  bool makeTrackFile = true;
  
  if(makeClusterFile)
    {
      TFile* fout_cluster = new TFile("PhotonEfficiency_pp_30GeV_2M_etaPhi.root","RECREATE");
      
      TGraphAsymmErrors* eff_cluster = new TGraphAsymmErrors(&h_Num, &h_Den);
      eff_cluster->Write("Efficiency");
      h_Den.Write("hTruth");
      h_Num.Write("hRecoEmbed");
      h_Reco.Write("hReco");
      h_Correlation.Write();
      h_Num2D.Write("hEtaPhi");
      fout_cluster->Close();
    }

  if(makeTrackFile)
    {
      TFile* fout_track = new TFile(Form("TrackOutput/17g6a3_%i_1GeV16GeV_500K_eventDiv_4layers_new.root",TrackBit),"RECREATE");
      
      TGraphAsymmErrors* eff = new TGraphAsymmErrors(&hRecoTruth, &hTruth);
      eff->SetTitle("; p_{T}^{true} ; #epsilon");
      eff->Write("Efficiency");
      //g_sigma->Write("g_sigma");
      //g_mean->Write("g_mean");

      hTruth.Write("hTruth");
      hRecoTruth.Write("hRecoTruth");
      hReco.Write("hReco");
      hTruth_eta.Write("hTruth_eta");
      hRecoTruth_eta.Write("hRecoEmbed_eta");
      hReco_eta.Write("hReco_eta");
      hTruth_phi.Write("hTruth_phi");
      hRecoTruth_phi.Write("hRecoEmbed_phi");
      hReco_phi.Write("hReco_phi");
      hFake.Write("hFake");
      hFake.Divide(&hRecof);
      hFake.Write("FakeRate");
      hZvertex.Write("hZvertex");
      
      hCorrelation.Write("hCorrelation");
      hCorrelation_cor.Write("hCorrelation_cor");
      hRecoTruth2D.Write("hRecoTruth_phiEta");
      hReco2D.Write("hReco_phiEta");
      hTruth2D.Write("hTruth_phiEta");

      hTrackQuality.Write();
      hTrackCut.Write();

      fout_track->Close();
    }
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}//end main

