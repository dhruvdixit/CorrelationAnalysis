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

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");
  
  //Config File
  FILE* config = fopen("Corr_config.yaml", "r");
  double DNN_min = 0;
  double DNN_max = 0;
  double Lambda0_cut = 0; 
  double pT_min = 0;
  double pT_max = 0;
  double Eta_max = 0;
  double Cluster_min = 0;
  float Cluster_DtoBad = 0;
  UChar_t Cluster_NLocal_Max = 0;
  double EcrossoverE_min = 0;
  int Track_Cut_Bit = 0;
  double iso_max = 0;
  double noniso_min = 0;
  double noniso_max = 0;
  double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  int n_eta_bins = 0;
  int n_phi_bins = 0;  
  std::string shower_shape = "DNN";

  // Zt bins
  //FIXME: Will have to likely set nztbins first, then initialize array
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
   
  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;


  // Loop through config file
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
      if (line[0] == '#') continue;
      
      char key[MAX_INPUT_LENGTH];
      char dummy[MAX_INPUT_LENGTH];
      char value[MAX_INPUT_LENGTH];
      
      // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
      key[0] = '\0';
      value[0] = '\0';
      sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
      
      //Read Config File: Detect Keys 
      if (strcmp(key, "DNN_min") == 0) {
          DNN_min = atof(value);
          std::cout << "DNN_min: " << DNN_min << std::endl; }

      else if (strcmp(key, "DNN_max") == 0) {
          DNN_max = atof(value);
          std::cout << "DNN_max: " << DNN_max << std::endl; }

      else if (strcmp(key, "Lambda0_cut") == 0){
	Lambda0_cut = atof(value);
	std::cout << "Lambda0_cut: " << Lambda0_cut << std::endl;}


      else if (strcmp(key, "pT_min") == 0) {
          pT_min = atof(value);
          std::cout << "pT_min: " << pT_min << std::endl; }

      else if (strcmp(key, "pT_max") == 0) {
          pT_max = atof(value);
          std::cout << "pT_max: " << pT_max << std::endl; }

      else if (strcmp(key, "Eta_max") == 0) {
          Eta_max = atof(value);
          std::cout << "Eta_max: " << Eta_max << std::endl;
      }
      else if (strcmp(key, "Cluster_min") == 0) {
          Cluster_min = atof(value);
          std::cout << "Cluster_min: " << Cluster_min << std::endl; }

      else if (strcmp(key, "Cluster_dist_to_bad_channel") == 0){
	Cluster_DtoBad = atof(value);
	std::cout << "Cluster_DtoBad: "<< Cluster_DtoBad << std::endl;}

      else if (strcmp(key, "Cluster_N_Local_Maxima") == 0){
        Cluster_NLocal_Max = atof(value);
	std::cout << "Cluster_NLocal_Max: "<< Cluster_NLocal_Max << std::endl;}

      else if (strcmp(key, "EcrossoverE_min") == 0) {
          EcrossoverE_min = atof(value);
          std::cout << "EcrossoverE_min; " << EcrossoverE_min << std::endl; }

      else if (strcmp(key, "iso_max") == 0) {
          iso_max = atof(value);
          std::cout << "iso_max: " << iso_max << std::endl; }

      else if (strcmp(key, "noniso_min") == 0) {
          noniso_min = atof(value);
          std::cout << "noniso_min: " << noniso_min << std::endl; }

      else if (strcmp(key, "noniso_max") == 0) {
          noniso_max = atof(value);
          std::cout << "noniso_max: " << noniso_max << std::endl; }

      else if (strcmp(key, "deta_max") == 0) {
          deta_max = atof(value);
          std::cout << "deta_max: " << deta_max << std::endl; }

      else if (strcmp(key, "N_Phi_Bins") == 0) {
	n_phi_bins = atoi(value);
	std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

      else if (strcmp(key, "N_Eta_Bins") == 0) {
        n_eta_bins = atoi(value);
	std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

      else if (strcmp(key, "Track_Cut_Bit") == 0) {
          Track_Cut_Bit = atoi(value);
          std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl; }

      else if (strcmp(key, "Zt_bins") == 0) {
          nztbins = -1;
          for (const char *v = value; *v != ']';) {
              while (*v != ']' && !isdigit(*v)) v++;
	      nztbins++;
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          ztbins = new float[nztbins + 1];
          int i = 0;
          for (const char *v = value; *v != ']' ;) {
              while (*v != ']' && !isdigit(*v)) v++;
              ztbins[i] = atof(v);
              i++;              
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
          for (int i = 0; i <= nztbins; i++)
	    std::cout << ztbins[i] << ", ";
          std::cout << "}\n"; 
      }

      else if (strcmp(key, "Pt_bins") == 0) {
          nptbins = -1;
          for (const char *v = value; *v != ']';) {
            while (*v != ']' && !isdigit(*v)) v++;
	    nptbins++;
	    while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
	
          ptbins = new float[nptbins + 1];
	  int i = 0;
	  for (const char *v = value; *v != ']' ;) {
	     while (*v != ']' && !isdigit(*v))  v++;
	     ptbins[i] = atof(v);
	     i++;
	     while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
	
	  std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
          for (int i = 0; i <= nptbins; i++)
	    std::cout << ptbins[i] << ", ";
	  std::cout << "}\n";
      }

      else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
          if (strcmp(value, "cluster_iso_tpc_04") == 0){
              determiner = CLUSTER_ISO_TPC_04;
              std::cout << "Isolation Variable: cluster_iso_tpc_04" << std::endl; }

          else if (strcmp(value, "cluster_iso_its_04") == 0){
              determiner = CLUSTER_ISO_ITS_04;
              std::cout << "Isolation Variable: cluster_iso_its_04" << std::endl; }

          else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_TPC_04_02;
              std::cout << "Isolation Variable: cluster_frixione_tpc_04_02" << std::endl; }

          else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_ITS_04_02;
              std::cout << "Isolation Variable: cluster_frixione_its_04_02" << std::endl; }

          else {
              std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
              exit(EXIT_FAILURE); }
      }

      else if (strcmp(key, "Shower_Shape") == 0){
	  shower_shape = value;
	  std::cout<<"Shower Shape: "<<shower_shape.data()<<std::endl;
	  //if (strcmp(shower_shape.data(),"Lambda")== 0) std::cout<<"test worked"<<std::endl;
      }

      else std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
  
  }
  //end Config Loop

  fclose(config);
  
  for (int i = 0; i <= nztbins; i++)
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  for (int i = 0; i <= nptbins; i++)
    std::cout << "pt bound: " << ptbins[i] << std::endl;


  //HISTOGRAMS
  TCanvas canvas("canvas", "");

  TH1D h_ntrig("h_ntrig", "", 2, -0.5,1.0);

  TH2D* Signal_pT_Dist = new TH2D("Signal_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",24,10,16,5,-0.5,2);
  TH2D* BKGD_pT_Dist = new TH2D("BKGD_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",24,10,16,5,-0.5,2);

  TH2D* Corr[nztbins*nptbins];
  TH2D* IsoCorr[nztbins*nptbins];
  TH2D* BKGD_IsoCorr[nztbins*nptbins];

  TH1D* H_Signal_Triggers[nptbins];
  TH1D* H_BKGD_Triggers[nptbins];
  TH1D* Triggers[nptbins];

  float N_Signal_Triggers = 0;
  float N_BKGD_Triggers = 0;
  
  //FIXME: Add to config file

    for (int ipt = 0; ipt <nptbins; ipt++) {
      H_Signal_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Photon Triggers", 2, -0.5,1.0);

      H_BKGD_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);

      Triggers[ipt] = new TH1D(
      Form("N_Triggers_pT%1.0f_%1.0f",ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);

      for (int izt = 0; izt<nztbins; izt++){

      Corr[izt+ipt*nztbins] = new TH2D(Form("Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [all] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      Corr[izt+ipt*nztbins]->Sumw2();
      Corr[izt+ipt*nztbins]->SetMinimum(0.);

      IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [Iso] Correlation", n_phi_bins,0,M_PI,n_eta_bins, -1.4, 1.4);

      IsoCorr[izt+ipt*nztbins]->Sumw2();
      IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);

      BKGD_IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [AntiIso] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      BKGD_IsoCorr[izt+ipt*nztbins]->Sumw2();
      BKGD_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);

    }//zt bins
  }//pt bins                           
  
    //for (int iarg = 1; iarg < argc; iarg++) {
    int iarg = 1;
    TString root_file = (TString)argv[iarg];
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;

    TFile *file = TFile::Open(root_file);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
    
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " fail " << std::endl;
	exit(EXIT_FAILURE);
      }
    }  

    //Tracks
    Double_t primary_vertex[3];
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    Float_t track_eta_emcal[NTRACK_MAX];
    Float_t track_phi_emcal[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
    UChar_t track_its_ncluster[NTRACK_MAX];
    Float_t track_its_chi_square[NTRACK_MAX];
    Float_t track_dca_xy[NTRACK_MAX];

    //Clusters
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
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];   
    Float_t cell_e[17664];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];

    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];

    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
    //Float_t eg_cross_section;
    //Int_t   eg_ntrial;
     
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);
  
    //track Addresses
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_eta_emcal", track_eta_emcal);
    _tree_event->SetBranchAddress("track_phi_emcal", track_phi_emcal);
    _tree_event->SetBranchAddress("track_quality", track_quality);
    _tree_event->SetBranchAddress("track_its_ncluster", &track_its_ncluster);
    _tree_event->SetBranchAddress("track_its_chi_square", &track_its_chi_square);
    _tree_event->SetBranchAddress("track_dca_xy", &track_dca_xy);

    //Cluster Addresses
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);

    //_tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    //_tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);

    Long64_t nentries = _tree_event->GetEntries();         
    std::cout << " Total Number of entries in TTree: " << nentries << std::endl;

    for(Long64_t ievent = 0; ievent < nentries ; ievent++){     
      //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      //double weight = (double)eg_cross_section/(double)eg_ntrial;

      for (ULong64_t n = 0; n < ncluster; n++) {
	if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue;   //select pt of photons
	if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue;              //cut edges of detector
	if( not(cluster_ncell[n]>Cluster_min)) continue;                    //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
	if( not(cluster_distance_to_bad_channel[n]>=Cluster_DtoBad)) continue; //removes clusters near bad channels
	//	if( not(cluster_nlocal_maxima[n] < Cluster_NLocal_Max)) continue; //require to have at most 2 local maxima.
	if( not(cluster_nlocal_maxima[n] < 3)) continue; //require to have at most 2 local maxima.

	float isolation;
	if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
	else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
	else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
	else isolation = cluster_frixione_its_04_02[n];
	
	Bool_t Signal = false;
	Bool_t Background = false;

	if (strcmp(shower_shape.data(),"Lambda")== 0) {
	  if ((cluster_lambda_square[n][0] < Lambda0_cut))
	    Signal = true;
	  
	  if ((cluster_lambda_square[n][0] > Lambda0_cut))
	    Background = true;
	}
	
	else if (strcmp(shower_shape.data(),"DNN")==0){
	  if ((cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1]<DNN_max))
	    Signal = true;
	  if (cluster_lambda_square[n][0] > Lambda0_cut)
	    Background = true;
	}
	

	//count triggers	
	if (isolation<iso_max){	
	  if (Signal){  
	    N_Signal_Triggers += 1;
	    for (double Iso_bin = -0.5; Iso_bin < iso_max; Iso_bin += 0.5)
	      if (isolation > Iso_bin && isolation < Iso_bin+0.5) Signal_pT_Dist->Fill(cluster_pt[n],isolation);
	    for (int ipt = 0; ipt < nptbins; ipt++)
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1])
		H_Signal_Triggers[ipt]->Fill(1);
	  }
	  //Low DNN Trigger BKGD
	  if (Background){
	    h_ntrig.Fill(0.5);
	    N_BKGD_Triggers += 1;
	    for (double Iso_bin = -0.5; Iso_bin < iso_max; Iso_bin += 0.5)
	      if (isolation > Iso_bin && isolation < Iso_bin+0.5) 
		BKGD_pT_Dist->Fill(cluster_pt[n],isolation);
	    for (int ipt = 0; ipt < nptbins; ipt++)
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]) 
		H_BKGD_Triggers[ipt]->Fill(1); 
	  }
	}
	//no iso, no dnn
	for (int ipt = 0; ipt < nptbins; ipt++)
	  Triggers[ipt]->Fill(1);

	//Track Loop
	for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {            
 	  if(track_pt[itrack] < 0.5) continue; //500 MeV Tracks
 	  if(track_pt[itrack] > 30) continue;
 	  if((track_quality[itrack]&Track_Cut_Bit)==0) continue; //select only tracks that pass selection 
	  if(abs(track_eta[itrack]) > 0.8) continue;
	  if( not(track_its_ncluster[itrack]>4)) continue;
	  if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
	  if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;

	  //Electron Veto for associated tracks outside of isolation cone
	  double dRmin = 0.02;
	  bool Track_HasMatch = false;
	  for (ULong64_t c = 0; c < ncluster; c++){
	    Float_t deta =  cluster_eta[c]-track_eta_emcal[itrack];
	    Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[c]-track_phi_emcal[itrack])/TMath::Pi();
	    float dR = sqrt(dphi*dphi + deta*deta);
	    if (dR < dRmin) {
	      Track_HasMatch = true;
	      break;
	    }
	  }
 	  if (Track_HasMatch) continue;

	  //Observables:
	  Double_t zt = track_pt[itrack]/cluster_pt[n];
	  Float_t DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - track_phi[itrack]));
	  Float_t DeltaEta = cluster_eta[n] - track_eta[itrack];
	  if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue; //Match Mixing Cut

	  for (int ipt = 0; ipt < nptbins; ipt++){
	    if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
	      for(int izt = 0; izt<nztbins ; izt++){
		if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
		  
		  //2 DNN Regions
		  if (isolation<iso_max){
		    //if (cluster_s_nphoton[n][1] > DNN_min && cluster_s_nphoton[n][1] < DNN_max)
		    //if ((cluster_lambda_square[n][0] < Lambda0_cut))   
		    if (Signal)
		      IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);

		    //if (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < 0.3) //sel deep photons       
		    //if ((cluster_lambda_square[n][0] > Lambda0_cut))                       
		    if (Background)
		      BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		  }		  
		  //No Iso, no DNN
		  Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);

		}//if in zt bin                                                                                               
	      } // end loop over zt bins                                                                         
	    }//end if in pt bin                                                                                                  
	  }//end pt loop bin   
	}//end loop over tracks
	//fprintf(stderr,"\n"); 
      }//end loop on clusters.
    } //end loop over events  
    //}//end loop over samples

  // Write to fout
  
  //TFile* fout = new TFile(Form("fout_Corr_config%s.root", opened_files.c_str()),"RECREATE");

  size_t lastindex = std::string(root_file).find_last_of(".");
  std::string rawname = std::string(root_file).substr(0, lastindex);
  fprintf(stderr,"%s: %d: Creating new file",__FILE__,__LINE__);

  TFile* fout;
  if (strcmp(shower_shape.data(),"Lambda")== 0) 
    fout = new TFile(Form("%s_SE_L0_Correlation.root",rawname.data()),"RECREATE");
  else if (strcmp(shower_shape.data(),"DNN")== 0)
    fout = new TFile(Form("%s_SE_NN_Correlation.root",rawname.data()),"RECREATE");
  else
    fout = new TFile(Form("%s_SE_Correlation.root",rawname.data()),"RECREATE");

  // TFile* fout = new TFile("Same_Event_Correlation_13defv1.root","RECREATE");
  h_ntrig.Write("ntriggers");
  std::cout<<"Clusters Passed Iosalation: "<<N_Signal_Triggers<<std::endl;
			  
  Signal_pT_Dist->Write();
  BKGD_pT_Dist->Write();

  for (int ipt = 0; ipt<nptbins; ipt++)
    H_Signal_Triggers[ipt]->Write();
  
  for (int ipt = 0; ipt<nptbins; ipt++)
    H_BKGD_Triggers[ipt]->Write();
    
  for (int ipt = 0; ipt < nptbins; ipt++)
    Triggers[ipt]->Write();

  for (int ipt = 0; ipt<nptbins; ipt++){
    for (int izt = 0; izt<nztbins; izt++)
      Corr[izt+ipt*nztbins]->Write();
    for (int izt = 0; izt<nztbins; izt++)
      IsoCorr[izt+ipt*nztbins]->Write();
    for (int izt = 0; izt<nztbins; izt++)
      BKGD_IsoCorr[izt+ipt*nztbins]->Write();
  }
  //Seperate zt loops for easier file reading
  fout->Close();     
  file->Close();  
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
