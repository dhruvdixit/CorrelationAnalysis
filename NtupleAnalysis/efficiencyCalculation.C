#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraphAsymmErrors.h"
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TVector2.h"
//#include "/root/atlasstyle-00-03-05/AtlasStyle.h"
//#include "/root/atlasstyle-00-03-05/AtlasStyle.C"
//#include "/root/atlasstyle-00-03-05/AtlasUtils.h"
//#include "/root/atlasstyle-00-03-05/AtlasUtils.C"
//#include "/root/atlasstyle-00-03-05/AtlasLabels.h"
//#include "/root/atlasstyle-00-03-05/AtlasLabels.C"
#include "TDatabasePDG.h"
#include "TEfficiency.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <bitset>

void Run(const int TrackBit, TString address, bool isMC, bool hasAliDir, bool triggered){

  gStyle->SetOptStat(0);

  //SetAtlasStyle();
  TString filename = address(address.Last('/')+1,address.Last('.')-address.Last('/')-1);
  TFile* f;
  f = TFile::Open(Form("/project/projectdirs/alice/NTuples/%s",address.Data()),"READ");

  cout << Form("/project/projectdirs/alice/NTuples/%s",address.Data()) << endl;
  if(!f){
      printf("Error: cannot open ntuple.root");
      return;
  }
  TTree* tree;
  if(hasAliDir)
    tree = (TTree*)f->Get("AliAnalysisTaskNTGJ/_tree_event");
  else
    tree = (TTree*)f->Get("_tree_event");
  if(!tree){ printf("Error: cannot find tree"); }
  const Int_t kMax = 5000;
  
  Double_t primary_vertex[3];
  Bool_t is_pileup_from_spd_5_08;
  int run_number;

  unsigned int ntrack;//
  unsigned int nmc_truth;//
  float mc_truth_pt[kMax];//
  float mc_truth_eta[kMax];//
  float mc_truth_phi[kMax];//
  short mc_truth_pdg_code[kMax];//
  char mc_truth_charge[kMax];//
  float track_pt[kMax];//
  float track_eta[kMax];//
  float track_phi[kMax];//
  float track_dca_xy[kMax];//
  float track_dca_z[kMax];//
  float track_its_chi_square[kMax];//
  

  UInt_t ncluster;
  Float_t cluster_e[kMax];
  Float_t cluster_e_cross[kMax];
  Float_t cluster_pt[kMax];
  Float_t cluster_eta[kMax];
  Float_t cluster_phi[kMax];
  Float_t cluster_iso_tpc_04[kMax];
  Float_t cluster_iso_its_04[kMax];
  Float_t cluster_iso_04_truth[kMax];
  Float_t cluster_frixione_tpc_04_02[kMax];
  Float_t cluster_frixione_its_04_02[kMax];
  Float_t cluster_s_nphoton[kMax][4];
  UChar_t cluster_nlocal_maxima[kMax];
  Float_t cluster_distance_to_bad_channel[kMax];
  
  unsigned short cluster_mc_truth_index[kMax][32];
  Int_t cluster_ncell[kMax];
  UShort_t  cluster_cell_id_max[kMax];
  Float_t cluster_lambda_square[kMax][2];
  Float_t cell_e[17664];
    
  //Jets reco
  UInt_t njet_ak04its;
  Float_t jet_ak04its_pt_raw[kMax];
  Float_t jet_ak04its_eta_raw[kMax];
  Float_t jet_ak04its_phi_raw[kMax];
  UInt_t njet_ak04tpc;
  Float_t jet_ak04tpc_pt_raw[kMax];
  Float_t jet_ak04tpc_eta_raw[kMax];
  Float_t jet_ak04tpc_phi_raw[kMax];


  char track_charge[kMax];//
  unsigned short track_mc_truth_index[kMax];//
  UChar_t track_its_ncluster[kMax];//
  unsigned char track_quality[kMax];//
  ULong64_t trigger_mask[2];

  float eg_cross_section;//
  int eg_ntrial;//
  
  tree->SetBranchAddress("primary_vertex", primary_vertex);
  tree->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
  tree->SetBranchAddress("run_number", &run_number);

  tree->SetBranchAddress("nmc_truth",&nmc_truth);//
  tree->SetBranchAddress("mc_truth_pt",mc_truth_pt);//
  tree->SetBranchAddress("mc_truth_eta",mc_truth_eta);//
  tree->SetBranchAddress("mc_truth_phi",mc_truth_phi);//
  tree->SetBranchAddress("mc_truth_charge",mc_truth_charge);//
  tree->SetBranchAddress("mc_truth_pdg_code",mc_truth_pdg_code);//

  tree->SetBranchAddress("track_its_ncluster", track_its_ncluster);//
  tree->SetBranchAddress("ntrack", &ntrack);//
  tree->SetBranchAddress("track_pt",track_pt);//
  tree->SetBranchAddress("track_eta",track_eta);//
  tree->SetBranchAddress("track_phi",track_phi);//
  tree->SetBranchAddress("track_mc_truth_index", track_mc_truth_index);
  tree->SetBranchAddress("track_quality", track_quality);//
  tree->SetBranchAddress("track_dca_xy", track_dca_xy);//
  tree->SetBranchAddress("track_dca_z", track_dca_z);//
  tree->SetBranchAddress("track_its_chi_square", track_its_chi_square);//
  if(!isMC)
    tree->SetBranchAddress("trigger_mask", trigger_mask);
  if(isMC)
    {
      cout << "Setting Xsection and trail branch address" << endl;
      tree->SetBranchAddress("eg_cross_section",&eg_cross_section);//
      tree->SetBranchAddress("eg_ntrial",&eg_ntrial);//
    }

  tree->SetBranchAddress("ncluster", &ncluster);
  tree->SetBranchAddress("cluster_e", cluster_e);
  tree->SetBranchAddress("cluster_e_cross", cluster_e_cross);
  tree->SetBranchAddress("cluster_pt", cluster_pt); // here
  tree->SetBranchAddress("cluster_eta", cluster_eta);
  tree->SetBranchAddress("cluster_phi", cluster_phi);
  tree->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
  tree->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
  tree->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
  tree->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
  tree->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
  tree->SetBranchAddress("cluster_iso_04_truth", cluster_iso_04_truth);
  tree->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
  tree->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
  tree->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);        
  tree->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
  tree->SetBranchAddress("cluster_ncell", cluster_ncell);
  tree->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
  tree->SetBranchAddress("cell_e", cell_e);
  

  tree->SetBranchAddress("njet_ak04its", &njet_ak04its);
  tree->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
  tree->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
  tree->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi_raw);
  tree->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);
  tree->SetBranchAddress("jet_ak04tpc_pt_raw", jet_ak04tpc_pt_raw);
  tree->SetBranchAddress("jet_ak04tpc_eta_raw", jet_ak04tpc_eta_raw);
  tree->SetBranchAddress("jet_ak04tpc_phi", jet_ak04tpc_phi_raw);

  const Double_t bins[10] = {  0.1,          0.16681005,   0.27825594,   0.46415888 ,  0.77426368, 1.29154967,   2.15443469,   3.59381366,   5.9948425,   10.0};

  /*const Double_t bins_track[50] = {
    0.15,    0.169,   0.19,    0.214,   0.241,   0.271,   0.306,   0.344,   0.387,   0.436,   
    0.491,   0.553,   0.622,   0.701,   0.789,   0.888,   1,       1.126,   1.267,   1.427,
    1.606,   1.808,   2.036,   2.292,   2.581,   2.906,   3.271,   3.683,   4.147,   4.669,   
    5.256,   5.918,   6.663,   7.501,   8.445,   9.508,   10.705,  12.052,  13.569,  15.277,
    17.2,    19.365,  21.802,  24.546,  27.636,  31.114,  35.03,   39.439,  44.403,  49.992};//*/

  //const Double_t bins_track[13] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0};//pPb
  //const Double_t bins_track[34] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 35.0, 40.0, 45.0, 50.0};//pPb

  //const Double_t bins_track[9] = {1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0, 13.0};//pp

  /*const Double_t bins_track[37] = {
    0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 
    0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20,
    1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 
    3.60, 4.00, 5.00, 6.00, 8.00, 10.00, 13.00};//pp*/

  /*const Double_t bins_track[59] = {
    0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,  0.50,  0.55,  0.60, 
    0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,  1.00,  1.10,  1.20,
    1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90,  2.00,  2.20,  2.40,
    2.60,  2.80,  3.00,  3.20,  3.40,  3.60,  3.80,  4.00,  4.50,  5.00,
    5.50,  6.00,  6.50,  7.00,  8.00,  9.00,  10.00, 11.00, 12.00, 13.00,
    14.00, 15.00, 16.00, 18.00, 20.00, 22.00, 24.00, 26.00, 30.00};//pPb*/

  /*const Double_t bins_track[32] = {
    1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 
    2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 
    4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 9.00, 10.00,
    11.00,12.00};//pPb*/

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
  
  const int nbinstrack = 21;
  /*Double_t trackbins[nbinstrack+1] = {
    0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,  0.50,  0.55,  0.60, 
    0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,  1.00,  1.10,  1.20,
    1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90,  2.00,  2.20,  2.40,
    2.60,  2.80,  3.00,  3.20,  3.40,  3.60,  3.80,  4.00,  4.50,  5.00,
    5.50,  6.00,  6.50,  7.00,  8.00,  9.00,  10.00, 11.00, 12.00, 13.00,
    14.00, 15.00, 16.00, 18.00, 20.00, 22.00, 24.00, 26.00, 30.00};//nbinsbstrack = 58*/
  /*Double_t trackbins[nbinstrack+1] = {
    0.15,  0.30,  0.50,  0.70,  0.90,  1.00,  1.50,  2.50,  3.00,  3.50,  
    4.00,  4.50,  5.00,  5.50,  6.00,  6.50,  7.00,  8.00,  9.00,  10.00, 
    11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00, 20.00, 22.00, 24.00, 
    26.00, 30.00};//nbinsbstrack = 32*/
  Double_t trackbins[nbinstrack+1] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.00, 24.00, 26.00, 30.00};//pPB nbinstrack = 21*/
  //Double_t trackbins[nbinstrack+1] = {1.0,2.0,3.0,4.0,5.0,6.0, 8.0, 10.0, 13.0, 20.0};//pp bining nbinstrack = 9
  /*Double_t trackbins[nbinstrack+1] = {};//{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};
  double ptmin = 1;
  double ptmax = 16;
  double ptstep = (ptmax-ptmin)/nbinstrack;
  for(int i = 0; i < nbinstrack+1; i++){
    trackbins[i] = ptmin + i*ptstep;
  }//*/

  //tracks
  auto hCorrelation   = new TH2F("hCorrelation", "", nbinstrack, trackbins, nbinstrack, trackbins);
  auto hRes_Pt   = new TH2F("hRes_Pt", "", 200, 0, 10.0, 80, -50, 50);
  auto hDen  = new TH1F("hDen", "", nbinstrack, trackbins);
  auto hNum  = new TH1F("hNum","", nbinstrack, trackbins);
  auto hNum2  = new TH1F("hNum2","", nbinstrack, trackbins);
  auto hFake = new TH1F("hFake", "", nbinstrack, trackbins);
  auto hReco = new TH1F("hReco","", nbinstrack,trackbins);
  auto hTrack_pt = new TH1F("hTrack_pt","", nbinstrack, trackbins);
  auto hNum2Deta = new TH2F("hNum2Deta","", nbinseta, etabins, nbinstrack, trackbins);
  auto hDen2Deta = new TH2F("hDen2Deta","", nbinseta, etabins, nbinstrack, trackbins);
  auto hEff2Deta = new TH2F("hEff2Deta","", nbinseta, etabins, nbinstrack, trackbins);
  auto hNum2Dphi = new TH2F("hNum2Dphi","", nbinsphi, phibins, nbinstrack, trackbins);
  auto hDen2Dphi = new TH2F("hDen2Dphi","", nbinsphi, phibins, nbinstrack, trackbins);
  auto hEff2Dphi = new TH2F("hEff2Dphi","", nbinsphi, phibins, nbinstrack, trackbins);
 
  auto hNum2D = new TH2F("hNum2D","", nbinsphi, phibins, nbinseta, etabins);
  auto hDen2D = new TH2F("hDen2D","", nbinsphi, phibins, nbinseta, etabins);


  auto hIsoCorrelation = new TH2F("hIsoCorrelation", "", 80, 0, 10.0, 80, 0, 10.0); 
  auto hIsoReco        = new TH1F("hIsoReco", "", 200, 0, 100);
  auto hIsoTrue        = new TH1F("hIsoTrue", "", 200, 0, 100);
  auto hCorrelation_dEtapT = new TH2F("hCorrelation_dEtapT", "" , 100, -0.1, 0.1, 100, 0.1, 10.0);
  auto hCorrelation_dPhipT = new TH2F("hCorrelation_dPhipT", "" , 100, -0.1, 0.1, 100, 0.1, 10.0);

  auto hDCA_xy = new TH1F("hDCA_xy", "", 500, -1.0, 1.0);
  auto hDCA_z = new TH1F("hDCA_z", "", 500, -1.0, 1.0);
  auto hDCA_xy_fake = new TH1F("hDCA_xy_fake", "", 500, -1.0, 1.0);
  auto hDCA_z_fake = new TH1F("hDCA_z_fake", "", 500, -1.0, 1.0);   
  
  auto hChi2 = new TH1F("hChi2","", 50, 0, 10.1);
  auto hChi2_fake = new TH1F("hChi2_fake","", 50, 0, 10.1);


  auto hITSclus = new TH1F("hITSclus", "", 7, -0.5, 6.5);
  auto hITSclus_fake = new TH1F("hITSclus_fake", "", 7, -0.5, 6.5);

  auto hIso_ITS = new TH1F("hIso_ITS","", 25, -10, 40);
  auto hIso_TPC = new TH1F("hIso_TPC","", 25, -10, 40);  
  auto hIso_Truth = new TH1F("hIso_Truth","", 25, -10, 40);

  auto hTrue_eta  = new TH1F("hTrue_eta", "", nbinseta, etabins);
  auto hReco_eta  = new TH1F("hReco_eta", "", nbinseta, etabins);
  auto hEta  = new TH1F("hEta", "", nbinseta, etabins);
  auto hEta_plus  = new TH1F("hEta_plus", "", nbinstrack,trackbins);
  auto hEta_minus  = new TH1F("hEta_minus", "", nbinstrack,trackbins);
  auto hTrue_phi = new TH1F("hTrue_phi","", nbinsphi, phibins);
  auto hReco_phi = new TH1F("hReco_phi","", nbinsphi, phibins);
  auto hPhi = new TH1F("hPhi","", nbinsphi, phibins);
  auto hReco2D = new TH2F("hReco2D","", nbinsphi, phibins, nbinseta, etabins);
  
  auto hTrackCut = new TH1F("hTrackCut", "", 10, -0.5, 9.5);
  auto hNumTracks = new TH1F("hNumTracks", "", 500, -0.5, 499.5);
  auto hEventCounts = new TH1F("hEventCounts","", 10, -0.5, 9.5);
  auto hZvertex = new TH1F("hZvertez", "", 60, -30, 30);

  hFake->Sumw2();  
  hTrack_pt->Sumw2();
  hEta->Sumw2();
  hPhi->Sumw2();
  hZvertex->Sumw2();

  hTrack_pt->SetTitle(";p_{T} GeV/c;dn/dp_{T}");
  hEta->SetTitle(";#eta;counts");
  hPhi->SetTitle(";#phi;counts");
  hEta_plus->SetTitle(";#p_{T} [GeV/c];counts");
  hEta_minus->SetTitle(";p_{T} [GeV/c];counts");

  //Photon
  auto h_Reco  = new TH1F("h_Reco","", 20, 10, 30);
  auto hCluster_pt = new TH1F("hCluster_pt", "", 20, 10, 30);

  h_Reco->Sumw2();
  hCluster_pt->Sumw2();

  hCluster_pt->SetTitle("; p_{T} (GeV) ; dN/dp_{T}");
  //Jets
  auto h_jetpt_reco_its = new TH1F("h_jetpt_reco_its", "reco jet reco pt", 30, 0, 30);
  auto h_jetEta_reco_its = new TH1F("h_jetEta_reco_its","", nbinseta, etabins);
  auto h_jetPhi_reco_its = new TH1F("h_jetPhi_reco_its","", nbinsphi, phibins);
  auto h_jet2D_reco_its = new TH2F("h_jet2D_reco_its","", nbinsphi, phibins, nbinseta, etabins);

  h_jetpt_reco_its->Sumw2();
  h_jetEta_reco_its->Sumw2();
  h_jetPhi_reco_its->Sumw2();
  h_jet2D_reco_its->Sumw2();

  h_jetpt_reco_its->SetTitle(";p_{T}^{reco};1/N_{evt}dN/dp_{T}");
  h_jetEta_reco_its->SetTitle(";#eta^{reco};counts");
  h_jetPhi_reco_its->SetTitle(";#phi^{reco};counts");
  h_jet2D_reco_its->SetTitle(";#phi^{reco};#eta^{reco}");

  int nevent = 0; 
  int numEvents_tracks = 0;
  int numEvents_clusters = 0;
  int numEvents_clusters2 = 0;
  int numEvents = 0;
  //const int TrackBit = 16; //ITSONLY==16; ITS--TPC with full-jet cuts

  const double maxEta = 0.8;
  const double jetptmin = 5.0;
  const bool doCutDCA = false; 
  const bool doCutChi2 = false;

  ULong64_t one1 = 1;
  ULong64_t triggerMask_13bc = (one1 << 6) | (one1 << 7) | (one1 << 8) | (one1 << 9) | (one1 << 38) | (one1 << 39) | (one1 << 40) | (one1 << 41) | (one1 << 42) | (one1 << 43) | (one1 << 41) | (one1 << 42) | (one1 << 43) | (one1 << 45) | (one1 << 46);//minbias 
  ULong64_t triggerMask_13d = (one1 << 18) | (one1 << 19) | (one1 << 20) | (one1 << 21);
  ULong64_t triggerMask_13ef = (one1 << 17) | (one1 << 18) | (one1 << 19) | (one1 << 20);

  ULong64_t triggerMask_17p[2];
  triggerMask_17p[0] = (one1 << 0) | (one1 << 3) | (one1 << 4) | (one1 << 5) | (one1 << 6) | (one1 << 26) | (one1 << 27) | (one1 << 28) | (one1 << 29);//minbais
  triggerMask_17p[1] = (one1 << 11);//minbias

  ULong64_t triggerMask_17q_group1_MB = (one1 << 6) | (one1 << 7) | (one1 << 8) | (one1 << 9) | (one1 << 21) | (one1 << 22) | (one1 << 23) | (one1 << 24) | (one1 << 34);//minbais
  ULong64_t triggerMask_17q_group1_EG2 = (one1 << 14) | (one1 << 27);//EG2 
  ULong64_t triggerMask_17q_group2_MB = (one1 << 6) | (one1 << 7) | (one1 << 8) | (one1 << 9) | (one1 << 19) | (one1 << 20) | (one1 << 21) | (one1 << 22) | (one1 << 30);//minbais
  ULong64_t triggerMask_17q_group2_EG2 = (one1 << 12) | (one1 << 23);//EG2
  ULong64_t triggerMask_17q_group3_MB = (one1 << 0) | (one1 << 3) | (one1 << 4) | (one1 << 5) | (one1 << 6) | (one1 << 26) | (one1 << 27) | (one1 << 28) | (one1 << 29) | (one1 << 61);//minbais
  ULong64_t triggerMask_17q_group3_EG2 = (one1 << 10) | (one1 << 45) | (one1 << 54);//EG2
 
  ULong64_t triggerMask_13d_MB = (one1 << 2) | (one1 << 3) | (one1 << 25) | (one1 << 26) | (one1 << 26) | (one1 << 36) | (one1 << 37) | (one1 << 38) | (one1 << 39) | (one1 << 40) | (one1 << 41) | (one1 << 42) | (one1 << 43) | (one1 << 44);//minbias 
  ULong64_t triggerMask_13e_MB = (one1 << 2) | (one1 << 3) | (one1 << 24) | (one1 << 25) | (one1 << 26) | (one1 << 36) | (one1 << 37) | (one1 << 38) | (one1 << 39) | (one1 << 40) | (one1 << 41) | (one1 << 42) | (one1 << 43) | (one1 << 44);//minbias
  ULong64_t triggerMask_13f_MB = (one1 << 2) | (one1 << 3) | (one1 << 24) | (one1 << 25) | (one1 << 26) | (one1 << 36) | (one1 << 37) | (one1 << 38) | (one1 << 39) | (one1 << 40) | (one1 << 41) | (one1 << 42) | (one1 << 43) | (one1 << 44);//minbias

  ULong64_t trigMask[2] = {0};
  
  if(filename(0,3) == "13b" || filename(0,3) == "13c")
    trigMask[0] = triggerMask_13bc;
  if(filename(0,3) == "13d")
    trigMask[0] = triggerMask_13d_MB;  
  if(filename(0,3) == "13f" || filename(0,3) == "13e")
    trigMask[0] = triggerMask_13ef;
  if(filename(0,3) == "17p")
    {trigMask[0] = triggerMask_17p[0]; trigMask[1] = triggerMask_17p[1];}   
  //  if(filename(0,3) == "17q")
  //  trigMask[0] = triggerMask_17q;

  //int 17q_group1[] = {282441, 282440, 282439, 282437, };
  

  //cout << triggerMask_13bc << "\t" << trigMask << endl;
  Long64_t totEvents = tree->GetEntries();
  Long64_t restrictEvents = 100;
  Long64_t numEntries = TMath::Min(totEvents,restrictEvents);
  std::cout << numEntries << std::endl;
  for (Long64_t ievent=0;ievent< numEntries ;ievent++) {
    tree->GetEntry(ievent);
    nevent += 1;
    if(nevent%10==0) //if(nevent%1000==0)
      {
	std::cout << nevent << std::endl;
	cout << run_number << endl;
      }

    
   bool eventChange = true;
        
    //Event Selection:
    if(not( TMath::Abs(primary_vertex[2])<10.0)) continue; //vertex z position
    if(primary_vertex[2] == 0.000000) continue;
    if(is_pileup_from_spd_5_08) continue; //removes pileup
    if(((trigMask[0] & trigger_mask[0]) == 0) && ((trigMask[1] & trigger_mask[1]) == 0)) continue; //trigger selection
    hZvertex->Fill(primary_vertex[2]);
    numEvents++;
    int eventFill = 0;

    //Loop over all the tracks    
    for (int n=0; n< ntrack; n++){ 
      hTrackCut->Fill(0);
      if((track_quality[n]&TrackBit)==0) continue; hTrackCut->Fill(1);//track quality cut
      if(TMath::Abs(track_eta[n])> maxEta) continue; hTrackCut->Fill(2);//eta cut
      if(track_pt[n] < 0.15) continue; hTrackCut->Fill(3);//pt cut
      if(track_its_chi_square[n]>36.0) continue; hTrackCut->Fill(4);
      if(TrackBit == 16)
	{
	  if(track_its_ncluster[n] < 4) continue; 
	  hTrackCut->Fill(5);//its cluster cut
	}
      //double DCAcut = 7*(27+50/TMath::Power(track_pt[n],1.01));//paper
      if(TMath::Abs(track_dca_xy[n]) > 2.4) continue; hTrackCut->Fill(6);
      if(TMath::Abs(track_dca_z[n]) > 3.2) continue; hTrackCut->Fill(7);
      
      if(track_eta[n] > 0)
	hEta_plus->Fill(track_pt[n]);
      if(track_eta[n] < 0)
	hEta_minus->Fill(track_pt[n]);
      
      eventFill = 1;
      if (eventChange && !isMC) {numEvents_tracks++; eventChange = false;}

      hReco->Fill(track_pt[n]);
      hTrack_pt->Fill(track_pt[n]);
      
      if(track_pt[n] < 1.0) continue;
      hEta->Fill(track_eta[n]);
      hPhi->Fill(track_phi[n]);
      hReco2D->Fill(track_phi[n], track_eta[n]);
      
    }//end ntrack loop
    
    hEventCounts->Fill(eventFill);

    bool eventChange2 = true;
    //Loop over clusters
    for(ULong64_t n=0; n< ncluster; n++)
      {
	if (eventChange2 && !isMC) {numEvents_clusters2++; eventChange2 = false;}

	//Photon Selection
	if( not(cluster_pt[n]>8)) {continue;} //select pt of photons
	if( not(cluster_ncell[n]>2)) continue;   //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>0.05)) continue; //removes "spiky" clusters
	if( not(cluster_nlocal_maxima[n]<= 2)) continue; //require to have at most 2 local maxima.
	if( not(cluster_distance_to_bad_channel[n]>=2.0)) continue;
	
	//Isolation and shower shape selection:
	if( not(cluster_iso_its_04[n] < 1.0)) continue;
	if( not(cluster_lambda_square[n][0]<0.27)) continue; //single-photon selection (as opposed to merged photon).
	
	
	if (eventChange && !isMC) {numEvents_clusters++; eventChange = false;}
	
	h_Reco->Fill(cluster_pt[n]);	
	hCluster_pt->Fill(cluster_pt[n]);
	hIso_ITS->Fill(cluster_iso_its_04[n]);
	hIso_TPC->Fill(cluster_iso_tpc_04[n]);
	hIso_Truth->Fill(cluster_iso_04_truth[n]);
      }
    
    //Loop over jets
    for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) {
      if(not (jet_ak04its_pt_raw[ijet]>jetptmin)) continue;
      if(not (TMath::Abs(jet_ak04its_eta_raw[ijet])  <0.5 ) ) continue;
      h_jetpt_reco_its->Fill(jet_ak04its_pt_raw[ijet]);
      h_jetEta_reco_its->Fill(jet_ak04its_eta_raw[ijet]);
      h_jetPhi_reco_its->Fill(jet_ak04its_phi_raw[ijet]);
      h_jet2D_reco_its->Fill(jet_ak04its_phi_raw[ijet], jet_ak04its_eta_raw[ijet]);
    }
    
  }//loop over events
  std::cout << " END LOOP  " << std::endl;

  auto gausfit = new TF1("gaus","gaus", -25,25);
  gausfit->SetLineColor(kRed);
  auto g_mean = new TGraphErrors();
  auto g_sigma = new TGraphErrors();
  
  auto c1 = new TCanvas();
   
  //Study of the ITS-only track pT resolution
  int nbins = 9;
  for(int i=0; i<nbins; i++){
    double minpt = bins[i];
    double maxpt = bins[i+1];
    double binwidth = maxpt-minpt;
    int minbin =  hRes_Pt->GetXaxis()->FindBin(minpt);
    int maxbin =  hRes_Pt->GetXaxis()->FindBin(maxpt);
    auto h1 = hRes_Pt->ProjectionY("h", minbin, maxbin);
    
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
    
    //c1->SaveAs(Form("PDFOUTPUT/projecting%i_TrackBit%i_%s.pdf", i, TrackBit,filename.Data()));
  }
  
  
  //Normalizing the bins and getting yaxsis to be 1/Nevt*dN/dptdeta
  cout << numEvents_tracks << endl;
  cout << filename(0,3) << "\tTotal Events: " << numEntries << "\tEvent selection: " << numEvents << "\tPre-Cluster selection: " << numEvents_clusters2 << "\tPostCluster selection: " << numEvents_clusters << endl;
  const double tot_eta = 1.6;
  for(int i = 1; i < hTrack_pt->GetNbinsX()+1; i++)
    {
      double dpt = hTrack_pt->GetBinWidth(i);
      //cout << dpt << endl;
      double content = hTrack_pt->GetBinContent(i);
      double temp = content/((double)numEvents*dpt*tot_eta);
      //double temp = content/dpt;
      hTrack_pt->SetBinContent(i, temp);
      
      double error = hTrack_pt->GetBinError(i);
      double tempErr = error/((double)numEvents*dpt*tot_eta);
      //double tempErr = error/dpt;
      hTrack_pt->SetBinError(i, tempErr);
    }
  
  //scaling the jets
  for(int i = 1; i < h_jetpt_reco_its->GetNbinsX()+1; i++)
    {
      double dpt = h_jetpt_reco_its->GetBinWidth(i);

      double content = h_jetpt_reco_its->GetBinContent(i);
      double temp = content/((double)numEvents*dpt*tot_eta);
      h_jetpt_reco_its->SetBinContent(i, temp);
      
      double error = h_jetpt_reco_its->GetBinError(i);
      double tempErr = error/((double)numEvents*dpt*tot_eta);
      h_jetpt_reco_its->SetBinError(i, tempErr);
    }

  
  auto c = new TCanvas();   
  
  //filename += "_pPb_1GeV30GeV_trig_pileup__zv_jets";
  filename += "junk_pPb_trig_pileup__zv_jets_1M";
  //filename += "_pp_minBiasTrigg_Allevents";
  auto fout = new TFile(Form("OutputData/fout_%i_%s.root",TrackBit, filename.Data()), "RECREATE");
  //   
    
  hTrackCut->GetXaxis()->SetBinLabel(1,"All Tracks");
  hTrackCut->GetXaxis()->SetBinLabel(2,"Track quality cut");
  hTrackCut->GetXaxis()->SetBinLabel(3,"Track #eta cut");
  hTrackCut->GetXaxis()->SetBinLabel(4,"pt cut");
  hTrackCut->GetXaxis()->SetBinLabel(5,"ITS #chi^{2} cut");
  hTrackCut->GetXaxis()->SetBinLabel(6,"ITS nCluster cut");
  hTrackCut->GetXaxis()->SetBinLabel(7,"DCAr cut");
  hTrackCut->GetXaxis()->SetBinLabel(8,"DCAz cut");

  
  hEventCounts->GetXaxis()->SetBinLabel(1, "Passing Event Selection only");
  hEventCounts->GetXaxis()->SetBinLabel(2, "Passing Track Selection");

  //writing tracking stuff
  // hTrack_pt->Write("hReco");
  // hEta->Write("hReco_eta");
  // hPhi->Write("hReco_phi");
  // hReco2D->Write("hReco_phiEta");
  // hEta_plus->Write("hPt_plusEta");
  // hEta_minus->Write("hPt_minusEta");
  // hTrackCut->Write("hTrackCut");
  
  // //writing photon info
  // hCluster_pt->Write("cluster_pt");
  // hEventCounts->Write("hEventCounts");
  
  // hZvertex->Write("hZvertex");


  //Writing out jet info
  h_jetpt_reco_its->Write("hReco");
  h_jetEta_reco_its->Write("hReco_eta");
  h_jetPhi_reco_its->Write("hReco_phi");
  h_jet2D_reco_its->Write("hReco_phiEta");
  fout->Close();
  
  
  hIso_ITS->SetTitle("; #Sigma_{R<0.4} p_{T} [GeV]; counts ");
  hIso_ITS->Draw("hist");
  //myText(0.18, 0.96, kRed, "ITS-only");
  //myText(0.18, 0.92, kBlack, "ITS+TPC");
  //myText(0.68, 0.85, kBlack, "p_{T}^{cluster} > 6 GeV");
  hIso_TPC->SetLineColorAlpha(2,0.5);
  hIso_TPC->Draw("histsame");
  //hIso_Truth->SetLineColorAlpha(4,0.5);
  // hIso_Truth->Draw("histsame");
  gPad->SetLogy(0); 
  //c->SaveAs(Form("PDFOUTPUT/Iso%s.pdf",filename.Data()));
  c->Clear();
  
  hIso_Truth->Draw("hist");
  c->SaveAs(Form("PDFOUTPUT/IsoTrue%s.pdf",filename.Data()));
  c->Clear();
  c->Close();
  c1->Close();
  
}

    
void efficiencyCalculation(){  
  
  //Run(16, "16c2");
  // Run(3, "16c2");
  
  //Run(16, "16c3b_small");
  //Run(3, "16c3b/16c3b_small");
  //Run(16, "17f4b");
  //Run(3, "17f4b");
  //Run(16, "17f4b_woSDD");
  //Run(3, "17f4b_woSDD");
  
  //Run(3, "17g8a_woSDD");
  //Run(16, "17g8a_woSDD");
  
  //Input to Run is as follow: Run(const int TrackBit, TString address, bool isMC, bool hasAliDir, bool triggered)
  
  ////////////////////////////////pPb////////////////////////////////////

  ////////////////////////////////MC////////////////////////////////////

  //Run(3, "MC/17g6a3/17g6a3_pthat2_clusterv2_small", true, false, false);
  //Run(16, "MC/17g6a3/17g6a3_pthat2_clusterv2_small", true, false, false);
  
  //Run(3, "MC/17g6a3/17g6a3_pthat2_latest.root", true, true, false);
  //Run(16, "MC/17g6a3/17g6a3_pthat2_latest.root", true, true, false);
  
  //Run(3, "MC/17g6a3/17g6a3_pthat6_v2_noClusCut.root", true, true, false);
  //Run(16, "MC/17g6a3/17g6a3_pthat6_v2_noClusCut.root", true, true, false);
  
  //Run(16, "MC/16c3c/16c3c_pthat2_v1.root", true, true, false);
  //Run(16, "MC/16c3c/16c3c_pthat3_v1.root", true, true, false);
  //Run(16, "MC/16c3c/16c3c_pthat6_v1.root", true, true, false);
  
  //Run(3, "MC/18b10a/18b10a_calo_pthat2.root", true, true, false);
  //Run(16, "MC/18b10a/18b10a_calo_pthat2.root", true, true, false);

  //Run(3, "MC/18b10b/18b10b_calo_pthat2.root", true, true, false);
  //Run(16, "MC/18b10b/18b10b_calo_pthat2.root", true, true, false);

  //////////////////////////////////Data///////////////////////////////////////
  //Run(3, "pPb/13b/13b_pass4_v2_1run", false, true, false);
  //Run(16, "pPb/13b/13b_pass4_v2_1run.root", false, true, false);
  
  //Run(3, "pPb/13b/13b_pass4_v2_3run.root", false, true, false);
  //Run(16, "pPb/13b/13b_pass4_v2_3run.root", false, true, false);
  
  //Run(3, "pPb/13b/13b_pass4_v1_3run.root", false, true, false);
  //Run(3, "pPb/13b/13b_pass4_v1_5run.root", false, true, false);
  //Run(16, "pPb/13b/13b_pass4_v1_5run.root", false, true, false);
  //Run(16, "pPb/13b/13b_trigTest.root", false, true, false);
  //Run(16, "pPb/13d/13d_trigTest.root", false, true, true);
  //Run(3, "pPb/13b/13b_3run.root", false, true, false);
  //Run(16, "pPb/13b/13b_3run.root", false, true, false);
  //Run(3, "pPb/13b/13b.root", false, true, false);
  //Run(16, "pPb/13b/13b.root", false, true, false);

  //Run(3, "pPb/13f/13f_pass4_v2_minbias_1run", false, true, true);
  //Run(16, "pPb/13f/13f_pass4_v2_minbias_1run", false, true, true);
  
  //Run(3, "pPb/13f/13f_pass4_minbias_v2_5run", false, true, true);
  //Run(16, "pPb/13f/13f_pass4_minbias_v2_5run", false, true, true);

  //Run(16, "pPb/13d/13d.root", false, true, true);
  //Run(16, "pPb/13e/13e.root", false, true, true);
  //Run(16, "pPb/13f/13f.root", false, true, true);
  //Run(16, "pPb/13def.root", false, true, true);
    
  //Run(3, "pPb/13c/13c_pass4_v2", false, true, false);
  //Run(16, "pPb/13c/13c_pass4_v2", false, true, false);
  
  //Run(3, "pPb/13c/13c_pass4_v1_0GevtrackSkim.root", false, true, false);
  //Run(16, "pPb/13c/13c_pass4_v1_0GevtrackSkim.root", false, true, false);

  //Run(3, "pp/17p/17p_CENTwSDD_v2_3runs_noClusCut.root", false, true, true);
  //Run(16, "pp/17p/17p_CENTwSDD_v2_3runs_noClusCut.root", false, true, true);

  //Run(3, "pp/17p/17p_CENT_wSDD_MB_v1_mini.root", false, true, true);
  //Run(16, "pp/17p/17p_CENT_wSDD_MB_v1_mini.root", false, true, true);
  //Run(16, "pp/17p/17p_CENT_wSDD_MB_v1.root", false, true, true);
  //Run(3, "pp/17p/17p_CENT_wSDD_v1_2run.root", false, true, true);
  //Run(16, "pp/17p/17p_CENT_wSDD_v1_2run.root", false, true, true);
  //Run(16, "pp/17p/17p.root", false, true, true);


  //Run(3, "pp/17q/17q_v1_CENTwSDD_noClusCut_3run.root", false, true, true);
  //Run(16, "pp/17q/17q_v1_CENTwSDD_noClusCut_3run.root", false, true, true);
  Run(16, "pp/17q/17q_wSDD.root", false, true, true);



  return;
}




