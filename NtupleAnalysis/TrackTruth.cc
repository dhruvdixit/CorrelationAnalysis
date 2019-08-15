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
#include <TStyle.h>


#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

using namespace std;

double SetPthatWeights(TString MCname, double Xsection, double ntrial)
{

  //13b2 and 17l4_cent weights
  if(MCname(0,4) == "13b2" || MCname(0,4) == "17l4")
    return 1.0;
   
  //17g6a3 weights
  if(MCname == "17g6a3_pthat1")
    return 4.47e-11;
  if(MCname == "17g6a3_pthat2")
    return 9.83e-11;
  if(MCname == "17g6a3_pthat3")
    return 1.04e-10;
  if(MCname == "17g6a3_pthat4")
    return 1.01e-10;
  if(MCname == "17g6a3_pthat5")
    return 6.93e-11;
  if(MCname == "17g6a3_pthat6")
    return 5.13e-11;
  if(MCname == "17g6a3_pthat7")
    return 3.03e-11;
  if(MCname == "17g6a3_pthat8")
    return 1.89e-11;
  
  //17g6a1 weights
  if(MCname == "17g6a1_pthat1")
    return 1.60e-11;
  if(MCname == "17g6a1_pthat2")
    return 2.72e-12;
  if(MCname == "17g6a1_pthat3")
    return 3.69e-13;
  if(MCname == "17g6a1_pthat4")
    return 6.14e-14;
  if(MCname == "17g6a1_pthat5")
    return 1.27e-14;
  

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
  if(MCname(0,2) == "18")
    return Xsection/ntrial;
 
  return 0.0;
  
}

// bool IsTracking(TString MCname)
// {
  
//   if(MCname == "13b2") return true;
//   if(MCname(0,2) == "18") return false;
//   if(MCname == "17l4") return true;
//   if(MCname == "16c3c") return true;
//   if(MCname == "17g6a3") return true;
//   if(MCname == "17g6a1") return false;
  
//   return false;

// }

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];
    
  dummyv[0] = strdup("main");

  gStyle->SetOptStat("");
  
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
  
  const int nbinstrack = 31;
  Double_t trackbins[nbinstrack+1] = {};
  //Double_t trackbins[nbinstrack+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.00, 24.00, 26.00, 30.00};//pPB nbinstrack = 21
  /*Double_t trackbins[nbinstrack+1] = {
    0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,  0.50,  0.55,  0.60, 
    0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,  1.00,  1.10,  1.20,
    1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90,  2.00,  2.20,  2.40,
    2.60,  2.80,  3.00,  3.20,  3.40,  3.60,  3.80,  4.00,  4.50,  5.00,
    5.50,  6.00,  6.50,  7.00,  8.00,  9.00,  10.00, 11.00, 12.00, 13.00,
    14.00, 15.00, 16.00, 18.00, 20.00, 22.00, 24.00, 26.00, 30.00};//nbinsbstrack = 58*/
  /*Double_t trackbins[nbinstrack+1] = {
    0.15,  0.25,  0.50,  0.75,  1.00,  2.00,  3.00,  4.00,  5.00,  6.00,  
    7.00,  8.00,  9.00,  10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00,
    18.00, 20.00, 22.00, 24.00, 26.00, 30.00};//nbinsbstrack = 25*/
  //Double_t trackbins[nbinstrack+1] = {1.0,2.0,3.0,4.0,5.0,6.0, 8.0, 10.0, 13.0, 20.0};//pp bining
  double ptmin = 0.5;
  double ptmax = 16;
  double ptstep = (ptmax-ptmin)/nbinstrack;
  for(int i = 0; i < nbinstrack+1; i++){
    trackbins[i] = ptmin + i*ptstep;
  }
  
  // const int nbinscluster = 55;
  // Double_t clusterbins[nbinscluster+1] = {};
  // //Double_t clusterbins[nbinscluster+1] = {12.0 , 14.06198608, 16.47828771, 19.30978769, 22.62783047, 26.51601976, 31.07232506, 36.4115502 , 42.668226  , 50.0};//geom binning
  // //Double_t clusterbins[nbinscluster+1]{12.0, 15.96389997, 21.23717518, 28.252345, 37.58480079, 50.0};
  // double Emin = 5;
  // double Emax = 60;
  // double Estep = (Emax-Emin)/nbinscluster;
  // for(int i = 0; i < nbinscluster+1; i++)
  //   {
  //     clusterbins[i] = Emin + i*Estep;
  //   }

  // TH1D h_Den("h_Den", "", nbinscluster, clusterbins);
  // TH1D h_Den_emcal("h_Den_emcal", "", nbinscluster, clusterbins);
  // TH1D h_Den_dcal("h_Den_dcal", "", nbinscluster, clusterbins);
  // TH1D h_Num("h_Num", "", nbinscluster, clusterbins);
  // TH1D h_Num_emcal("h_Num_emcal", "", nbinscluster, clusterbins);
  // TH1D h_Num_dcal("h_Num_dcal", "", nbinscluster, clusterbins);
  // TH1D h_Reco("h_Reco", "", nbinscluster, clusterbins);
  // TH1D h_Reco_emcal("h_Reco_emcal", "", nbinscluster, clusterbins);
  // TH1D h_Reco_dcal("h_Reco_dcal", "", nbinscluster, clusterbins);
  // TH1D hClusterCut("hClusterCut", "", 20, -0.5, 19.5);
  // TH1D hCluster_iso_04_truth("hCluster_iso_04_truth", "", 20 , 0, 18);
  // TH2D h_Correlation("h_Correlation", "", nbinscluster, clusterbins, nbinscluster, clusterbins);
  // //TH2D h_Num2D("h_Num2D","", nbinsphi, phibins, nbinseta, etabins);
  // //TH2D h_Den2D("h_Den2D","", nbinsphi, phibins, nbinseta, etabins);
 
  // TH2D h_Num2D("h_Num2D","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  // TH2D h_Den2D("h_Den2D","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  // TH2D h_Num2D_b("h_Num2D_b","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  // TH2D h_Den2D_b("h_Den2D_b","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  
  // h_Den.Sumw2();
  // h_Num.Sumw2();
  // h_Den_emcal.Sumw2();
  // h_Num_emcal.Sumw2();
  // h_Den_dcal.Sumw2();
  // h_Num_dcal.Sumw2();
  // h_Reco.Sumw2();
  // h_Reco_emcal.Sumw2();
  // h_Reco_dcal.Sumw2();
  // h_Correlation.Sumw2();
  
  // h_Den.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  // h_Num.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  // h_Den_emcal.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  // h_Num_emcal.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  // h_Den_dcal.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  // h_Num_dcal.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  // h_Reco.SetTitle("reco photons filled with pt reco; p_{T}^{Reco,pt} [GeV]; 1/N_{event} dN/dp_{T}");
  // h_Reco_emcal.SetTitle("reco clustes filled with pt reco; p_{T}^{Reco,pt} [GeV];1/N_{event} dN/dp_{T}");
  // h_Reco_dcal.SetTitle("reco clustes filled with pt reco; p_{T}^{Reco,pt} [GeV];1/N_{event} dN/dp_{T}");
  // h_Correlation.SetTitle("; True p_{T} [GeV]; Reconstructed p_{T} [GeV]");
  // h_Num2D.SetTitle(";#phi_{true}; #eta_{true}");
  // h_Den2D.SetTitle(";#phi_{true}; #eta_{true}");
  // h_Num2D_b.SetTitle(";#phi_{true}; #eta_{true}");
  // h_Den2D_b.SetTitle(";#phi_{true}; #eta_{true}");


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
  TH2F hRes_Pt("hRes_Pt", "", 300, ptmin, ptmax, 80, -50, 50);
 

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
  
  // //jets
  // TH1D hBR_jetpt("hBR_jetpt",   "Associated jet pt spectrum (reco), bkg region", 30, 0, 30);
  // TH1D hSR_jetpt("hSR_jetpt", "Associated jet pt spectrum (reco), signal region", 30, 0, 30);
  // TH1D hBR_jetpt_truth("hBR_jetpt_truth",   "Associated jet pt spectrum (truth), bkg region", 30, 0, 30);
  // TH1D hSR_jetpt_truth("hSR_jetpt_truth", "Associated jet pt spectrum (truth), signal region", 30, 0, 30);
  // TH1D h_jetpt_truth("h_jetpt_truth", "truth jet pt", 30, 0, 30);
  // TH1D h_jetpt_truthreco("h_jetpt_truthreco", "reco jet truth pt (numerator of efficiency)", 30, 0, 30);
  // TH1D h_jetpt_reco("h_jetpt_reco", "reco jet reco pt", 30, 0, 30);
  // TH2D h_jetpt_correlation("h_jetpt_correlation", "jet response matrix", 30, 0, 30, 30, 0, 30);
  // TH2F h_jetRes_Pt("h_jetRes_Pt", "", 30, 0, 30, 320, -250, 150);

  // h_jetpt_truth.Sumw2();
  // h_jetpt_truthreco.Sumw2();
  // h_jetpt_reco.Sumw2();
  // hBR_jetpt.Sumw2();
  // hSR_jetpt.Sumw2();
  // hBR_jetpt_truth.Sumw2();
  // hSR_jetpt_truth.Sumw2();

  // h_jetRes_Pt.SetTitle("Jet Resolution Response;p_{T}^{true} [GeV/c];(p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%]");

  const int TrackBit = 3;//3 for TPC+ITS, 16 for ITS only
  TString ntupleName = "junk";
  TString MCname = "junk";
  int numEvents, numEvents_tracks, numEvents_clusters;
  numEvents = numEvents_tracks = numEvents_clusters = 0;
  double sumWeight = 0.0;

  //const int numMC = argc;
  //double aveXsectionArray[] = {0.0};
  
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
    TString temp = (TString)argv[iarg];
    ntupleName = temp(temp.Last('/')+1,temp.First('_')-temp.Last('/')+6);
    MCname = ntupleName(0,ntupleName.Length()-7);


    // Get all the TTree variables from the file to open, I guess
    TTree *_tree_event = NULL;
    if(file->Get("AliAnalysisTaskNTGJ"))
      _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>  (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));

    if (_tree_event == NULL) {
      std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
      _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " fail " << std::endl;
	exit(EXIT_FAILURE);
      }
    } 

        
    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }
            
    //you define variables
    Double_t primary_vertex[3];
    Float_t ue_estimate_its_const;

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
        
    // UInt_t ncluster;
    // Float_t cluster_e[NTRACK_MAX];
    // Float_t cluster_e_cross[NTRACK_MAX];
    // Float_t cluster_pt[NTRACK_MAX];
    // Float_t cluster_eta[NTRACK_MAX];
    // Float_t cluster_phi[NTRACK_MAX];
    // //Float_t cluster_iso_tpc_04[NTRACK_MAX];
    // Float_t cluster_iso_its_04[NTRACK_MAX];
    // Float_t cluster_iso_its_04_ue[NTRACK_MAX];
    // Float_t cluster_iso_04_truth[NTRACK_MAX];
    // //Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    // Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    // Float_t cluster_s_nphoton[NTRACK_MAX][4];
    // UChar_t cluster_nlocal_maxima[NTRACK_MAX];
    // Float_t cluster_distance_to_bad_channel[NTRACK_MAX];   
 
    // unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    // Int_t cluster_ncell[NTRACK_MAX];
    // UShort_t  cluster_cell_id_max[NTRACK_MAX];
    // Float_t cluster_lambda_square[NTRACK_MAX][2];
    // Float_t cell_e[17664];
        
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

    // //Jets reco 
    // UInt_t njet_ak04its;
    // Float_t jet_ak04its_pt_raw[NTRACK_MAX];
    // Float_t jet_ak04its_eta_raw[NTRACK_MAX];
    // Float_t jet_ak04its_phi[NTRACK_MAX];
    
    // Float_t jet_ak04its_pt_truth[NTRACK_MAX];
    // Float_t jet_ak04its_eta_truth[NTRACK_MAX];
    // Float_t jet_ak04its_phi_truth[NTRACK_MAX];
  
    // //The z_reco is defined as the fraction of the true jet that ended up in this reco jet 
    // //There are two entries and indices, the first is the best. 
    // Int_t   jet_ak04its_truth_index_z_reco[NTRACK_MAX][2];
    // Float_t jet_ak04its_truth_z_reco[NTRACK_MAX][2];
    // Float_t jet_ak04its_ptd_raw[NTRACK_MAX];
    // Float_t jet_ak04its_width_sigma[NTRACK_MAX][2];
    // UShort_t jet_ak04its_multiplicity[NTRACK_MAX];

    // //Truth Jets
    // UInt_t njet_truth_ak04;
    // Float_t jet_truth_ak04_pt[NTRACK_MAX];
    // Float_t jet_truth_ak04_eta[NTRACK_MAX];
    // Float_t jet_truth_ak04_phi[NTRACK_MAX];    
    
        
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);

    //tracks
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
        
    //clusters
    // _tree_event->SetBranchAddress("ncluster", &ncluster);
    // _tree_event->SetBranchAddress("cluster_e", cluster_e);
    // _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    // _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
    // _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    // _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    // _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
    // _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    // _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    // _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    // _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);
    // _tree_event->SetBranchAddress("cluster_iso_04_truth",cluster_iso_04_truth);
    // _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    // _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);        
    // _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);

    // _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    // _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    // _tree_event->SetBranchAddress("cell_e", cell_e);
        
    // //MC
    _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
    _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);        
    _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",mc_truth_first_parent_pdg_code);
    _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
   
    // //jets
    // _tree_event->SetBranchAddress("njet_ak04its", &njet_ak04its);
    // _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
    // _tree_event->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
    // _tree_event->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi);
    // _tree_event->SetBranchAddress("jet_ak04its_pt_truth", jet_ak04its_pt_truth);
    // _tree_event->SetBranchAddress("jet_ak04its_eta_truth", jet_ak04its_eta_truth);
    // _tree_event->SetBranchAddress("jet_ak04its_phi_truth", jet_ak04its_phi_truth);

    // //quark-gluon discriminator variables
    // _tree_event->SetBranchAddress("jet_ak04its_ptd_raw", jet_ak04its_ptd_raw);
    // _tree_event->SetBranchAddress("jet_ak04its_width_sigma", jet_ak04its_width_sigma);
    // _tree_event->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_ak04its_multiplicity);



    // _tree_event->SetBranchAddress("jet_ak04its_truth_index_z_reco",     jet_ak04its_truth_index_z_reco);
    // _tree_event->SetBranchAddress("jet_ak04its_truth_z_reco", jet_ak04its_truth_z_reco);    

    // //truth jets
    // _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04);
    // _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_truth_ak04_pt);
    // _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_truth_ak04_phi);
    // _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_truth_ak04_eta); 

    const double maxEta = 0.8;
    Long64_t totEvents = _tree_event->GetEntries();
    Long64_t restrictEvents = 2000000;
    Long64_t numEntries = TMath::Min(totEvents,restrictEvents);
    cout << numEntries << endl;
    double aveXsection = 0.0;
    // Loop over events
    for(Long64_t ievent = 0; ievent < numEntries ; ievent++){
      _tree_event->GetEntry(ievent);
      
      bool eventChange = true;
      aveXsection += (double)eg_cross_section;

      //Selecting pthat weights
      double weight = SetPthatWeights(ntupleName, (double)eg_cross_section, (double)eg_ntrial);
      if(ievent%10000 == 0){
	cout << ievent << "\t" << totEvents << endl;
	cout << weight << endl;
	}
      sumWeight += weight;
      

      //event selection
      if(not(TMath::Abs(primary_vertex[2])<10.0)) continue; //vertex z position
      //if(primary_vertex[2] == 0.000000000) continue;
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
	  if(track_pt[n] > 1){
	    hRecoTruth_eta.Fill(mc_truth_eta[index]);
	    hRecoTruth_phi.Fill(mc_truth_phi[index]);
	    hRecoTruth2D.Fill(mc_truth_phi[index], mc_truth_eta[index]);
	  }
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
	if(mc_truth_pt[nTru] > 1)
	  {
	    hTruth_eta.Fill(mc_truth_eta[nTru]);
	    hTruth_phi.Fill(mc_truth_phi[nTru]);
	    hTruth2D.Fill(mc_truth_phi[nTru], mc_truth_eta[nTru]);
	  }
      }//end loop over MC particles       
    
    }//end loop over events
    
    //cout << (TString)argv[iarg] << "\t" << "average Xsection:\t" << aveXsection << "\tNumber of events:\t" << numEntries << "\tAveXsection:\t" << aveXsection/numEntries << endl;
    
  }//end loop over ntuples
  cout << numEvents_tracks << endl;
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
  hTrackCut.GetXaxis()->SetBinLabel(5,"ITS #chi^{2} cut");
  hTrackCut.GetXaxis()->SetBinLabel(6,"ITS nCluster cut");
  hTrackCut.GetXaxis()->SetBinLabel(7,"DCAr cut");
  hTrackCut.GetXaxis()->SetBinLabel(8,"DCAz cut");

  TF1* gausfit = new TF1("gaus","gaus", -25,25);
  gausfit->SetLineColor(kRed);
  TGraphErrors* g_mean = new TGraphErrors();
  TGraphErrors* g_sigma = new TGraphErrors();
  
  auto c1 = new TCanvas();
   
  //Study of the ITS-only track pT resolution
  //const Double_t bins[10] = {  0.1,          0.16681005,   0.27825594,   0.46415888 ,  0.77426368, 1.29154967,   2.15443469,   3.59381366,   5.9948425,   10.0};
  
  int nbins = nbinstrack;
  for(int i=0; i<nbins; i++){
    double minpt = trackbins[i];
    double maxpt = trackbins[i+1];
    double binwidth = maxpt-minpt;
    int minbin =  hRes_Pt.GetXaxis()->FindBin(minpt);
    int maxbin =  hRes_Pt.GetXaxis()->FindBin(maxpt);
    TH1D* h1 = hRes_Pt.ProjectionY("h", minbin, maxbin);
    
    h1->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    h1->Draw();  
    h1->GetYaxis()->SetNdivisions(5);
    h1->GetXaxis()->SetNdivisions(5);
    h1->GetYaxis()->SetTitle("counts");
    h1->Fit(gausfit,"RN0");
    h1->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    //gausfit->Draw("same");
    //myText(0.18, 0.8, kBlack, Form("%2.1f < p_{T}^{truth} < %2.1f GeV", minpt, maxpt));
    //myText(0.18, 0.74, kRed, Form("#mu = %2.1f [%]", gausfit->GetParameter(1)));
    //myText(0.18, 0.68, kRed, Form("#sigma = %2.1f [%]", gausfit->GetParameter(2))); 
    g_sigma->SetPoint(g_sigma->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(2));     
    g_sigma->SetPointError(g_sigma->GetN()-1, binwidth/2.0, gausfit->GetParError(2));
    g_mean->SetPoint(g_mean->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(1));
    g_mean->SetPointError(g_mean->GetN()-1, binwidth/2.0, gausfit->GetParError(1));
    
    //if(forTracking)
      //c1->SaveAs(Form("TrackOutput/PDFOUTPUT/%s_projecting%i_TrackBit%i.C", MCname.Data(), i, TrackBit));
    }//*/
  
  g_sigma->SetTitle("Relative resolution vs p_{T} ; p_{T}^{true} [GeV]; #sigma(p_{T})/p_{T} [%]"); 
  g_mean->SetTitle("; p_{T}^{true} [GeV]; Relative bias [%]");

  TGraphErrors* g_mean_jet = new TGraphErrors();
  TGraphErrors* g_sigma_jet = new TGraphErrors();

  bool makeXsectionFile = false;

  TFile* fout_track = new TFile(Form("TrackOutput/%s_%i_31bins_halfGeV16GeV_4L.root", MCname.Data(), TrackBit),"RECREATE");
  
  TGraphAsymmErrors* eff = new TGraphAsymmErrors(&hRecoTruth, &hTruth);
  eff->SetTitle("; p_{T}^{true} ; #epsilon");
  eff->Write("Efficiency");
  g_sigma->Write("g_sigma");
  g_mean->Write("g_mean");
  
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


  if(makeXsectionFile)
    {
      
    }
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}//end main

