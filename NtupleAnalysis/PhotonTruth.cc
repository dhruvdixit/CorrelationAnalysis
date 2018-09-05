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
    return 1.0;//Xsection/ntrial;
   
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
  //if(MCname == "16c3c_pthat1")
  //  return 3.941701e-03;
  //if(MCname == "16c3c_pthat2")
  //  return 2.001984e-03;
  //if(MCname == "16c3c_pthat3")
  //  return 9.862765e-04 ;
  //if(MCname == "16c3c_pthat4")
  //  return 9.862765e-04 ;

  //18b10ab, 18g7a
  if(MCname(0,2) == "18" || MCname(0,2) == "16")
    return Xsection/ntrial;
 
  return 0.0;
  
}

bool IsTracking(TString MCname)
{
  
  if(MCname == "13b2") return true;
  if(MCname(0,2) == "18") return false;
  if(MCname == "17l4") return true;
  if(MCname == "16c3c") return true;
  if(MCname == "17g6a3") return true;
  if(MCname == "17g6a1") return false;
  
  return false;

}

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
  
  const int nbinstrack = 21;
  //Double_t trackbins[nbinstrack+1] = {};
  Double_t trackbins[nbinstrack+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.00, 24.00, 26.00, 30.00};//pPB nbinstrack = 21
  //Double_t trackbins[nbinstrack+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 12.0, 20.0, 30.00};//pPB nbinstrack = 10
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
  //Double_t trackbins[nbinstrack+1] = {1.0,2.0,3.0,4.0,5.0,6.0, 8.0, 10.0, 13.0, 20.0};//pp bining nbinstrack = 9
  double ptmin = 1;
  double ptmax = 16;
  double ptstep = (ptmax-ptmin)/nbinstrack;
  //for(int i = 0; i < nbinstrack+1; i++){
  //  trackbins[i] = ptmin + i*ptstep;
  //}
  
  const int nbinscluster = 55;
  Double_t clusterbins[nbinscluster+1] = {};
  //Double_t clusterbins[nbinscluster+1] = {12.0 , 14.06198608, 16.47828771, 19.30978769, 22.62783047, 26.51601976, 31.07232506, 36.4115502 , 42.668226  , 50.0};//geom binning
  //Double_t clusterbins[nbinscluster+1]{12.0, 15.96389997, 21.23717518, 28.252345, 37.58480079, 50.0};
  double Emin = 5;
  double Emax = 60;
  double Estep = (Emax-Emin)/nbinscluster;
  for(int i = 0; i < nbinscluster+1; i++)
    {
      clusterbins[i] = Emin + i*Estep;
    }

  TH1D h_Den("h_Den", "", nbinscluster, clusterbins);
  TH1D h_Den_emcal("h_Den_emcal", "", nbinscluster, clusterbins);
  TH1D h_Den_dcal("h_Den_dcal", "", nbinscluster, clusterbins);
  TH1D h_Num("h_Num", "", nbinscluster, clusterbins);
  TH1D h_Num_emcal("h_Num_emcal", "", nbinscluster, clusterbins);
  TH1D h_Num_dcal("h_Num_dcal", "", nbinscluster, clusterbins);
  TH1D h_Reco("h_Reco", "", nbinscluster, clusterbins);
  TH1D h_Reco_emcal("h_Reco_emcal", "", nbinscluster, clusterbins);
  TH1D h_Reco_dcal("h_Reco_dcal", "", nbinscluster, clusterbins);
  TH1D h_YesISO("h_YesISO", "", nbinscluster, clusterbins);
  TH1D h_NoISO("h_NoISO", "", nbinscluster, clusterbins);
  TH1D hClusterCut("hClusterCut", "", 20, -0.5, 19.5);
  TH1D hCluster_iso_04_truth("hCluster_iso_04_truth", "", 20 , 0, 18);
  TH2D h_Correlation("h_Correlation", "", nbinscluster, clusterbins, nbinscluster, clusterbins);
  //TH2D h_Num2D("h_Num2D","", nbinsphi, phibins, nbinseta, etabins);
  //TH2D h_Den2D("h_Den2D","", nbinsphi, phibins, nbinseta, etabins);
 
  TH2D h_Num2D("h_Num2D","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  TH2D h_Den2D("h_Den2D","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  TH2D h_Num2D_b("h_Num2D_b","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  TH2D h_Den2D_b("h_Den2D_b","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  
  h_Den.Sumw2();
  h_Num.Sumw2();
  h_Den_emcal.Sumw2();
  h_Num_emcal.Sumw2();
  h_Den_dcal.Sumw2();
  h_Num_dcal.Sumw2();
  h_Reco.Sumw2();
  h_Reco_emcal.Sumw2();
  h_Reco_dcal.Sumw2();
  h_YesISO.Sumw2();
  h_NoISO.Sumw2();
  h_Correlation.Sumw2();
  
  h_Den.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  h_Num.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  h_Den_emcal.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  h_Num_emcal.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  h_Den_dcal.SetTitle("truth photons; p_{T}^{truth} [GeV]; entries");
  h_Num_dcal.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV]; entries");
  h_Reco.SetTitle("reco photons filled with pt reco; p_{T}^{Reco,pt} [GeV]; 1/N_{event} dN/dp_{T}");
  h_Reco_emcal.SetTitle("reco clustes filled with pt reco; p_{T}^{Reco,pt} [GeV];1/N_{event} dN/dp_{T}");
  h_Reco_dcal.SetTitle("reco clustes filled with pt reco; p_{T}^{Reco,pt} [GeV];1/N_{event} dN/dp_{T}");
  h_YesISO.SetTitle(";p_{T}^{Reco,pt} [GeV]; 1/N_{event} dN/dp_{T}");
  h_NoISO.SetTitle(";p_{T}^{Reco,pt} [GeV]; 1/N_{event} dN/dp_{T}");
  h_Correlation.SetTitle("; True p_{T} [GeV]; Reconstructed p_{T} [GeV]");
  h_Num2D.SetTitle(";#phi_{true}; #eta_{true}");
  h_Den2D.SetTitle(";#phi_{true}; #eta_{true}");
  h_Num2D_b.SetTitle(";#phi_{true}; #eta_{true}");
  h_Den2D_b.SetTitle(";#phi_{true}; #eta_{true}");


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
  TH2D hRecoTruth2DPtEta("hRecoTruth2DPtEta","", nbinstrack, trackbins, nbinseta, etabins);
  TH2D hReco2DPtEta("hReco2DPtEta","", nbinstrack, trackbins, nbinseta, etabins);
  TH2D hTruth2DPtEta("hTruth2DPtEta","", nbinstrack, trackbins, nbinseta, etabins);
  TH2F hCorrelation("hCorrelation", "", nbinstrack, trackbins, nbinstrack, trackbins);
  TH2F hCorrelation_cor("hCorrelation_cor", "", nbinstrack, trackbins, nbinstrack, trackbins);
  
  TH2F hRes_Pt("hRes_Pt", "", nbinstrack, trackbins, 80, -50, 50);
 

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
  hTruth2D.Sumw2();
  hRecoTruth2D.Sumw2();
  hReco2D.Sumw2();
  hTruth2DPtEta.Sumw2();
  hRecoTruth2DPtEta.Sumw2();
  hReco2DPtEta.Sumw2();
  hCorrelation.Sumw2();
  hCorrelation_cor.Sumw2();
  hZvertex.Sumw2();


  hTruth.SetTitle("; p_{T}^{true} [GeV/c]; entries");
  hRecoTruth.SetTitle("; p_{T}^{Reco,Truth} [GeV/c]; entries");
  hReco.SetTitle("; p_{T}^{Reco} [GeV/c]; entries");
  hFake.SetTitle("; p_{T}^{reco} [GeV/c]; Fake Rate");
  hCorrelation.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
  hCorrelation_cor.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
  hRecoTruth2D.SetTitle(";#phi;#eta");
  hReco2D.SetTitle(";#phi;#eta");
  hTruth2D.SetTitle(";#phi;#eta");
  hRecoTruth2DPtEta.SetTitle(";p_{T};#eta");
  hReco2DPtEta.SetTitle(";p_{T};#eta");
  hTruth2DPtEta.SetTitle(";p_{T};#eta");
  hZvertex.SetTitle(";Z_{v} [cm]; counts");
  


  //jets
  //tpc jets
  TH1D h_jetpt_truth_tpc("h_jetpt_truth_tpc", "truth jet pt", 30, 0, 30);
  TH1D h_jetpt_truthreco_tpc("h_jetpt_truthreco_tpc", "reco jet truth pt (numerator of efficiency)", 30, 0, 30);
  TH1D h_jetpt_reco_tpc("h_jetpt_reco_tpc", "reco jet reco pt", 30, 0, 30);
  
  TH2D h_jetpt_correlation_tpc("h_jetpt_correlation_tpc", "jet response matrix", 30, 0, 30, 60, -30, 30);
  TH2D h_jetpt_correlation2_tpc("h_jetpt_correlation2_tpc", "jet response matrix", 30, 0, 30, 60,-30, 30);
  TH2F h_jetRes_Pt_tpc("h_jetRes_Pt_tpc", "", 30, 0, 30, 150, -450, 150);//bin_width = 4

  //its jets
  TH1D h_jetpt_truth("h_jetpt_truth", "truth jet pt", 30, 0, 30);
  TH1D h_jetpt_truth_its("h_jetpt_truth_its", "truth jet pt its", 30, 0, 30);
  TH1D h_jetpt_truthreco_its("h_jetpt_truthreco_its", "reco jet truth pt (numerator of efficiency)", 30, 0, 30);
  TH1D h_jetpt_reco_its("h_jetpt_reco_its", "reco jet reco pt", 30, 0, 30);
  TH1F h_jetEta_truth_its("h_jetEta_truth_its","", nbinseta, etabins);
  TH1F h_jetEta_truthreco_its("h_jetEta_truthreco_its","", nbinseta, etabins);
  TH1F h_jetEta_reco_its("h_jetEta_reco_its","", nbinseta, etabins);
  TH1F h_jetPhi_truth_its("h_jetPhi_truth_its","", nbinsphi, phibins);
  TH1F h_jetPhi_truthreco_its("h_jetPhi_truthreco_its","", nbinsphi, phibins);
  TH1F h_jetPhi_reco_its("h_jetPhi_reco_its","", nbinsphi, phibins);
  
  TH2D h_jet2D_truth_its("h_jet2D_truth_its","", nbinsphi, phibins, nbinseta, etabins);
  TH2D h_jet2D_truthreco_its("h_jet2D_truthreco_its","", nbinsphi, phibins, nbinseta, etabins);
  TH2D h_jet2D_reco_its("h_jet2D_reco_its","", nbinsphi, phibins, nbinseta, etabins);
  TH2D h_jetpt_correlation_its("h_jetpt_correlation_its", "jet response matrix", 30, 0, 30, 60, -30, 30);
  TH2D h_jetpt_correlation2_its("h_jetpt_correlation2_its", "jet response matrix", 30, 0, 30, 60,-30, 30);
  TH2F h_jetRes_Pt_its("h_jetRes_Pt_its", "", 30, 0, 30, 150, -450, 150);//bin_width = 4
  TH2F h_jetRes_Phi_its("h_jetRes_Phi_its", "", nbinsphi, phibins, 800, -TMath::Pi()/2, TMath::Pi()/2);//bin_width = 4


  h_jetpt_truth_its.Sumw2();
  h_jetpt_truthreco_its.Sumw2();
  h_jetpt_reco_its.Sumw2();
  h_jetpt_truth_tpc.Sumw2();
  h_jetpt_truthreco_tpc.Sumw2();
  h_jetpt_reco_tpc.Sumw2();
  h_jetEta_truth_its.Sumw2();
  h_jetEta_truthreco_its.Sumw2();
  h_jetEta_reco_its.Sumw2();
  h_jetPhi_truth_its.Sumw2();
  h_jetPhi_truthreco_its.Sumw2();
  h_jetPhi_reco_its.Sumw2();
  h_jet2D_truth_its.Sumw2();
  h_jet2D_truthreco_its.Sumw2();
  h_jet2D_reco_its.Sumw2();
  h_jetRes_Phi_its.Sumw2();

 
  h_jetpt_correlation_its.SetTitle("Jet Response Matrix;p_{T}^{true};p_{T}^{reco}");
  h_jetpt_correlation2_its.SetTitle("Jet Response Matrix - p_{T}^{reco} > 5 GeV/c;p_{T}^{true};p_{T}^{reco}");
  h_jetRes_Pt_its.SetTitle("Jet Resolution Response;p_{T}^{true} [GeV/c];(p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%]");
  h_jetRes_Phi_its.SetTitle("Jet #phi Resolution;#phi^{true};#Delta#phi(#phi^{reco}-#phi^{true});counts");
  h_jetpt_correlation_tpc.SetTitle("Jet Response Matrix;p_{T}^{true};p_{T}^{reco}");
  h_jetpt_correlation2_tpc.SetTitle("Jet Response Matrix - p_{T}^{reco} > 5 GeV/c;p_{T}^{true};p_{T}^{reco}");
  h_jetRes_Pt_tpc.SetTitle("Jet Resolution Response;p_{T}^{true} [GeV/c];(p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%]");
  h_jetEta_truth_its.SetTitle(";#eta;counts");
  h_jetEta_truthreco_its.SetTitle(";#eta;counts");
  h_jetEta_reco_its.SetTitle(";#eta;counts");
  h_jetPhi_truth_its.SetTitle(";#phi;counts");
  h_jetPhi_truthreco_its.SetTitle(";#phi;counts");
  h_jetPhi_reco_its.SetTitle(";#phi;counts");
  h_jet2D_truth_its.SetTitle(";#phi;#eta");
  h_jet2D_truthreco_its.SetTitle(";#phi;#eta");
  h_jet2D_reco_its.SetTitle(";#phi;#eta");
  

  const int TrackBit = 3;//3 for TPC+ITS, 16 for ITS only
  TString ntupleName = "junk";
  TString MCname = "junk";
  int numEvents, numEvents_tracks, numEvents_clusters;
  numEvents = numEvents_tracks = numEvents_clusters = 0;
  double sumWeight = 0.0;
  TString topic = "none";
  topic = (TString)argv[1];
  cout << topic << endl;
  Float_t jetptmin = 8.0;
  //const int numMC = argc;
  //double aveXsectionArray[] = {0.0};
  
  TApplication application("", &dummyc, dummyv);
  TCanvas* canvas = new TCanvas();
  
  //Looping over ntuples
  for (int iarg = 2; iarg < argc; iarg++) {
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
    //forTracking = IsTracking(MCname);
    //cout << "\nIs this for tracking?" << forTracking << endl;


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
        
    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX];
    //Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_iso_its_04_ue[NTRACK_MAX];
    Float_t cluster_iso_04_truth[NTRACK_MAX];
    //Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
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

    //Jets reco 
    UInt_t njet_ak04its;
    Float_t jet_ak04its_pt_raw[NTRACK_MAX];
    Float_t jet_ak04its_eta_raw[NTRACK_MAX];
    Float_t jet_ak04its_phi_raw[NTRACK_MAX];    
    UInt_t njet_ak04tpc;
    Float_t jet_ak04tpc_pt_raw[NTRACK_MAX];
    Float_t jet_ak04tpc_eta_raw[NTRACK_MAX];
    Float_t jet_ak04tpc_phi_raw[NTRACK_MAX];
    
    Float_t jet_ak04its_pt_truth[NTRACK_MAX];
    Float_t jet_ak04its_eta_truth[NTRACK_MAX];
    Float_t jet_ak04its_phi_truth[NTRACK_MAX];
    Float_t jet_ak04tpc_pt_truth[NTRACK_MAX];
    Float_t jet_ak04tpc_eta_truth[NTRACK_MAX];
    Float_t jet_ak04tpc_phi_truth[NTRACK_MAX];

    Int_t   jet_ak04its_truth_index_z_reco[NTRACK_MAX][2];
    Float_t jet_ak04its_truth_z_reco[NTRACK_MAX][2];
    Float_t jet_ak04its_ptd_raw[NTRACK_MAX];
    Float_t jet_ak04its_width_sigma[NTRACK_MAX][2];
    UShort_t jet_ak04its_multiplicity[NTRACK_MAX];
    Int_t   jet_ak04tpc_truth_index_z_reco[NTRACK_MAX][2];
    Float_t jet_ak04tpc_truth_z_reco[NTRACK_MAX][2];
    Float_t jet_ak04tpc_ptd_raw[NTRACK_MAX];
    Float_t jet_ak04tpc_width_sigma[NTRACK_MAX][2];
    UShort_t jet_ak04tpc_multiplicity[NTRACK_MAX];

    //Truth Jets
    UInt_t njet_truth_ak04;
    Float_t jet_truth_ak04_pt[NTRACK_MAX];
    Float_t jet_truth_ak04_eta[NTRACK_MAX];
    Float_t jet_truth_ak04_phi[NTRACK_MAX];    
    
        
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
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    //_tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);
    _tree_event->SetBranchAddress("cluster_iso_04_truth",cluster_iso_04_truth);
    //_tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);        
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
        
    //MC
    _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
    _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);        
    _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",mc_truth_first_parent_pdg_code);
    _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
   
    //jets
    _tree_event->SetBranchAddress("njet_ak04its", &njet_ak04its);
    _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
    _tree_event->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi_raw);
    _tree_event->SetBranchAddress("jet_ak04its_pt_truth", jet_ak04its_pt_truth);
    _tree_event->SetBranchAddress("jet_ak04its_eta_truth", jet_ak04its_eta_truth);
    _tree_event->SetBranchAddress("jet_ak04its_phi_truth", jet_ak04its_phi_truth);
    _tree_event->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);
    _tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", jet_ak04tpc_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04tpc_eta_raw", jet_ak04tpc_eta_raw);
    _tree_event->SetBranchAddress("jet_ak04tpc_phi", jet_ak04tpc_phi_raw);
    _tree_event->SetBranchAddress("jet_ak04tpc_pt_truth", jet_ak04tpc_pt_truth);
    _tree_event->SetBranchAddress("jet_ak04tpc_eta_truth", jet_ak04tpc_eta_truth);
    _tree_event->SetBranchAddress("jet_ak04tpc_phi_truth", jet_ak04tpc_phi_truth);

    //quark-gluon discriminator variables
    _tree_event->SetBranchAddress("jet_ak04its_ptd_raw", jet_ak04its_ptd_raw);
    _tree_event->SetBranchAddress("jet_ak04its_width_sigma", jet_ak04its_width_sigma);
    _tree_event->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_ak04its_multiplicity);
    _tree_event->SetBranchAddress("jet_ak04tpc_ptd_raw", jet_ak04tpc_ptd_raw);
    _tree_event->SetBranchAddress("jet_ak04tpc_width_sigma", jet_ak04tpc_width_sigma);
    _tree_event->SetBranchAddress("jet_ak04tpc_multiplicity_raw", jet_ak04tpc_multiplicity);

    _tree_event->SetBranchAddress("jet_ak04its_truth_index_z_reco",     jet_ak04its_truth_index_z_reco);
    _tree_event->SetBranchAddress("jet_ak04its_truth_z_reco", jet_ak04its_truth_z_reco);    
    _tree_event->SetBranchAddress("jet_ak04tpc_truth_index_z_reco",     jet_ak04tpc_truth_index_z_reco);
    _tree_event->SetBranchAddress("jet_ak04tpc_truth_z_reco", jet_ak04tpc_truth_z_reco);    

    //truth jets
    _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04);
    _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_truth_ak04_pt);
    _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_truth_ak04_phi);
    _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_truth_ak04_eta); 

    const double maxEta = 0.8;
    Long64_t totEvents = _tree_event->GetEntries();
    Long64_t restrictEvents = 200000;
    Long64_t numEntries = TMath::Min(totEvents,restrictEvents);
    cout << numEntries << endl;
    double aveXsection = 0.0;
    // Loop over events
    for(Long64_t ievent = 0; ievent < numEntries ; ievent++){
      _tree_event->GetEntry(ievent);
      
      bool eventChange = true;
      aveXsection += (double)eg_cross_section;      

      //event selection
      if(not(TMath::Abs(primary_vertex[2])<10.0) || primary_vertex[2] == 0.0000000 ) continue; //vertex z position
      hZvertex.Fill(primary_vertex[2]);
      numEvents++;

      //Selecting pthat weights
      double weight = SetPthatWeights(ntupleName, (double)eg_cross_section, (double)eg_ntrial);
      if(ievent%10000 == 0)
	cout << weight << endl;
      sumWeight += weight;

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
	  hReco2DPtEta.Fill(track_pt[n], track_eta[n],weight);
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
	    hRecoTruth2DPtEta.Fill(mc_truth_pt[index], mc_truth_eta[index], weight);
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
	    hTruth2DPtEta.Fill(mc_truth_pt[nTru], mc_truth_eta[nTru],weight);
	  }
      }//end loop over MC particles



      eventChange = true;
      //loop over clusters
      for (ULong64_t n = 0; n < ncluster; n++) {
        //Photon Selection
	hClusterCut.Fill(0);
	h_Num2D_b.Fill(cluster_phi[n], cluster_eta[n]);
        //if( not(cluster_pt[n]>12)) continue; hClusterCut.Fill(1);                        //select pt of photons
        if( not(cluster_ncell[n]>2)) continue; hClusterCut.Fill(2);                      //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>0.05)) continue; hClusterCut.Fill(3);    //removes "spiky" clusters
        if( not(cluster_nlocal_maxima[n]< 3)) continue; hClusterCut.Fill(4);             //require to have at most 2 local maxima.
	//if( not(cluster_distance_to_bad_channel[n]>=2.0)) continue; hClusterCut.Fill(5); //require cluster to be away from a bad channel by 2 cells
	if(TMath::Abs(cluster_eta[n]) > 0.67) continue; hClusterCut.Fill(8);              //Avoiding cells on the edge    
        //Isolation and shower shape selection:                                           
        double isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n];             //remove UE subtraction
	isolation = isolation - ue_estimate_its_const*0.4*0.4*TMath::Pi();               //Use rhoxA subtraction
	//if(not (isolation < 1.0)) continue; hClusterCut.Fill(6);                         //isolation cut r= 0.4 and pt > 1
	h_NoISO.Fill(cluster_pt[n],weight);
	if(isolation < 1.0)
	  h_YesISO.Fill(cluster_pt[n],weight); 
        //if( not(cluster_lambda_square[n][0]<0.27)) continue; hClusterCut.Fill(7);        //single-photon selection (as opposed to merged photon).
	
 
	// Access the corresonding mc_truth particle; skip if index is 65535, which is invalid, or the truth particle pT is less than 10, or the mc_truth_pdg_code is not 22 (it's not a photon)
	
	Bool_t isTruePhoton = false;
        Float_t truth_pt = -999.0;
        Float_t truth_eta = -999.0;
        Float_t truth_phi = -999.0;
	for(int counter = 0 ; counter<32; counter++){
	  unsigned short index = cluster_mc_truth_index[n][counter];                   

          if(isTruePhoton) break;
          if(index==65535) continue;
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
	  //cout << "is true photon" << endl;
	  //EMCAL&DCAL recotruth hists
	  h_Num.Fill(truth_pt,weight);
	  bool inEMCalDcalTruth, inEMCalDcalReco;
	  inEMCalDcalTruth = inEMCalDcalReco = false;
	  if(1.396 < truth_phi && truth_phi < 3.264) 
	    {
	      //cout << "photon is in emcal" << endl;
	      h_Num_emcal.Fill(truth_pt,weight);
	      h_Num2D.Fill(truth_phi, truth_eta);
	      inEMCalDcalTruth = true;
	    }
	  //if(!((TMath::Abs(truth_eta[n]) < 0.2) && (-1.745 < truth_phi[n] && truth_phi[n] < -0.698)) && (-1.745 < truth_phi && truth_phi < -0.576)) //continue; hClusterCut.Fill(9);
	  if(TMath::Abs(truth_eta) > 0.2 && -1.745 < truth_phi && truth_phi < -0.576) //continue; hClusterCut.Fill(9);
	    {
	      h_Num_dcal.Fill(truth_pt,weight);
	      h_Num2D.Fill(truth_phi, truth_eta);
	      inEMCalDcalTruth = true;
	    }

	  //EMCAL&DCAL reco hists
	  h_Reco.Fill(cluster_pt[n],weight); 
	  if(1.396 < cluster_phi[n] && cluster_phi[n] < 3.264) 
	    {
	      h_Reco_emcal.Fill(cluster_pt[n],weight);
	      inEMCalDcalReco = true;
	    }  
	  //if(!((TMath::Abs(truth_eta[n]) < 0.2) && (-1.745 < truth_phi[n] && truth_phi[n] < -0.698)) && (-1.745 < truth_phi && truth_phi < -0.576)) //continue; hClusterCut.Fill(9);
	  if(TMath::Abs(cluster_eta[n]) > 0.2 && -1.745 < cluster_phi[n] && cluster_phi[n] < -0.576) //continue; hClusterCut.Fill(9);
	      //if(-1.745 < cluster_phi[n] && cluster_phi[n] < -0.576) 
	    {
	      h_Reco_dcal.Fill(cluster_pt[n], weight);
	      inEMCalDcalReco = true;
	    }

	  if(inEMCalDcalTruth && inEMCalDcalReco)
	    h_Correlation.Fill(truth_pt, cluster_pt[n],weight);
	  
	  if (eventChange) {numEvents_clusters++; eventChange = false;}
	} 
      }//end loop on clusters
       


      //loop over truth particles
      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
	if(mc_truth_pdg_code[nmc]==22 && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==22)
	  {
	    
	    hCluster_iso_04_truth.Fill(cluster_iso_04_truth[nmc]);
	    h_Den2D_b.Fill(mc_truth_phi[nmc],mc_truth_eta[nmc]);
	    
	    if(TMath::Abs(mc_truth_eta[nmc]) > 0.67) continue;
	    if(1.396 < mc_truth_phi[nmc] && mc_truth_phi[nmc] < 3.264) 
	      {
		h_Den.Fill(mc_truth_pt[nmc],weight);
		h_Den_emcal.Fill(mc_truth_pt[nmc],weight);
		h_Den2D.Fill(mc_truth_phi[nmc],mc_truth_eta[nmc]);
	      }
	    //if((TMath::Abs(mc_truth_eta[nmc]) < 0.2) && (-1.745 < mc_truth_phi[nmc] && mc_truth_phi[nmc] < -0.698)) continue;
	    //if(-1.745 < mc_truth_phi[nmc] && mc_truth_phi[nmc] < -0.576) 
	    if(TMath::Abs(mc_truth_eta[nmc]) > 0.2 && -1.745 < mc_truth_phi[nmc] && mc_truth_phi[nmc] < -0.576) //continue; hClusterCut.Fill(9);
	      {
		h_Den.Fill(mc_truth_pt[nmc],weight);
		h_Den_dcal.Fill(mc_truth_pt[nmc],weight);
		h_Den2D.Fill(mc_truth_phi[nmc],mc_truth_eta[nmc]);
	      }  
	    
	    
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

      //loop over truth jets
      for (ULong64_t ijet = 0; ijet < njet_truth_ak04; ijet++) {
	if(not (jet_truth_ak04_pt[ijet]>jetptmin)) continue;
	if(not(TMath::Abs(jet_truth_ak04_eta[ijet])<0.5)) continue;
	h_jetpt_truth.Fill(jet_truth_ak04_pt[ijet], weight);
	h_jetEta_truth_its.Fill(jet_truth_ak04_eta[ijet]);
	h_jetPhi_truth_its.Fill(jet_truth_ak04_phi[ijet]);
	h_jet2D_truth_its.Fill(jet_truth_ak04_phi[ijet], jet_truth_ak04_eta[ijet]);
      }//end loop over truth jets
      
      std::set<int> temp; //to store truth indices associated with reco jets
      std::vector<pair<int, float>> tempRecoPt; //to store truth indices associated with reco jets and recojet pt
      for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) { 
	if(not (jet_ak04its_pt_raw[ijet]>jetptmin)) continue;
	if(not (TMath::Abs(jet_ak04its_eta_raw[ijet])  <0.5 ) ) continue;
	//if((jet_ak04its_pt_raw[ijet]>jetptmin) && (TMath::Abs(jet_ak04its_eta_raw[ijet])  <0.5 ))
	h_jetpt_reco_its.Fill(jet_ak04its_pt_raw[ijet], weight);
	h_jetEta_reco_its.Fill(jet_ak04its_eta_raw[ijet]);
	h_jetPhi_reco_its.Fill(jet_ak04its_phi_raw[ijet]);
	h_jet2D_reco_its.Fill(jet_ak04its_phi_raw[ijet], jet_ak04its_eta_raw[ijet]);
	  
	//temp.insert(jet_ak04its_truth_index_z_reco[ijet][0]);
	//std::pair<int, float> recoIndexPt (jet_ak04its_truth_index_z_reco[ijet][0],jet_ak04its_pt_raw[ijet]);
	//tempRecoPt.emplace_back(jet_ak04its_truth_index_z_reco[ijet][0],jet_ak04its_pt_raw[ijet]);
	//int index = jet_ak04its_truth_index_z_reco[ijet][0];
	//if(index>0){
	//if(not(TMath::Abs(jet_truth_ak04_eta[ijet])<0.5)) continue;
	h_jetpt_truthreco_its.Fill(jet_ak04its_pt_truth[ijet],weight);
	h_jetpt_correlation_its.Fill(jet_ak04its_pt_truth[ijet],jet_ak04its_pt_raw[ijet], weight);
	h_jetEta_truthreco_its.Fill(jet_ak04its_eta_truth[ijet]);
	h_jetPhi_truthreco_its.Fill(jet_ak04its_phi_truth[ijet]);
	h_jet2D_truthreco_its.Fill(jet_ak04its_phi_truth[ijet], jet_ak04its_eta_truth[ijet]);
	//if(jet_ak04its_pt_raw[ijet] > 5.0)
	//{
	h_jetpt_correlation2_its.Fill(jet_ak04its_pt_truth[ijet], jet_ak04its_pt_raw[ijet], weight);	
	h_jetRes_Pt_its.Fill(jet_ak04its_pt_truth[ijet], 100*(jet_ak04its_pt_raw[ijet]-jet_ak04its_pt_truth[ijet])/(jet_ak04its_pt_truth[ijet]),weight);
	h_jetRes_Phi_its.Fill(jet_ak04its_phi_truth[ijet], (jet_ak04its_phi_raw[ijet]-jet_ak04its_phi_truth[ijet]));
	//}
	
      } //end loop over reco its jets
      
      /*for (ULong64_t ijet = 0; ijet < njet_ak04tpc; ijet++) { 
	//if(not (jet_ak04tpc_pt_raw[ijet]>jetptmin)) continue;
	if(not (TMath::Abs(jet_ak04tpc_eta_raw[ijet])  <0.5 ) ) continue;
	h_jetpt_reco_tpc.Fill(jet_ak04tpc_pt_raw[ijet], weight);
	//temp.insert(jet_ak04tpc_truth_index_z_reco[ijet][0]);
	//std::pair<int, float> recoIndexPt (jet_ak04tpc_truth_index_z_reco[ijet][0],jet_ak04tpc_pt_raw[ijet]);
	//tempRecoPt.emplace_back(jet_ak04tpc_truth_index_z_reco[ijet][0],jet_ak04tpc_pt_raw[ijet]);
	int index = jet_ak04tpc_truth_index_z_reco[ijet][0];
	if(index>0){
	  if(not(TMath::Abs(jet_truth_ak04_eta[index])<0.5)) continue;
          h_jetpt_truthreco_tpc.Fill(jet_truth_ak04_pt[index],weight);
	  h_jetpt_correlation_tpc.Fill(jet_truth_ak04_pt[index],jet_ak04tpc_pt_raw[ijet], weight);
	  if(jet_ak04tpc_pt_raw[ijet] > 5.0)
	    {	  
	      h_jetpt_correlation2_tpc.Fill(jet_truth_ak04_pt[index], jet_ak04tpc_pt_raw[ijet], weight);
	      h_jetRes_Pt_tpc.Fill(jet_truth_ak04_pt[index], 100*(jet_ak04tpc_pt_raw[ijet]-jet_truth_ak04_pt[index])/(jet_truth_ak04_pt[index]),weight);
	    }
	}
      } //end loop over reco tpc jets*/
    
    }//end over events
    
    cout << (TString)argv[iarg] << "\t" << "average Xsection:\t" << aveXsection << "\tNumber of events:\t" << numEntries << "\tAveXsection:\t" << aveXsection/numEntries << endl;
    
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

  for(int x = 1; x < hTruth2DPtEta.GetNbinsX()+1; x++)
    {
      for(int y = 1; y < hTruth2DPtEta.GetNbinsY()+1; y++)
	{
	  double dx , dy, content, temp, contentErr, tempErr;
	  dx = dy = content = temp = contentErr = tempErr = 0.0;

	  dx = hTruth2DPtEta.GetXaxis()->GetBinWidth(x);
	  dy = hTruth2DPtEta.GetYaxis()->GetBinWidth(y);
	  content = hTruth2DPtEta.GetBinContent(x, y);
	  contentErr = hTruth2DPtEta.GetBinError(x, y);
	  temp = content/(dx*dy);
	  hTruth2DPtEta.SetBinContent(x, y, temp);
	  tempErr = tempErr/(dx*dy);
	  hTruth2DPtEta.SetBinError(x, y, tempErr);

	  dx = hRecoTruth2DPtEta.GetXaxis()->GetBinWidth(x);
	  dy = hRecoTruth2DPtEta.GetYaxis()->GetBinWidth(y);
	  content = hRecoTruth2DPtEta.GetBinContent(x, y);
	  contentErr = hRecoTruth2DPtEta.GetBinError(x, y);
	  temp = content/(dx*dy);
	  hRecoTruth2DPtEta.SetBinContent(x, y, temp);
	  tempErr = tempErr/(dx*dy);
	  hRecoTruth2DPtEta.SetBinError(x, y, tempErr);	  
	}
    }
  
  
  for(int i = 1; i < h_Reco.GetNbinsX()+1; i++)
    {
      double dpt, content, temp, error, tempErr;
      dpt = h_Reco.GetBinWidth(i);
      content = h_Reco.GetBinContent(i);
      temp = content/(dpt*numEvents);
      h_Reco.SetBinContent(i, temp);
      error = h_Reco.GetBinError(i);
      tempErr = error/(dpt*numEvents);
      h_Reco.SetBinError(i, tempErr);

      dpt = h_Reco_emcal.GetBinWidth(i);
      content = h_Reco_emcal.GetBinContent(i);
      temp = content/(dpt*numEvents);
      h_Reco_emcal.SetBinContent(i, temp);
      error = h_Reco_emcal.GetBinError(i);
      tempErr = error/(dpt*numEvents);
      h_Reco_emcal.SetBinError(i, tempErr);

      dpt = h_Reco_dcal.GetBinWidth(i);
      content = h_Reco_dcal.GetBinContent(i);
      temp = content/(dpt*numEvents);
      h_Reco_dcal.SetBinContent(i, temp);
      error = h_Reco_dcal.GetBinError(i);
      tempErr = error/(dpt*numEvents);
      h_Reco_dcal.SetBinError(i, tempErr);

      dpt = h_YesISO.GetBinWidth(i);
      content = h_YesISO.GetBinContent(i);
      temp = content/(dpt*numEvents);
      h_YesISO.SetBinContent(i, temp);
      error = h_YesISO.GetBinError(i);
      tempErr = error/(dpt*numEvents);
      h_YesISO.SetBinError(i, tempErr);

      dpt = h_NoISO.GetBinWidth(i);
      content = h_NoISO.GetBinContent(i);
      temp = content/(dpt*numEvents);
      h_NoISO.SetBinContent(i, temp);
      error = h_NoISO.GetBinError(i);
      tempErr = error/(dpt*numEvents);
      h_NoISO.SetBinError(i, tempErr);

    }
  
  h_jetpt_truth.Scale(1/sumWeight);
  h_jetpt_truth_its.Scale(1/sumWeight);
  h_jetpt_truthreco_its.Scale(1/sumWeight);
  h_jetpt_reco_its.Scale(1/sumWeight);
  h_jetpt_correlation_its.Scale(1/sumWeight);
  h_jetpt_correlation2_its.Scale(1/sumWeight);
  h_jetRes_Pt_its.Scale(1/sumWeight);
  h_jetpt_truth_tpc.Scale(1/sumWeight);
  h_jetpt_truthreco_tpc.Scale(1/sumWeight);
  h_jetpt_reco_tpc.Scale(1/sumWeight);
  h_jetpt_correlation_tpc.Scale(1/sumWeight);
  h_jetpt_correlation2_tpc.Scale(1/sumWeight);
  h_jetRes_Pt_tpc.Scale(1/sumWeight);
  
  hTrackCut.GetXaxis()->SetBinLabel(1,"All Tracks");
  hTrackCut.GetXaxis()->SetBinLabel(2,"Track quality cut");
  hTrackCut.GetXaxis()->SetBinLabel(3,"Track #eta cut");
  hTrackCut.GetXaxis()->SetBinLabel(4,"pt cut");
  hTrackCut.GetXaxis()->SetBinLabel(5,"ITS #chi^{2} cut");
  hTrackCut.GetXaxis()->SetBinLabel(6,"ITS nCluster cut");
  hTrackCut.GetXaxis()->SetBinLabel(7,"DCAr cut");
  hTrackCut.GetXaxis()->SetBinLabel(8,"DCAz cut");

  hClusterCut.GetXaxis()->SetBinLabel(1,"All Clusters");
  hClusterCut.GetXaxis()->SetBinLabel(2,"pt cut");
  hClusterCut.GetXaxis()->SetBinLabel(3,"number of cells cut");
  hClusterCut.GetXaxis()->SetBinLabel(4,"exoticity cut");
  hClusterCut.GetXaxis()->SetBinLabel(5,"local maxima cut");
  hClusterCut.GetXaxis()->SetBinLabel(6,"distance to bad channel");
  hClusterCut.GetXaxis()->SetBinLabel(7,"isolation cut");
  hClusterCut.GetXaxis()->SetBinLabel(8,"lambda cut");
  hClusterCut.GetXaxis()->SetBinLabel(9,"#eta cut");
  hClusterCut.GetXaxis()->SetBinLabel(10,"PHOS cut");


  TF1* gausfit = new TF1("gaus","gaus", -25,25);
  gausfit->SetLineColor(kRed);
  TGraphErrors* g_mean = new TGraphErrors();
  TGraphErrors* g_sigma = new TGraphErrors();
  
  auto c1 = new TCanvas();
   
  //Study of the ITS-only track pT resolution
  //const Double_t bins[10] = {  0.1,          0.16681005,   0.27825594,   0.46415888 ,  0.77426368, 1.29154967,   2.15443469,   3.59381366,   5.9948425,   10.0};
  //TFile* fout_trackRes = new TFile(Form("TrackOutput/ResProj/projections_%s_%i_%ibins_1GeV30GeV_2M_ResCut.root", MCname.Data(), TrackBit, nbinstrack),"RECREATE");
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
    
    //h1->Write(Form("trackResProjection%i_%s_TrackBit%i", i, MCname.Data(), TrackBit));
    //if(forTracking)
      //c1->SaveAs(Form("TrackOutput/PDFOUTPUT/%s_projecting%i_TrackBit%i.C", MCname.Data(), i, TrackBit));
    }//*/
  
  //fout_trackRes->Close();
  
  g_sigma->SetTitle("Relative resolution vs p_{T} ; p_{T}^{true} [GeV]; #sigma(p_{T})/p_{T} [%]"); 
  g_mean->SetTitle("; p_{T}^{true} [GeV]; Relative bias [%]");



  TGraphErrors* g_mean_jet_its = new TGraphErrors();
  TGraphErrors* g_sigma_jet_its = new TGraphErrors();
  TGraphErrors* g_mean_jet_tpc = new TGraphErrors();
  TGraphErrors* g_sigma_jet_tpc = new TGraphErrors();

   
  //Study of the ITS-only track pT resolution
  //const Double_t bins[10] = {  0.1,          0.16681005,   0.27825594,   0.46415888 ,  0.77426368, 1.29154967,   2.15443469,   3.59381366,   5.9948425,   10.0};
  //TFile* fout_jetRes = new TFile(Form("JetOutput/ResProjections/Projections_%s_0GeV30GeV_100K_wTpc_ResCut.root", MCname.Data()),"RECREATE");
  nbins = 30;
  for(int i=0; i<nbins; i++){
    double minpt = i;
    double maxpt = i+1;
    double binwidth = maxpt-minpt;
    int minbin =  h_jetRes_Pt_its.GetXaxis()->FindBin(minpt);
    int maxbin =  h_jetRes_Pt_its.GetXaxis()->FindBin(maxpt);
    TH1D* h1 = h_jetRes_Pt_its.ProjectionY("h", minbin, maxbin);
    TH1D* h2 = h_jetRes_Pt_tpc.ProjectionY("h2", minbin, maxbin);
    
    h1->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    h1->Draw();  
    h1->GetYaxis()->SetNdivisions(5);
    h1->GetXaxis()->SetNdivisions(5);
    h1->GetYaxis()->SetTitle("counts");
    h1->Fit("gaus","N0", "", -150, 150);
    h1->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    g_sigma_jet_its->SetPoint(g_sigma_jet_its->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(2));     
    g_sigma_jet_its->SetPointError(g_sigma_jet_its->GetN()-1, binwidth/2.0, gausfit->GetParError(2));
    g_mean_jet_its->SetPoint(g_mean_jet_its->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(1));
    g_mean_jet_its->SetPointError(g_mean_jet_its->GetN()-1, binwidth/2.0, gausfit->GetParError(1)); 

    
    h2->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    h2->Draw();  
    h2->GetYaxis()->SetNdivisions(5);
    h2->GetXaxis()->SetNdivisions(5);
    h2->GetYaxis()->SetTitle("counts");
    h2->Fit("gaus","N0", "", -150, 150);
    h2->SetTitle("; (p_{T}^{reco}-p_{T}^{true})/p_{T}^{true} [%] ; counts");
    g_sigma_jet_tpc->SetPoint(g_sigma_jet_tpc->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(2));     
    g_sigma_jet_tpc->SetPointError(g_sigma_jet_tpc->GetN()-1, binwidth/2.0, gausfit->GetParError(2));
    g_mean_jet_tpc->SetPoint(g_mean_jet_tpc->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(1));
    g_mean_jet_tpc->SetPointError(g_mean_jet_tpc->GetN()-1, binwidth/2.0, gausfit->GetParError(1));
   
    //h1->Write(Form("jetResProjection%i_%s_TrackBit16", i, MCname.Data()));
    //h2->Write(Form("jetResProjection%i_%s_TrackBit3", i, MCname.Data()));
  }//*/
  //fout_jetRes->Close();
  
  g_sigma_jet_its->SetTitle("Relative resolution vs p_{T} ; p_{T}^{true} [GeV]; #sigma(p_{T})/p_{T} [%]"); 
  g_mean_jet_its->SetTitle("; p_{T}^{true} [GeV]; Relative bias [%]");
  g_sigma_jet_tpc->SetTitle("Relative resolution vs p_{T} ; p_{T}^{true} [GeV]; #sigma(p_{T})/p_{T} [%]"); 
  g_mean_jet_tpc->SetTitle("; p_{T}^{true} [GeV]; Relative bias [%]");

  
  bool makeJetFile = topic.Contains("jet");
  bool makeClusterFile = topic.Contains("cluster");
  bool makeTrackFile =  topic.Contains("track");
  bool makeXsectionFile = false;

  if(makeJetFile)
    {
      TFile* fout_jet = new TFile(Form("JetOutput/MC/%s_Jets_0GeV30GeV_200K_8GeVcut_phiRes.root", MCname.Data()),"RECREATE");
      
      h_jetpt_truth.Write("hTruth");
      h_jetpt_truthreco_its.Write("hRecoTruth_its");
      h_jetpt_reco_its.Write("hReco_its");
      h_jetpt_correlation_its.Write("hCorrelation_its");  
      h_jetpt_correlation2_its.Write("hCorrelation2_its");
      h_jetRes_Pt_its.Write("hRes_Pt_its");
      h_jetRes_Phi_its.Write("hRes_phi_its");
      h_jetEta_truth_its.Write("hTruth_eta_its");
      h_jetEta_truthreco_its.Write("hRecoTruth_eta_its");
      h_jetEta_reco_its.Write("hReco_eta_its");
      h_jetPhi_truth_its.Write("hTruth_phi_its");
      h_jetPhi_truthreco_its.Write("hRecoTruth_phi_its");
      h_jetPhi_reco_its.Write("hReco_phi_its");
      h_jet2D_truth_its.Write("hTruth_phiEta_its");
      h_jet2D_truthreco_its.Write("hRecoTruth_phiEta_its");
      h_jet2D_reco_its.Write("hReco_phiEta_its");

      /*
      TGraphAsymmErrors* eff_jet_its = new TGraphAsymmErrors(&h_jetpt_truthreco_its, &h_jetpt_truth);
      eff_jet_its->Write("Efficiency_its");
      g_sigma_jet_its->Write("g_sigma_its");
      g_mean_jet_its->Write("g_mean_its");

      TGraphAsymmErrors* eff_jet_tpc = new TGraphAsymmErrors(&h_jetpt_truthreco_tpc, &h_jetpt_truth);
      eff_jet_tpc->Write("Efficiency_tpc");
      g_sigma_jet_tpc->Write("g_sigma_tpc");
      g_mean_jet_tpc->Write("g_mean_tpc");

      h_jetpt_truth.Write("hTruth");
      h_jetpt_truthreco_tpc.Write("hRecoTruth_tpc");
      h_jetpt_reco_tpc.Write("hReco_tpc");
      h_jetpt_correlation_tpc.Write("hCorrelation_tpc");  
      h_jetpt_correlation2_tpc.Write("hCorrelation2_tpc");
      h_jetRes_Pt_tpc.Write("hRes_Pt_tpc");*/

    }
  
  if(makeClusterFile)
    {
      TFile* fout_cluster = new TFile(Form("PhotonOutput/MC/%s_5GeV60GeV_100Kevents_phosDCAlgap_ISOeff_noDisBadChn_noShwrShp_noISO.root", MCname.Data()),"RECREATE");
      
      /*TGraphAsymmErrors* eff_cluster = new TGraphAsymmErrors(&h_Num, &h_Den);
      TGraphAsymmErrors* eff_cluster_emcal = new TGraphAsymmErrors(&h_Num_emcal, &h_Den_emcal);
      TGraphAsymmErrors* eff_cluster_dcal = new TGraphAsymmErrors(&h_Num_dcal, &h_Den_dcal);
      eff_cluster->Write("Efficiency");
      eff_cluster_emcal->Write("Efficiency_emcal");
      eff_cluster_dcal->Write("Efficiency_dcal");*/
      h_Den.Write("hTruth");
      h_Num.Write("hRecoTruth");
      h_Reco.Write("hReco");
      h_Den_emcal.Write("hTruth_emcal");
      h_Num_emcal.Write("hRecoTruth_emcal");
      h_Reco_emcal.Write("hReco_emcal");
      h_Den_dcal.Write("hTruth_dcal");
      h_Num_dcal.Write("hRecoTruth_dcal");
      h_Reco_dcal.Write("hReco_dcal");
      h_YesISO.Write("h_YesISO");
      h_NoISO.Write("h_NoISO");
      hCluster_iso_04_truth.Write("hCluster_iso_04_truth");
      h_Correlation.Write();
      h_Num2D_b.Write("hEtaPhi_RecoTruth_b");
      h_Den2D_b.Write("hEtaPhi_Truth_b");
      h_Num2D.Write("hEtaPhi_RecoTruth");
      h_Den2D.Write("hEtaPhi_Truth");
      fout_cluster->Close();
    }

  if(makeTrackFile)
    {
      TFile* fout_track = new TFile(Form("TrackOutput/%s_%i_%ibins_1GeV30GeV_100Levents.root", MCname.Data(), TrackBit, nbinstrack),"RECREATE");
      
      TGraphAsymmErrors* eff = new TGraphAsymmErrors(&hRecoTruth, &hTruth);
      eff->SetTitle("; p_{T}^{true} ; #epsilon");
      eff->Write("Efficiency");
      g_sigma->Write("g_sigma");
      g_mean->Write("g_mean");

      hTruth.Write("hTruth");
      hRecoTruth.Write("hRecoTruth");
      hReco.Write("hReco");
      hTruth_eta.Write("hTruth_eta");
      hRecoTruth_eta.Write("hRecoTruth_eta");
      hReco_eta.Write("hReco_eta");
      hTruth_phi.Write("hTruth_phi");
      hRecoTruth_phi.Write("hRecoTruth_phi");
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
      hRecoTruth2DPtEta.Write("hRecoTruth_ptEta");
      hReco2DPtEta.Write("hReco_ptEta");
      hTruth2DPtEta.Write("hTruth_ptEta");

      hTrackQuality.Write();
      hTrackCut.Write();

      fout_track->Close();
    }

  if(makeXsectionFile)
    {
      
    }
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}//end main

