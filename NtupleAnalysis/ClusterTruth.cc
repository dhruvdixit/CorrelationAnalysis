
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

bool ptDepShowerShapeCut(Float_t clus_pt, Float_t lambda2)
{
  //bool passShowerShape = true;

  if(lambda2 < 0.1)
    return false;
  
  if((clus_pt > 10 && clus_pt < 12) && lambda2 > 0.40)
    return false;

  if((clus_pt > 12 && clus_pt < 14) && lambda2 > 0.35)
    return false;

  if((clus_pt > 14 && clus_pt < 16) && lambda2 > 0.32)
    return false;

  if((clus_pt > 16 && clus_pt < 60) && lambda2 > 0.30)
    return false;


  return true;
}

Float_t Get_Purity_ErrFunction(Float_t pT_GeV, std::string deviation = "") {

  Float_t purity_val = 0;

  //Non-platue assumption
  // Float_t par[3] = {0.548247710,
  //                   8.794543375,
  //                   12.7423900};

  //Old
//   Float_t par[3] = {0.54225742923,
//                     8.09242373515,
//                     11.8085154181};

  Float_t par[3] = {0.539253098581,
		    8.84942587038,
		    12.360736025};

  if (strcmp(deviation.data(),"Plus")==0){
    par[0] = 0.60750016509;
    par[1] = 7.05184155403;
    par[2] = 13.6116163603;
  }

  if (strcmp(deviation.data(),"Minus")==0){
    par[0] = 0.479958593235;
    par[1] = 9.05392932723;
    par[2] = 10.2061359452;
  }


  purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
  return purity_val;
}

double SetPthatWeights(TString MCname, double Xsection, double ntrial)
{

  //General purpose MC need no weights
  if(MCname(0,4) == "13b2" || MCname(0,4) == "17l4" || MCname(0,4) == "17l3" || MCname(0,4) == "16k5")
    return 1;
   
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


int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  //int dummyc = 1;
  //char **dummyv = new char *[1];
    
  //dummyv[0] = strdup("main");

  gStyle->SetOptStat("");
  
  TString topic = "none";
  //topic = (TString)argv[1];
  //cut << topic << endl;

  bool doJets = topic.Contains("jet");
  bool doClusters = true;//topic.Contains("cluster");
  bool doTracks =  topic.Contains("track");

  TH1F hZvertex("hZvertex","",60, -30, 30);
  hZvertex.Sumw2();  
  hZvertex.SetTitle(";Z_{v} [cm]; counts");
  

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
    
  const int nbinscluster = 14;
  //Double_t clusterbins[nbinscluster+1] = {12.0 , 14.06198608, 16.47828771, 19.30978769, 22.62783047, 26.51601976, 31.07232506, 36.4115502 , 42.668226  , 50.0};//geom binning
  Double_t clusterbins[nbinscluster+1] = {5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 12.00, 14.00, 16.00, 18.00, 20.00, 25.00, 30.00, 40.00, 60.00};//nbinscluster = 14, Erwann binning
  //Double_t clusterbins[nbinscluster+1] = {};
  

  TH1D h_Den("h_Den", "", nbinscluster, clusterbins);
  TH1D h_Den_emcal("h_Den_emcal", "", nbinscluster, clusterbins);
  TH1D h_Den_dcal("h_Den_dcal", "", nbinscluster, clusterbins);
  TH1D h_Num("h_Num", "", nbinscluster, clusterbins);
  TH1D h_Num_emcal("h_Num_emcal", "", nbinscluster, clusterbins);
  TH1D h_Num_dcal("h_Num_dcal", "", nbinscluster, clusterbins);
  TH1D h_Reco("h_Reco", "", nbinscluster, clusterbins);
  TH1D h_RecoAll("h_RecoAll", "", nbinscluster, clusterbins);
  TH1D h_Reco_emcal("h_Reco_emcal", "", nbinscluster, clusterbins);
  TH1D h_Reco_dcal("h_Reco_dcal", "", nbinscluster, clusterbins);
  TH1D h_YesISO("h_YesISO", "", nbinscluster, clusterbins);
  TH1D h_NoISO("h_NoISO", "", nbinscluster, clusterbins);
  TH1D h_YesSSC("h_YesSSC", "", nbinscluster, clusterbins);
  TH1D h_NoSSC("h_NoSSC", "", nbinscluster, clusterbins);
  TH1D h_Efficiency("h_Efficiency", "", nbinscluster, clusterbins);
  TH1D h_migEfficiency("h_migEfficiency", "", nbinscluster, clusterbins);
  TH1D h_recoEfficiency("h_recoEfficiency", "", nbinscluster, clusterbins);
  TH1D h_ISOEfficiency("h_ISOEfficiency", "", nbinscluster, clusterbins);
  TH1D h_SSCEfficiency("h_SSCEfficiency", "", nbinscluster, clusterbins);
  TH1D h_Fake("h_Fake", "", nbinscluster, clusterbins);
  TH1D h_FakeRate("h_FakeRate", "", nbinscluster, clusterbins);
  TH1D hClusterCut("hClusterCut", "", 20, -0.5, 19.5);
  TH1D hCluster_iso_04_truth("hCluster_iso_04_truth", "", 20 , 0, 18);
  
  TH2D h_Correlation("h_Correlation", "", nbinscluster, clusterbins, nbinscluster, clusterbins);
  TH2D h_Correlation_corr("h_Correlation_corr", "", nbinscluster, clusterbins, nbinscluster, clusterbins);
  //TH2D h_Num2D("h_Num2D","", nbinsphi, phibins, nbinseta, etabins);
  //TH2D h_Den2D("h_Den2D","", nbinsphi, phibins, nbinseta, etabins);
  
  TH2D h_Num2D("h_Num2D","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  TH2D h_Den2D("h_Den2D","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  TH2D h_Num2D_b("h_Num2D_b","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  TH2D h_Den2D_b("h_Den2D_b","", 416, -1.9, 3.3, 96, -0.8, 0.8);
  
  
  /*////////////////////////////////////
    
    
  ////////////////////////////////////*/

  TH1D* hClusterCutFlow[10];
  for(int i = 0; i < 10; i++)
    hClusterCutFlow[i] = new TH1D(Form("hClusterCutFlow%i", i), "", 20, -0.5, 19.5);
  

  h_Den.Sumw2();
  h_Num.Sumw2();
  h_Den_emcal.Sumw2();
  h_Num_emcal.Sumw2();
  h_Den_dcal.Sumw2();
  h_Num_dcal.Sumw2();
  h_Reco.Sumw2();
  h_RecoAll.Sumw2();
  h_Reco_emcal.Sumw2();
  h_Reco_dcal.Sumw2();
  h_YesISO.Sumw2();
  h_NoISO.Sumw2();
  h_YesSSC.Sumw2();
  h_NoSSC.Sumw2();
  h_Fake.Sumw2();
  h_Efficiency.Sumw2();
  h_migEfficiency.Sumw2();
  h_recoEfficiency.Sumw2();
  h_ISOEfficiency.Sumw2();
  h_SSCEfficiency.Sumw2();
  h_FakeRate.Sumw2();
  h_Correlation.Sumw2();
  h_Correlation_corr.Sumw2();
  
  h_Den.SetTitle("truth photons; p_{T}^{truth} [GeV/c]; entries");
  h_Num.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV/c]; entries");
  h_Den_emcal.SetTitle("truth photons; p_{T}^{truth} [GeV/c]; entries");
  h_Num_emcal.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV/c]; entries");
  h_Den_dcal.SetTitle("truth photons; p_{T}^{truth} [GeV/c]; entries");
  h_Num_dcal.SetTitle("reco photons filled with truthpt reco; p_{T}^{Reco,truthpt} [GeV/c]; entries");
  h_Reco.SetTitle("true reco photons filled with pt reco; p_{T}^{Reco,pt} [GeV/c]; 1/N_{event} dN/dp_{T}");
  h_RecoAll.SetTitle("reco photons filled with pt reco; p_{T}^{Reco,pt} [GeV/c]; 1/N_{event} dN/dp_{T}");
  h_Reco_emcal.SetTitle("reco clustes filled with pt reco; p_{T}^{Reco,pt} [GeV/c];1/N_{event} dN/dp_{T}");
  h_Reco_dcal.SetTitle("reco clustes filled with pt reco; p_{T}^{Reco,pt} [GeV/c];1/N_{event} dN/dp_{T}");
  h_Fake.SetTitle("fake photons filled with pt true; p_{T}^{Reco,truthpt} [GeV/c]; 1/N_{event} dN/dp_{T}");
  h_YesISO.SetTitle(";p_{T}^{Reco,pt} [GeV/c];  dN/dp_{T}");
  h_NoISO.SetTitle(";p_{T}^{Reco,pt} [GeV/c];  dN/dp_{T}");
  h_YesSSC.SetTitle(";p_{T}^{Reco,pt} [GeV/c];  dN/dp_{T}");
  h_NoSSC.SetTitle(";p_{T}^{Reco,pt} [GeV/c];  dN/dp_{T}");
  h_Efficiency.SetTitle("Efficiency = #epsilon_{reco} x #epsilon_{iso} x #epsilon_{ssc}; E_{T} [GeV]; #epsilon");
  h_recoEfficiency.SetTitle("Reconstruction Efficiency; E_{T} [GeV]; #epsilon_{reco}");
  h_migEfficiency.SetTitle("E_{T} Bin Migration Efficiency (recoE/trueE); E_{T} [GeV]; #epsilon_{mig}");
  h_ISOEfficiency.SetTitle("Isolation Efficiency; E_{T} [GeV]; #epsilon_{iso}");
  h_SSCEfficiency.SetTitle("Single Photon Efficiency; E_{T} [GeV]; #epsilon_{ssc}");
  h_FakeRate.SetTitle("Fake Rate; E_{T} [GeV]; #epsilon_{fake}");
  
  h_Correlation.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
  h_Correlation_corr.SetTitle("; True p_{T} [GeV/c]; Reconstructed p_{T} [GeV/c]");
  h_Num2D.SetTitle(";#phi_{true}; #eta_{true}");
  h_Den2D.SetTitle(";#phi_{true}; #eta_{true}");
  h_Num2D_b.SetTitle(";#phi_{true}; #eta_{true}");
  h_Den2D_b.SetTitle(";#phi_{true}; #eta_{true}");
  
  hClusterCut.GetXaxis()->SetBinLabel(1,"All Clusters");
  hClusterCut.GetXaxis()->SetBinLabel(2,"ncell>=2");
  hClusterCut.GetXaxis()->SetBinLabel(3,"exoticity>0.03");
  hClusterCut.GetXaxis()->SetBinLabel(4,"num local maxima <=2");
  hClusterCut.GetXaxis()->SetBinLabel(5,"distance to bad channel>2");
  hClusterCut.GetXaxis()->SetBinLabel(6,"isolation");
  hClusterCut.GetXaxis()->SetBinLabel(7,"shower shape");
  hClusterCut.GetXaxis()->SetBinLabel(8,"#eta");
  hClusterCut.GetXaxis()->SetBinLabel(9,"#phi");
  hClusterCut.GetXaxis()->SetBinLabel(10,"accepted cluster");
  
  TString ntupleName = "junk";
  TString MCname = "junk";
  int numEvents, numEvents_tracks, numEvents_clusters;
  numEvents = numEvents_tracks = numEvents_clusters = 0;
  double sumWeight = 0.0;
  
  Float_t jetptmin = 8.0;
  //const int numMC = argc;
  //double aveXsectionArray[] = {0.0};
  
  //TApplication application("", &dummyc, dummyv);
  //TCanvas* canvas = new TCanvas();
  
  //Looping over ntuples
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile* file = TFile::Open((TString)argv[iarg]);
        
    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
    TString temp = (TString)argv[iarg];
    ntupleName = temp(temp.Last('/')+1,temp.First('_')-temp.Last('/')+6);
    MCname = ntupleName(0,ntupleName.Length()-7);
    //cout << MCname << endl;

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
    //}//end loop over ntuples
            
    //you define variables
    Double_t primary_vertex[3];
    Float_t ue_estimate_its_const;
    Float_t eg_cross_section;
    Int_t   eg_ntrial;

    UInt_t ntrack;
        
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

            
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
    _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
   
    //tracks
    _tree_event->SetBranchAddress("ntrack", &ntrack);
        
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
    

    const double maxEta = 0.52;
    Long64_t totEvents = _tree_event->GetEntries();
    cout << totEvents << endl;
    Long64_t restrictEvents = 50000;
    Long64_t numEntries = TMath::Min(totEvents,restrictEvents);
    cout << numEntries << endl;
    double aveXsection = 0.0;
    // Loop over events
    for(Long64_t ievent = 0; ievent < numEntries ; ievent++){
      _tree_event->GetEntry(ievent);
      
      if(ievent%100000 == 0)
	cout << "Event #" << ievent << endl;
      
      bool eventChange = true;
      aveXsection += (double)eg_cross_section;      

      //Event Selection
      if(not(TMath::Abs(primary_vertex[2])<10.0)) continue; //vertex z position
      if(primary_vertex[2] == 0.000000) continue;
      //if(ntrack < 0) continue;
      hZvertex.Fill(primary_vertex[2]);
      numEvents++;

      //Selecting pthat weights
      double weight = SetPthatWeights(ntupleName, (double)eg_cross_section, (double)eg_ntrial);
      if(ievent%10000 == 0)
	{
	  cout << weight << endl;
	  cout << ntupleName.Data() << endl;
	  cout << ievent << endl;
	}
      sumWeight += weight;


      eventChange = true;
      //loop over clusters
      //if(doClusters)
	{
	  for (ULong64_t n = 0; n < ncluster; n++) {
	    //Photon Selection
	    hClusterCut.Fill(0);
	    h_Num2D_b.Fill(cluster_phi[n], cluster_eta[n]);
	    
	    double isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n];             //remove UE subtraction
	    isolation = isolation - ue_estimate_its_const*0.4*0.4*TMath::Pi();               //Use rhoxA subtraction
	    bool isoPassed = false;
	    
	    //Photon Selection
	    ULong64_t clusterCutBits = 0;
	    ULong64_t clusterCutPassed = 0;
	    //if( not(cluster_pt[n]>8)) {continue;} //select pt of photons
	    if( (cluster_ncell[n]>=2))
	      clusterCutBits |= (1 << 0); 
	    clusterCutPassed |= (1 << 0); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(1);//removes clusters with 1 or 2 cells
	    if( (cluster_e_cross[n]/cluster_e[n]>0.05))                          
	      clusterCutBits |= (1 << 1); 
	    clusterCutPassed |= (1 << 1); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(2);//removes "spiky" clusters
	    if( (cluster_nlocal_maxima[n]<= 2))                                  
	      clusterCutBits |= (1 << 2); 
	    clusterCutPassed |= (1 << 2); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(3);//require to have at most 2 local maxima.
	    if( (cluster_distance_to_bad_channel[n]>=1.0)) 
	      clusterCutBits |= (1 << 3); 
	    clusterCutPassed |= (1 << 3); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(4);//distnace to bad channels
	    
	    //Eta+Phi cuts
	    if(TMath::Abs(cluster_eta[n] < 0.67)) 
	      clusterCutBits |= (1 << 6); 
	    clusterCutPassed |= (1 << 6); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(7);//eta cut
	    if(1.396 < cluster_phi[n] && cluster_phi[n] < 3.28) 
	      clusterCutBits |= (1 << 7); 
	    clusterCutPassed |= (1 << 7); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(8);//phi cut*/
	    
	    //Isolation and shower shape selection:
	    //if( (ptDepShowerShapeCut(cluster_pt[n], cluster_lambda_square[n][0])))
	    if(( 0.1 < cluster_lambda_square[n][0]) &&  ( 0.3 > cluster_lambda_square[n][0]))
	      h_YesSSC.Fill(cluster_pt[n],weight);//clusterCutBits |= (1 << 5); 
	    //clusterCutPassed |= (1 << 5); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(6); //single-photon selection, not merged
	    h_NoSSC.Fill(cluster_pt[n],weight);
	    if( (isolation < 1.5))
	      h_YesISO.Fill(cluster_pt[n],weight);//clusterCutBits |= (1 << 4); 
	    //clusterCutPassed |= (1 << 4); if(clusterCutBits == clusterCutPassed) hClusterCut.Fill(5);//isolation r = 0.4 and energy < 2
	    h_NoISO.Fill(cluster_pt[n],weight);

	    if(clusterCutBits != clusterCutPassed) { continue;}
	    //cout << "cluster rejected" << endl;
	    
	    
	    h_RecoAll.Fill(cluster_pt[n], weight);

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
	    
	    if(!isTruePhoton)
	      h_Fake.Fill(cluster_pt[n], weight);

	    if(isTruePhoton){
	      //fill in this histogram only photons that can be traced to a generated non-decay photon.	
	      //cout << "is true photon" << endl;

	      h_Num.Fill(truth_pt,weight);
	      h_Reco.Fill(cluster_pt[n], weight);
	      h_Correlation.Fill(truth_pt, cluster_pt[n], weight);
	      h_Correlation_corr.Fill(truth_pt, cluster_pt[n], weight);
	      
	      if (eventChange) {numEvents_clusters++; eventChange = false;}
	    } 
	  }//end loop on clusters
	  
	  
	  
	  //loop over truth particles for clusters
	  for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {

	    if(not (TMath::Abs(mc_truth_eta[nmc] < 0.667))) continue;
	    if(not (1.396 < mc_truth_phi[nmc] && mc_truth_phi[nmc] < 3.28)) continue;
	    if(mc_truth_pdg_code[nmc]==22 && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==22)
	      {
		
		hCluster_iso_04_truth.Fill(cluster_iso_04_truth[nmc]);
		h_Den2D_b.Fill(mc_truth_phi[nmc],mc_truth_eta[nmc]);
		h_Den.Fill(mc_truth_pt[nmc],weight);		
		
	      }
	  } //end loop over mc truth particle
	}//end doCluster
    }//end over events
    
    cout << (TString)argv[iarg] << "\t" << "average Xsection:\t" << aveXsection << "\tNumber of events:\t" << numEvents << "\tAveXsection:\t" << aveXsection/numEntries << endl;
    
  }//end loop over ntuples
  cout << numEvents_tracks << endl;
  cout << numEvents_clusters << endl;
  cout << sumWeight << endl;

  
  for(int i = 1; i < h_Den.GetNbinsX()+1; i++)
    {
      double dEt = h_Den.GetBinWidth(i);
      cout << dEt << "\t" << endl;
      double contentDen = h_Den.GetBinContent(i);
      double errorDen = h_Den.GetBinError(i);
      double newContDen = contentDen/dEt;
      double newErrDen = errorDen/dEt;
      h_Den.SetBinContent(i, newContDen);
      h_Den.SetBinError(i, newErrDen);

      double contentNum = h_Num.GetBinContent(i);
      double errorNum = h_Num.GetBinError(i);
      double newContNum = contentNum/dEt;
      double newErrNum = errorNum/dEt;
      h_Num.SetBinContent(i, newContNum);
      h_Num.SetBinError(i, newErrNum);

      double contentReco = h_Reco.GetBinContent(i);
      double errorReco = h_Reco.GetBinError(i);
      double newContReco = contentReco/dEt;
      double newErrReco = errorReco/dEt;
      h_Reco.SetBinContent(i, newContReco);
      h_Reco.SetBinError(i, newErrReco);

      
      double contentRecoAll = h_RecoAll.GetBinContent(i);
      double errorRecoAll = h_RecoAll.GetBinError(i);
      double newContRecoAll = contentRecoAll/dEt;
      double newErrRecoAll = errorRecoAll/dEt;
      h_Reco.SetBinContent(i, newContRecoAll);
      h_Reco.SetBinError(i, newErrRecoAll);

      double contentFake = h_Fake.GetBinContent(i);
      double errorFake = h_Fake.GetBinError(i);
      double newContFake = contentFake/dEt;
      double newErrFake = errorFake/dEt;
      h_Fake.SetBinContent(i, newContFake);
      h_Fake.SetBinError(i, newErrFake);

      double contentYesISO = h_YesISO.GetBinContent(i);
      double errorYesISO = h_YesISO.GetBinError(i);
      double newContYesISO = contentYesISO/dEt;
      double newErrYesISO = errorYesISO/dEt;
      h_YesISO.SetBinContent(i, newContYesISO);
      h_YesISO.SetBinError(i, newErrYesISO);
	
      double contentNoISO = h_NoISO.GetBinContent(i);
      double errorNoISO = h_NoISO.GetBinError(i);
      double newContNoISO = contentNoISO/dEt;
      double newErrNoISO = errorNoISO/dEt;
      h_NoISO.SetBinContent(i, newContNoISO);
      h_NoISO.SetBinError(i, newErrNoISO);

      double contentYesSSC = h_YesSSC.GetBinContent(i);
      double errorYesSSC = h_YesSSC.GetBinError(i);
      double newContYesSSC = contentYesSSC/dEt;
      double newErrYesSSC = errorYesSSC/dEt;
      h_YesSSC.SetBinContent(i, newContYesSSC);
      h_YesSSC.SetBinError(i, newErrYesSSC);
	
      double contentNoSSC = h_NoSSC.GetBinContent(i);
      double errorNoSSC = h_NoSSC.GetBinError(i);
      double newContNoSSC = contentNoSSC/dEt;
      double newErrNoSSC = errorNoSSC/dEt;
      h_NoSSC.SetBinContent(i, newContNoSSC);
      h_NoSSC.SetBinError(i, newErrNoSSC);


      h_recoEfficiency.SetBinContent(i, newContNum/newContDen);
      h_recoEfficiency.SetBinError(i, TMath::Sqrt((TMath::Power(newErrNum, 2))+(TMath::Power(newErrDen,2))));

      h_ISOEfficiency.SetBinContent(i, newContYesISO/newContNoISO);
      h_ISOEfficiency.SetBinError(i, TMath::Sqrt((TMath::Power(newErrYesISO, 2))+(TMath::Power(newErrNoISO,2))));

      h_SSCEfficiency.SetBinContent(i, newContYesSSC/newContNoSSC);
      h_SSCEfficiency.SetBinError(i, TMath::Sqrt((TMath::Power(newErrYesSSC, 2))+(TMath::Power(newErrNoSSC,2))));

      h_migEfficiency.SetBinContent(i, newContReco/newContNum);
      h_migEfficiency.SetBinError(i, TMath::Sqrt((TMath::Power(newErrReco, 2))+(TMath::Power(newErrNum,2))));

      h_FakeRate.SetBinContent(i, newContFake/newContRecoAll);
      h_FakeRate.SetBinError(i, TMath::Sqrt((TMath::Power(newErrFake, 2))+(TMath::Power(newErrRecoAll,2))));


      double recoEffRelError = h_recoEfficiency.GetBinError(i)/h_recoEfficiency.GetBinContent(i); 
      double isoEffRelError = h_ISOEfficiency.GetBinError(i)/h_ISOEfficiency.GetBinContent(i); 
      double sscEffRelError = h_SSCEfficiency.GetBinError(i)/h_SSCEfficiency.GetBinContent(i); 
      double effRelError = TMath::Sqrt(TMath::Power(recoEffRelError, 2)+TMath::Power(isoEffRelError, 2)+TMath::Power(sscEffRelError, 2));
      h_Efficiency.SetBinContent(i, h_recoEfficiency.GetBinContent(i)*h_ISOEfficiency.GetBinContent(i)*h_SSCEfficiency.GetBinContent(i));
      h_Efficiency.SetBinError(i, h_Efficiency.GetBinContent(i)*effRelError);
    }
  

  for(int x = 1; x < h_Correlation_corr.GetNbinsX()+1; x++){
    for(int y = 1; y < h_Correlation_corr.GetNbinsY()+1; y++){
      double dx , dy, content, temp, contentErr, tempErr;
      dx = dy = content = temp = contentErr = tempErr = 0.0;
      dx = h_Correlation_corr.GetXaxis()->GetBinWidth(x);
      dy = h_Correlation_corr.GetYaxis()->GetBinWidth(y);
      content = h_Correlation_corr.GetBinContent(x, y);
      contentErr = h_Correlation_corr.GetBinError(x, y);
      
      temp = content/(dx*dy);
      h_Correlation_corr.SetBinContent(x, y, temp);
      tempErr = tempErr/(dx*dy);
      h_Correlation_corr.SetBinError(x, y, tempErr);
      
    }
  }
  


  TH1D eventSelection("eventSelection","", 8, -0.5, 7.5);
  TH1D normalizer("normalizer", "normalizer", 10, -0.5, 9.5);
  
  const double tot_eta = 1.04;  
  
  normalizer.SetBinContent(1,numEvents);
  normalizer.SetBinContent(2,tot_eta);

  
  
  TFile* fout_cluster = new TFile(Form("PhotonOutput/MC/%s_50Kevents_erwanbinning_noNorm_newIsoDef_fullEMCalAcceptance_isoEffsscEffRecoEffSeperate_cutFlow_binMigration_hCorr_corr_stdCuts.root", MCname.Data()),"RECREATE");
  
  h_Den.Write("hTruth");
  h_Num.Write("hRecoTruth");
  h_Reco.Write("hReco");

  h_RecoAll.Write("hRecoAll");
  h_Fake.Write("hFake");
  h_FakeRate.Write("hFakeRate");

  h_Den_emcal.Write("hTruth_emcal");
  h_Num_emcal.Write("hRecoTruth_emcal");
  h_Reco_emcal.Write("hReco_emcal");
  h_Den_dcal.Write("hTruth_dcal");
  h_Num_dcal.Write("hRecoTruth_dcal");
  h_Reco_dcal.Write("hReco_dcal");

  h_Efficiency.Write("hEfficiency");
  h_recoEfficiency.Write("hRecoEfficiency");
  h_migEfficiency.Write("hMigEfficiency");

  h_YesISO.Write("h_YesISO");
  h_NoISO.Write("h_NoISO");
  h_ISOEfficiency.Write("hISOEfficiency");

  h_YesSSC.Write("h_YesSSC");
  h_NoSSC.Write("h_NoSSC");
  h_SSCEfficiency.Write("hSSCEfficiency");

  hCluster_iso_04_truth.Write("hCluster_iso_04_truth");
  h_Correlation.Write("hCorrelation");
  h_Correlation_corr.Write("hCorrelation_dNdEt");
  h_Correlation_corr.Scale(1/sumWeight);
  h_Correlation_corr.Write("hCorrelation_dNdEt_weightScale");
 
  h_Num2D_b.Write("hEtaPhi_RecoTruth_b");
  h_Den2D_b.Write("hEtaPhi_Truth_b");
  h_Num2D.Write("hEtaPhi_RecoTruth");
  h_Den2D.Write("hEtaPhi_Truth");
  hZvertex.Write("hZvertex");
  hClusterCut.Write("hClusterCut");
  normalizer.Write("normalizer");
  fout_cluster->Close();//*/

  cout << "finish writing to file" << endl;//*/

  return EXIT_SUCCESS;
}//end main

