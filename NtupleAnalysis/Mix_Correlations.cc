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
#include "H5Cpp.h"

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

using namespace H5;
int main(int argc, char *argv[])
{
    if (argc < 3) {
      std::cout<<"Temporary Syntax for Mixing is [Command] [mixed_root_file] [hdf5_file] (repeat for mult. samples)"<<std::endl;
      exit(EXIT_FAILURE);
    }

    int dummyc = 1;
    char **dummyv = new char *[1];

    dummyv[0] = strdup("main");

    for (int iarg = 1; iarg < argc; iarg+=2) {
      //FIXME: Current sytanx is [root] [hdf5], [root] [hdf5], (repeat)
      // int iarg = 1;
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);

        if (file == NULL) {
	  std::cout << " fail" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();

        TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

	if (_tree_event == NULL) {
	  std::cout << " fail " << std::endl;
	  exit(EXIT_FAILURE);
        }  
        //_tree_event->Print();

	//TApplication application("", &dummyc, dummyv);
	
   
        TCanvas canvas("canvas", "");
	TH1D histogram0("histogram0", "", 16, 8.0, 16.0);
        //TH2D histogram1("histogram1", "", 30, -1.5, 1.5, 18, -0.5, 1.5);
        //TH1D histogram2("histogram2", "", 18, -0.5,1.5);
        TH1D histogram3("histogram3", "", 18, -0.5,1.5);
        TH1D h_ntrig("h_ntrig", "", 2, -0.5,1.0);
         
        TH1D htrack_pt("htrack_pt","", 100, 0.0,15.0);
        TH1D htrack_eta("htrack_eta","",100, -0.8, 0.8);
        TH1D htrack_phi("htrack_phi","",100, -TMath::Pi(), TMath::Pi());

        TH1D hcluster_pt("hcluster_pt","", 100, 0.0,15.0);
	TH1D hcluster_eta("hcluster_eta","",100, -0.8, 0.8);
        TH1D hcluster_phi("hcluster_phi","",100, -TMath::Pi(), TMath::Pi());

        TH1D hcluster_sumiso("hcluster_sumiso","", 100, -5.0, 25.0);
	TH1D hcluster_sumisoNoUE("hcluster_sumisoNoUE","", 100, -5.0, 25.0);
 
	const int nztbins = 7;
        const float ztbins[nztbins+1] = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};

	TH1F* h_dPhi_iso[nztbins];
        TH1F* h_dPhi_noniso[nztbins];
        TH1F* h_dPhiLargedEta_iso[nztbins];
        TH1F* h_dPhiLargedEta_noniso[nztbins];

	TH2D Corr = TH2D("Correlation", "GS Mixed #gamma-H [all] Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7);
        Corr.Sumw2();
        Corr.SetMinimum(0.);

	TH2D IsoCorr = TH2D("Iso_Correlation", "GS Mixed #gamma-H [Iso] Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7);
	IsoCorr.Sumw2();
	IsoCorr.SetMinimum(0.);

	TH2D AntiIsoCorr = TH2D("Anti_Iso_Correlation", "GS Mixed #gamma-H [AntiIso] Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7);
        AntiIsoCorr.Sumw2();
        AntiIsoCorr.SetMinimum(0.);

        const int nphibins = 18;
	for (int izt = 0; izt<nztbins; izt++){
          h_dPhi_iso[izt] = new TH1F(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"",nphibins,-0.5,1.5);
          h_dPhi_noniso[izt] = new TH1F(Form("dPhi_noniso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", nphibins,-0.5,1.5);
	  h_dPhiLargedEta_iso[izt] = new TH1F(Form("dPhiLargedEta_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"",nphibins,-0.5,1.5);
          h_dPhiLargedEta_noniso[izt] = new TH1F(Form("dPhiLargedEta_noniso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", nphibins,-0.5,1.5);

          h_dPhi_iso[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
          h_dPhi_noniso[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
	  h_dPhiLargedEta_iso[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
          h_dPhiLargedEta_noniso[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
    
          h_dPhi_iso[izt]->Sumw2();
          h_dPhi_noniso[izt]->Sumw2();
	  h_dPhiLargedEta_iso[izt]->Sumw2();
          h_dPhiLargedEta_noniso[izt]->Sumw2();
	}

        //variables 
        Double_t primary_vertex[3];
        UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
	Float_t track_eta_emcal[NTRACK_MAX];
	Float_t track_phi_emcal[NTRACK_MAX];
	UChar_t track_quality[NTRACK_MAX];

	UInt_t ncluster;
        Float_t cluster_e[NTRACK_MAX];
	Float_t cluster_e_cross[NTRACK_MAX];
        Float_t cluster_pt[NTRACK_MAX];
        Float_t cluster_eta[NTRACK_MAX];
        Float_t cluster_phi[NTRACK_MAX]; 
        Float_t cluster_iso_tpc_04[NTRACK_MAX];
        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        Float_t cluster_lambda_square[NTRACK_MAX][2];   
        Float_t cell_e[17664];

	Long64_t Mix_Events[10];

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
 
	_tree_event->SetBranchStatus("*mc*", 0);
  
        _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
	_tree_event->SetBranchAddress("track_eta_emcal",track_eta_emcal);
	_tree_event->SetBranchAddress("track_phi_emcal",track_phi_emcal);
        _tree_event->SetBranchAddress("track_quality", track_quality);

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

 	_tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        _tree_event->SetBranchAddress("cell_e", cell_e);

	_tree_event->SetBranchAddress("Mix_Events", Mix_Events);
	//_tree_event->SetBranchAddress("LimitUse_Mixed_Events", Mix_Events);
	 
 	
 	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
    
        const float isomax = 2.0; //2 GeV is maximum energy
        const float nonisomin = 5.0;

// 	UInt_t ntrack_max = 517; 
// 	UInt_t ncluster_max = 125;

	UInt_t ntrack_max = 0; 
	UInt_t ncluster_max = 0;
	
	fprintf(stderr, "\r%s:%d: %s\n", __FILE__, __LINE__, "Determining ntrack_max and ncluster_max needed for hdf5 hyperslab");
	for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
	  _tree_event->GetEntry(i);
	  ntrack_max = std::max(ntrack_max, ntrack);
	  ncluster_max = std::max(ncluster_max, ncluster);
	  fprintf(stderr, "\r%s:%d: %llu", __FILE__, __LINE__, i);
	}
	ncluster_max = 23;
	fprintf(stderr, "\n%s:%d: maximum tracks:%i maximum clusters:%i\n", __FILE__, __LINE__, ntrack_max,ncluster_max);


	//FIXME:: when mixing with a nother data set, or hdf5 file from other data set, make sure ntrackmax is in bounds


	//open hdf5

	const H5std_string hdf5_file_name(argv[iarg+1]);
	//FIXME: Can parse the string s.t. remove .root and simply add .hdf5
	//FIXME: This will obviously require stringent naming conventions

	const H5std_string track_ds_name( "track" );
	H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY );
	DataSet track_dataset = h5_file.openDataSet( track_ds_name );
	DataSpace track_dataspace = track_dataset.getSpace();

	const H5std_string cluster_ds_name( "cluster" );
	DataSet cluster_dataset = h5_file.openDataSet( cluster_ds_name );
	DataSpace cluster_dataspace = cluster_dataset.getSpace();

	//Define array hyperslab will be read into

	float track_data_out[1][ntrack_max][7];
	float cluster_data_out[1][ncluster_max][5];

	//Define hyperslab size and offset in  FILE;
	hsize_t track_offset[3] = {0, 0, 0};
	hsize_t track_count[3] = {1, ntrack_max, 7};
	hsize_t cluster_offset[3] = {0, 0, 0};
	hsize_t cluster_count[3] = {1, ncluster_max, 5};

	track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
	cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
	fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

	//Define the memory dataspace to place hyperslab
	const int RANK_OUT = 3;
	hsize_t track_dimsm[3] = {1, ntrack_max, 7};      //memory space dimensions
	DataSpace track_memspace( RANK_OUT, track_dimsm );
	hsize_t cluster_dimsm[3] = {1, ncluster_max, 5};      //memory space dimensions
	DataSpace cluster_memspace( RANK_OUT, cluster_dimsm );

	//Define memory offset for hypreslab->array
	hsize_t track_offset_out[3] = {0};
	hsize_t cluster_offset_out[3] = {0};

	//define 2D memory hyperslab
	hsize_t track_count_out[3] = {1, ntrack_max, 7};    // size of the hyperslab in memory
	hsize_t cluster_count_out[3] = {1, ncluster_max, 5};

	//define space in memory for hyperslab, then write from file to memory
	track_memspace.selectHyperslab( H5S_SELECT_SET, track_count_out, track_offset_out );
	track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
	fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "track dataset read into array: OK");

	cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
	cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
	fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "cluster dataset read into array: OK");

	Long64_t nentries = _tree_event->GetEntries();
	int clusters_passed = 0;
	for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){     
	  //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
             _tree_event->GetEntry(ievent);
	     fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

	      for (ULong64_t n = 0; n < ncluster; n++) {
	          if( not(cluster_pt[n]>8 and cluster_pt[n]<16)) continue; //select pt of photons
                  if( not(cluster_s_nphoton[n][1]>0.75 and cluster_s_nphoton[n][1]<0.85)) continue; //select deep-photons
                  if( not(TMath::Abs(cluster_eta[n])<0.6)) continue; //cut edges of detector
                  if( not(cluster_ncell[n]>2)) continue;   //removes clusters with 1 or 2 cells 
                  if( not(cluster_e_cross[n]/cluster_e[n]>0.03)) continue; //removes "spiky" clusters
		  //electron veto
                  float SumIso =0.0;
                  bool hasMatch = false;

		  const int TrackCutBit =16;
		  for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
		    if((track_quality[itrack]&TrackCutBit)==0) continue; //select only tracks that pass selection 16
		    if( not(track_pt[itrack]>0.15)) continue;
		    Float_t deta =  cluster_eta[n]-track_eta[itrack];
		    Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack]);
		    Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
		    //Calculate the isolation energy    
		    if(dR<0.4)  SumIso = SumIso + track_pt[itrack];
		    //Veto charged-particles
		    if(not (track_pt[itrack]>5)) continue; //only veto if track has > 5 GeV.
		    if(TMath::Abs(deta)<0.015) hasMatch = true; //removes clusters near a high-pt track
		    if(TMath::Abs(dphi)<0.015) hasMatch = true; //removes clusters neat a high-pt track
		  }
		  
		  if(hasMatch) continue; //remove clusters near a high-pt track
		  //Cluster selection ended
		  hcluster_sumiso.Fill(SumIso);
		  hcluster_sumisoNoUE.Fill(cluster_iso_tpc_04[n]);
		  hcluster_pt.Fill(cluster_pt[n]);
		  hcluster_eta.Fill(cluster_eta[n]);
		  hcluster_phi.Fill(cluster_phi[n]);
		  
		  const float Isolation = cluster_iso_tpc_04[n];
		  //const float Isolation = SumIso;
		  if(Isolation<isomax){
		    histogram0.Fill(cluster_pt[n]); //isolated deep-photon pt spectra
		    h_ntrig.Fill(0);
		  }
		  if(Isolation>nonisomin){ 
		    h_ntrig.Fill(0.5);
		  }
		  
		  for (Long64_t imix = 0; imix < 10; imix++){
		    Long64_t mix_event = Mix_Events[imix];	
		    if (mix_event == ievent) continue;
		    	    
		    //read mixed event into hyperslab Place in track loop if high pt cut.
		    track_offset[0]=mix_event;
		    track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
		    track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );

		    cluster_offset[0]=mix_event;
		    cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
		    cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
		    
		    for (ULong64_t itrack = 0; itrack < ntrack_max; itrack++) {            
		      //FIXED was looping nutil ntrack, a property of original, not mixed, event. Needed isnan check!
		      //select only tracks that pass selection 3
		      if (std::isnan(track_data_out[0][itrack][1])) break;
		      if ((int(track_data_out[0][itrack][4]+0.5)&TrackCutBit)==0) continue;		      
		      if (track_data_out[0][itrack][1] < 0.15) continue;
                      htrack_pt.Fill(track_data_out[0][itrack][1]);
                      htrack_eta.Fill(track_data_out[0][itrack][2]);
                      htrack_phi.Fill(track_data_out[0][itrack][3]);
		      
		      //veto charged particles from mixed event tracks
		      bool Mix_HasMatch = false;
		      for (unsigned int l = 0; l < ncluster_max; l++){
			if (std::isnan(cluster_data_out[0][l][0])) break;
			if (TMath::Abs(cluster_data_out[0][l][2] - track_data_out[0][itrack][5]) < 0.015) {
			  Mix_HasMatch = true;
			  break; }
			if (TMath::Abs(cluster_data_out[0][l][3] - track_data_out[0][itrack][6]) < 0.015) {
			  Mix_HasMatch = true;
			  break; }
		      }
		      if (Mix_HasMatch) continue;
		      //FIXME: See if need to skip when emcal_eta = -999, or track_emcal_phi = 0.026....

		      //FIXME: Lazy implementation from past code. Will use this repositories ∆'s soon
		      Float_t DeltaPhi = cluster_phi[n] - track_data_out[0][itrack][3];
		      if (DeltaPhi < -M_PI/2){DeltaPhi += 2*M_PI;}  //if less then -pi/2 add 2pi
		      if (DeltaPhi > 3*M_PI/2){DeltaPhi =DeltaPhi -2*M_PI;}
		      Float_t DeltaEta = cluster_eta[n] - track_data_out[0][itrack][2];
		      if ((TMath::Abs(DeltaPhi) < 0.01) && (TMath::Abs(DeltaEta) < 0.01)) continue;
		      
		      Corr.Fill(DeltaPhi,DeltaEta);
		      if(Isolation<isomax) IsoCorr.Fill(DeltaPhi,DeltaEta);
		      if(Isolation>nonisomin) AntiIsoCorr.Fill(DeltaPhi,DeltaEta);

		      Double_t zt = track_pt[itrack]/cluster_pt[n];
		      Float_t deta =  cluster_eta[n]-track_data_out[0][itrack][2];;
		      Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_data_out[0][itrack][3]);
                      dphi = dphi/TMath::Pi();
		      //if(!(TMath::Abs(deta)<0.6)) continue;
                      if(dphi<-0.5) dphi +=2;
		      
		      for(int izt = 0; izt<nztbins ; izt++){
                        if(zt>ztbins[izt] and  zt<ztbins[izt+1])
			  {
			    if(Isolation< isomax){
			      if(TMath::Abs(deta)<0.6) h_dPhi_iso[izt]->Fill(dphi);
			      else if(TMath::Abs(deta)>0.8 and TMath::Abs(deta)<1.4) h_dPhiLargedEta_iso[izt]->Fill(dphi);
			    }
                	    else if(Isolation> nonisomin){
			      if(TMath::Abs(deta)<0.6) h_dPhi_noniso[izt]->Fill(dphi);
			      else if(TMath::Abs(deta)>0.8 and TMath::Abs(deta)<1.4) h_dPhiLargedEta_noniso[izt]->Fill(dphi);
			    }
			  } 
		      } //end loop over zt bins
		    }//end loop over tracks
		  }//end loop over mixed events
		  if(Isolation< isomax) clusters_passed += 1;
	      }//end loop on clusters. 
	    if (ievent % 25000 == 0) {
	       hcluster_sumiso.Draw("e1x0");
               hcluster_sumisoNoUE.SetLineColor(2);
               hcluster_sumisoNoUE.Draw("e1x0same");
               canvas.Update();
               //std::cout << "Event # " << ievent << " / " << _tree_event->GetEntries() << std::endl;
            }
	} //end loop over events
	std::cout<<clusters_passed<<std::endl;
	TFile* mix_fout = new TFile("mix_fout.root","RECREATE");
	histogram0.Write("DeepPhotonSpectra");
        h_ntrig.Write("ntriggers");
   
	for (int izt = 0; izt<nztbins; izt++){
	  h_dPhi_iso[izt]->SetMinimum(0.0);
	  h_dPhi_iso[izt]->Write();
	  h_dPhi_noniso[izt]->SetMinimum(0.0);
	  h_dPhi_noniso[izt]->Write();	  
          h_dPhiLargedEta_iso[izt]->SetMinimum(0.0);
          h_dPhiLargedEta_iso[izt]->Write();
          h_dPhiLargedEta_noniso[izt]->SetMinimum(0.0);
          h_dPhiLargedEta_noniso[izt]->Write();
	}
        htrack_pt.Write("track_pt");
        htrack_eta.Write("track_eta");
        htrack_phi.Write("track_phi");

        hcluster_pt.Write("cluster_pt");
        hcluster_eta.Write("cluster_eta");
        hcluster_phi.Write("cluster_phi");

        hcluster_sumiso.Write("cluster_sumiso");
	hcluster_sumisoNoUE.Write("cluster_tpc04iso");
	mix_fout->Close();     
	  
	TString NewFileName = "GS_gamma_hadron.root";
	TFile *MyFile = new TFile(NewFileName,"RECREATE");
	Corr.Write();
	IsoCorr.Write();
	AntiIsoCorr.Write();
	MyFile->Print();

    }//end loop over samples

    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}