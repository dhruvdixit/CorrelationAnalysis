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

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

float c_eta_array[17664];
float c_phi_array[17664];


unsigned int GetSuperModule(const unsigned int n){
   unsigned int sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;


  return sm;
}

void cell_5_5(unsigned int n_5_5[], const unsigned int n,
              const unsigned int ld = 5)
{
  const unsigned int sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;
        const unsigned int nphi =
	  sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        n_5_5[0 * ld + 0] = n - 2 * nphi - 4;
        n_5_5[0 * ld + 1] = n - 2 * nphi - 2;
        n_5_5[0 * ld + 2] = n - 2 * nphi;
        n_5_5[0 * ld + 3] = n - 2 * nphi + 2;
        n_5_5[0 * ld + 4] = n - 2 * nphi + 4;
        if (n % 2 == 0) {
	  n_5_5[1 * ld + 0] = n - 3;
          n_5_5[1 * ld + 1] = n - 1;
          n_5_5[1 * ld + 2] = n + 1;
          n_5_5[1 * ld + 3] = n + 3;
          n_5_5[1 * ld + 4] = n + 5;
        }
        else {
          n_5_5[1 * ld + 0] = n - 2 * nphi - 5;
          n_5_5[1 * ld + 1] = n - 2 * nphi - 3;
          n_5_5[1 * ld + 2] = n - 2 * nphi - 1;
          n_5_5[1 * ld + 3] = n - 2 * nphi + 1;
          n_5_5[1 * ld + 4] = n - 2 * nphi + 3;
        }
        n_5_5[2 * ld + 0] = n - 4;
        n_5_5[2 * ld + 1] = n - 2;
        n_5_5[2 * ld + 2] = n;
        n_5_5[2 * ld + 3] = n + 2;
        n_5_5[2 * ld + 4] = n + 4;
        if (n % 2 == 0) {
	  n_5_5[3 * ld + 0] = n + 2 * nphi - 3;
          n_5_5[3 * ld + 1] = n + 2 * nphi - 1;
          n_5_5[3 * ld + 2] = n + 2 * nphi + 1;
          n_5_5[3 * ld + 3] = n + 2 * nphi + 3;
          n_5_5[3 * ld + 4] = n + 2 * nphi + 5;
        }
        else {
          n_5_5[3 * ld + 0] = n - 5;
          n_5_5[3 * ld + 1] = n - 3;
          n_5_5[3 * ld + 2] = n - 1;
          n_5_5[3 * ld + 3] = n + 1;
          n_5_5[3 * ld + 4] = n + 3;
        }
        n_5_5[4 * ld + 0] = n + 2 * nphi - 4;
        n_5_5[4 * ld + 1] = n + 2 * nphi - 2;
        n_5_5[4 * ld + 2] = n + 2 * nphi;
        n_5_5[4 * ld + 3] = n + 2 * nphi + 2;
        n_5_5[4 * ld + 4] = n + 2 * nphi + 4;
}


int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    int dummyc = 1;
    char **dummyv = new char *[1];

    dummyv[0] = strdup("main");

    for (int iarg = 1; iarg < argc; iarg++) {
      std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);

        if (file == NULL) {
	  std::cout << " fail" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
	std::cout<<"About to get TTree" << std::endl;


        TTree *_tree_event = NULL;
        _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
	  _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));

	  if (_tree_event == NULL) {
	      std::cout << " fail " << std::endl;
	      exit(EXIT_FAILURE);
	      //}
        }  
        //_tree_event->Print();


   	UInt_t ncluster;
        UInt_t cluster_nmc_truth[NTRACK_MAX];
        Float_t cluster_e[NTRACK_MAX];
        Float_t cluster_pt[NTRACK_MAX];
        Float_t cluster_eta[NTRACK_MAX];
        Float_t cluster_phi[NTRACK_MAX];
          
	UShort_t  cluster_cell_id_max[NTRACK_MAX];

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

	//new branches:
        Float_t cluster_NN1[NTRACK_MAX];
        Float_t cluster_NN2[NTRACK_MAX];
        Float_t cluster_Lambda[NTRACK_MAX];
        Float_t cluster_isHardPhoton[NTRACK_MAX];
        Float_t cluster_minMass[NTRACK_MAX];
        Int_t cluster_pi0tagged[NTRACK_MAX];
        ////
	Double_t primary_vertex[3];

        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];

        Float_t cluster_lambda_square[NTRACK_MAX][2];   
    
	_tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
	_tree_event->SetBranchAddress("cluster_eta", cluster_eta);
	_tree_event->SetBranchAddress("cluster_phi", cluster_phi);
	_tree_event->SetBranchAddress("cluster_e", cluster_e);
        _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
        _tree_event->SetBranchAddress("cluster_lambda_square",  cluster_lambda_square);
        _tree_event->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth);
	_tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);

	Float_t cell_e[17664];
        Float_t cell_eta[17664];
        Float_t cell_phi[17664];

        _tree_event->SetBranchAddress("cell_e", cell_e);
        _tree_event->SetBranchAddress("cell_eta", cell_eta);
        _tree_event->SetBranchAddress("cell_phi", cell_phi);

	_tree_event->SetBranchAddress("nmc_truth",&nmc_truth);
        _tree_event->SetBranchAddress("mc_truth_pt",mc_truth_pt);
        _tree_event->SetBranchAddress("mc_truth_eta",mc_truth_eta);
        _tree_event->SetBranchAddress("mc_truth_phi",mc_truth_phi);
        _tree_event->SetBranchAddress("mc_truth_charge",mc_truth_charge);
        _tree_event->SetBranchAddress("mc_truth_pdg_code",mc_truth_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code", mc_truth_first_parent_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_first_parent_e", mc_truth_first_parent_e);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pt", mc_truth_first_parent_pt);
        _tree_event->SetBranchAddress("mc_truth_first_parent_eta", mc_truth_first_parent_eta);
        _tree_event->SetBranchAddress("mc_truth_first_parent_phi", mc_truth_first_parent_phi);
	_tree_event->SetBranchAddress("mc_truth_status",  mc_truth_status);

        _tree_event->SetBranchStatus("*",0);
        _tree_event->SetBranchStatus("*cluster*",1);
	_tree_event->SetBranchStatus("*cell*",1);
	_tree_event->SetBranchStatus("track*",1);
        _tree_event->SetBranchStatus("*mc*", 1);
        _tree_event->SetBranchStatus("*vertex*",1);
        _tree_event->SetBranchStatus("*trigger*",1);
        _tree_event->SetBranchStatus("*period*",1);
        _tree_event->SetBranchStatus("*bunch*",1);
        _tree_event->SetBranchStatus("*number*",1);
        _tree_event->SetBranchStatus("*time*",1);
        _tree_event->SetBranchStatus("*pileup*",1);
	_tree_event->SetBranchStatus("*event*",1);

        _tree_event->SetBranchStatus("*centrality*",1);
	_tree_event->SetBranchStatus("*multiplicity*",1);
        _tree_event->SetBranchStatus("run_number",1);
        _tree_event->SetBranchStatus("*Mix*",1);
	_tree_event->SetBranchStatus("*jet*",1);
	_tree_event->SetBranchStatus("*muon*",0);


 	std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

	TFile *newfile = new TFile("13dsmall.root","recreate");
	TTree *newtree = _tree_event->CloneTree(0);
        newtree->Branch("cluster_NN1", cluster_NN1, "cluster_NN1[ncluster]/F");
        newtree->Branch("cluster_NN2", cluster_NN2, "cluster_NN2[ncluster]/F");

	newtree->Branch("cluster_Lambda", cluster_Lambda, "cluster_Lambda[ncluster]/F");
        newtree->Branch("cluster_isHardPhoton", cluster_isHardPhoton, "cluster_isHardPhoton[ncluster]/F");
        newtree->Branch("cluster_minMass", cluster_minMass, "cluster_minMass[ncluster]/F");
        newtree->Branch("cluster_pi0tagged",cluster_pi0tagged,"cluster_pi0tagged[ncluster]/I");

	Float_t cluster_b5x5[NTRACK_MAX];
	Float_t cluster_b5x5_lin[NTRACK_MAX];
	unsigned int cluster_SuperModule[NTRACK_MAX];
	newtree->Branch("cluster_b5x5", cluster_b5x5, "cluster_b5x5[ncluster]/F");
        newtree->Branch("cluster_b5x5_lin", cluster_b5x5_lin, "cluster_b5x5_lin[ncluster]/F");
        newtree->Branch("cluster_SuperModule", cluster_SuperModule, "cluster_SuperModule[ncluster]/i");
       
	_tree_event->GetEntry(0);
        for(int i=0; i<17664; i++){
          c_eta_array[i] = cell_eta[i];
          c_phi_array[i] = cell_phi[i]; 
        }
        const Long64_t nevents = _tree_event->GetEntries();
	//const Long64_t nevents = 500;  
	for(Long64_t ievent = 0; ievent < nevents ; ievent++){     
	  _tree_event->GetEntry(ievent);
          
          bool KeepEvent = false;
	  for (ULong64_t n = 0; n < ncluster; n++) {
            Float_t IsTrueHardPhoton = -666.0;
	    cluster_NN1[n] = cluster_s_nphoton[n][1]; 
            cluster_NN2[n] = cluster_s_nphoton[n][2];
            cluster_Lambda[n] = cluster_lambda_square[n][0];
	    for(UInt_t counter = 0 ; counter<cluster_nmc_truth[n]; counter++)
	      {
                unsigned short index = cluster_mc_truth_index[n][counter];
                if(index!=65535){
		  if(mc_truth_pdg_code[index]==22 && mc_truth_first_parent_pdg_code[index]==22) IsTrueHardPhoton = 1.0;
                }
	      }//end loop over indices
            cluster_isHardPhoton[n] = IsTrueHardPhoton;
             //pi0 veto: cluster_minMass
            float minmass = 999;
            Int_t pi0tag = 0;

            TLorentzVector v1;
	    v1.SetPtEtaPhiM(cluster_pt[n], cluster_eta[n], cluster_phi[n], 0.0);
	    for (ULong64_t m = 0; m < ncluster; m++) {
	      if(cluster_pt[m]<0.7) continue;
	      if(m==n) continue;
	      TLorentzVector v2;
	      v2.SetPtEtaPhiM(cluster_pt[m], cluster_eta[m],   cluster_phi[m], 0.0);
	      TLorentzVector pi0 = v1+v2;
              if(pi0.M()<minmass) minmass = pi0.M();
              if(pi0.M()<0.160 and pi0.M()>0.100) pi0tag=1;
	    }
	    cluster_pi0tagged[n] =pi0tag; 
	    cluster_minMass[n] = minmass;
            cluster_SuperModule[n] = GetSuperModule(cluster_cell_id_max[n]);

	    if(cluster_pt[n]>10.0){
              KeepEvent = true;
            }

	    unsigned int cell_id_5_5[25];
	    cell_5_5(cell_id_5_5, cluster_cell_id_max[n]);

            
	    double eta_center = 0;
	    double phi_center = 0;
	    double sumw        = 0;
	    for (size_t i = 0; i < 25; i++) {
	      auto c_id = cell_id_5_5[i];
	      if( not(c_id <17664)) continue;
	      if( not(cell_e[c_id] >0.100)) continue;
	      eta_center =+ cell_e[c_id]* c_eta_array[c_id];
	      phi_center =+ cell_e[c_id]* c_phi_array[c_id];
	      sumw       = + cell_e[c_id];
	    }
	    eta_center = eta_center/sumw;
	    phi_center = phi_center/sumw;
	    float ce  = 0;
	    float cep = 0;
	    float cp  = 0;
	    float m   = 0;
	    for(size_t i=0; i<25; i++){
	      auto c_id = cell_id_5_5[i];
              if( not(c_id <17664)) continue;
	      if( not(cell_e[c_id] >0.100)) continue;
              float w = TMath::Log(cell_e[c_id]/cluster_e[n]);
	      ce  += (c_eta_array[c_id] - eta_center) * (c_eta_array[c_id]- eta_center) * w;
	      cep += (c_eta_array[c_id] - eta_center) * (c_phi_array[c_id] - phi_center) * w;
	      cp  += (c_phi_array[c_id] - phi_center) * (c_phi_array[c_id] - phi_center) * w;
	      m   += w;
	    }
	    ce /= m;
	    cep /= m;
	    cp /= m;
	    cluster_b5x5[n] = 1000*0.5 * (ce + cp + TMath::Sqrt( TMath::Power(ce - cp, 2.0) + 4*TMath::Power(cep,2.0)));



	    ce  = 0;
            cep = 0;
            cp  = 0;
            m   = 0;
            for(size_t i=0; i<25; i++){
              auto c_id = cell_id_5_5[i];
              if( not(c_id <17664)) continue;
              if( not(cell_e[c_id] >0.100)) continue;
              float w = cell_e[c_id];
              ce  += (c_eta_array[c_id] - eta_center) * (c_eta_array[c_id]- eta_center) * w;
              cep += (c_eta_array[c_id] - eta_center) * (c_phi_array[c_id] - phi_center) * w;
              cp  += (c_phi_array[c_id] - phi_center) * (c_phi_array[c_id] - phi_center) * w;
              m   += w;
            }
            ce /= m;
            cep /= m;
            cp /= m;
            cluster_b5x5_lin[n] = 1000*0.5 * (ce + cp + TMath::Sqrt( TMath::Power(ce - cp, 2.0) + 4*TMath::Power(cep,2.0)));
	  
	  }//end looop clusters 
	  if (ievent % 25000 == 0) { std::cout << " event " << ievent << std::endl;}
	  if(ievent==0){
	    newtree->Fill();
	    //std::cout <<" cell e, eta , phi " << cell_e[0] << " " << cell_eta[0] << std::endl;
	  }
	  if(TMath::Abs(primary_vertex[2])>10) continue;
	  if(KeepEvent){
            newtree->Fill();
          }
	}

	//newtree->Print();
	std::cout << "#events passing skimming: " << newtree->GetEntries() << std::endl;
	newtree->AutoSave();
	delete file;
	delete newfile;
	  
    }//end loop over samples

    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
