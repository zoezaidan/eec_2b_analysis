//Create the response matrix for 1D, 2D and 3D unfolding and possibly does some systematic uncertainty estimation
//(that part of the code was given to me by Lida, I did not use it). The functions create_response_nD for n=1,2,3 
//are basically the same, except that 1D does not have pair matching. The most commented function is 
//create_response_3D, so if something is not clear in the other two, you can check there :) 
//There are also two functions to check the dr and eec bin migration (to decide whether unfolding is needed also
//for the dr and eec)
//
//As of last spring the latest version of RooUnfold has a bug, you should be using an older version in order for it
//to work. What I do to use that version is write the following commands in the terminal (works only for the polui 
//machines):
// source /cvmfs/cms.cern.ch/cmsset_default.sh
// cd CMSSW_10_6_30_patch1/
// cmsenv
// source /data_CMS/cms/meuli/forZoe/Unfolding/setup.sh
//then you are able to use ROOT as usual, but some things don't work in it, so I only use it whenever I strictly 
//need any of the RooUnfold versions.

#include "binning_histos_small.h"
#include <iostream>
#include <vector>
#include <cstdlib> // For std::abs
#include <cmath>    // For std::abs (integer types)
#include "TFile.h"
#include "tTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

//Skips MC events with too large event weight
bool skipMC(double jtpt, double refpt, double pthat) {
    if (!(refpt>0)) return true;    
    if (pthat<0.35*jtpt) return true;
    return false;
}

//Set of functions given by Lida but not used, needed for systematics/some tests

//For 1D unfolding
void fill_jk_resampling(std::vector<TH1D *> histos, double num, double x, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,w);
    }
}

void fill_jk_resampling_response(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco,  
                                 double x_gen, 
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,x_gen,w);
    }
}

//For 2D unfolding
void fill_jk_resampling(std::vector<TH2D *> histos, double num, double x, double y, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,y,w);
    }
}

void fill_jk_resampling_response(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco, double y_reco,
                                 double x_gen, double y_gen,
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,y_reco,x_gen,y_gen,w);
    }
}

//For 3D unfolding
void fill_jk_resampling(std::vector<TH3D *> histos, double num, double x, double y, double z, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,y,z,w);
    }
}

void fill_jk_resampling_response(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco, double y_reco, double z_reco,
                                 double x_gen, double y_gen, double z_gen,
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,y_reco,z_reco,x_gen,y_gen,z_gen,w);
    }
}

//___________________________________________________________________________________________________
//_________________________________Create response matrices__________________________________________
//___________________________________________________________________________________________________

//Creates a 1D response matrix and purity/efficiency corrections from a tree
void create_response_1D(TString &filename,  TString &dataset, TString &label, TString &folder, bool &btag, Int_t n, Float_t &pT_low, Float_t &pT_high)
{    //"trees_" + label + "_" + dataset + "_" + pT_selection + "gen.root";   

    TString flav = label; 
    std::cout << "flav:" << flav << std::endl;

    //Create the fout name depending on the selection
    TString fout_name = "histos_response_1D_";

    if(!btag) label += "_notag"; 

    fout_name += TString(Form("n%i_", n)) + dataset + "_" + label + ".root";

    //Builds the full file path and prints which file is being processed
    TString fullpath = folder + filename;
    std::cout << "Processing input file: " << fullpath << std::endl;


    std::cout << "fin: " << folder+filename << std::endl;
    TFile *fin = new TFile(folder+filename);

    TString tree_name = "tree";
    std::cout << "tree: " << tree_name << std::endl;
    TTree *tree = (TTree *) fin->Get(tree_name);

    
    Int_t ndr_reco, ndr_gen, jtHadFlav, jtNbHad;
    Float_t weight, pthat, jpt_reco, jpt_gen, mb_reco, mb_gen, jt_eta_reco, jt_eta_gen, discr;
    Float_t dr_reco[4000], dr_gen[4000], eec_reco[4000], eec_gen[4000];
   

    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("pthat", &pthat);
    tree->SetBranchAddress("ndr_reco", &ndr_reco);
    tree->SetBranchAddress("ndr_gen", &ndr_gen);
    tree->SetBranchAddress("jpt_reco", &jpt_reco);
    tree->SetBranchAddress("jpt_gen", &jpt_gen);
    tree->SetBranchAddress("mB_reco", &mb_reco);
    tree->SetBranchAddress("mB_gen", &mb_gen);
    tree->SetBranchAddress("dr_reco", &dr_reco);
    tree->SetBranchAddress("dr_gen", &dr_gen);
    tree->SetBranchAddress("eec_reco", &eec_reco);
    tree->SetBranchAddress("eec_gen", &eec_gen);
    tree->SetBranchAddress("jt_eta_reco", &jt_eta_reco);
    tree->SetBranchAddress("jt_eta_gen", &jt_eta_gen);
    tree->SetBranchAddress("discr", &discr);
    tree->SetBranchAddress("jtHadFlav", &jtHadFlav);
    tree->SetBranchAddress("jtNbHad", &jtNbHad);


    Int_t dr_bins = bins_dr;

    // random number generator for jackknife resampling (Lida)
    const double range_from = 0;
    const double range_to = 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(range_from, range_to);

    
    // Declare histograms
    TH2D *h_half0_purity_numerator_eecpt = new TH2D("h_half0_purity_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half0_purity_denominator_eecpt = new TH2D("h_half0_purity_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half0_efficiency_numerator_eecpt = new TH2D("h_half0_efficiency_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half0_efficiency_denominator_eecpt = new TH2D("h_half0_efficiency_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_purity_numerator_eecpt = new TH2D("h_half1_purity_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_purity_denominator_eecpt = new TH2D("h_half1_purity_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_efficiency_numerator_eecpt = new TH2D("h_half1_efficiency_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_efficiency_denominator_eecpt = new TH2D("h_half1_efficiency_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);

    RooUnfoldResponse *response_half0_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half0", "response for 2d: eec and jet pt"); 
    RooUnfoldResponse *response_half1_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half1", "response for 2d: eec and jet pt"); 
    RooUnfoldResponse *response_full_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_full_eecpt", "response for 2d: eec and jet pt"); 

    // rg jk resampling histogram definiton (Lida)
    TH2D *h0_purity_numerator_eecpt, 
        *h1_purity_numerator_eecpt, 
        *h2_purity_numerator_eecpt, 
        *h3_purity_numerator_eecpt, 
        *h4_purity_numerator_eecpt, 
        *h5_purity_numerator_eecpt, 
        *h6_purity_numerator_eecpt, 
        *h7_purity_numerator_eecpt, 
        *h8_purity_numerator_eecpt, 
        *h9_purity_numerator_eecpt;

    std::vector<TH2D *> histos_purity_numerator_eecpt = {
        h0_purity_numerator_eecpt, 
        h1_purity_numerator_eecpt, 
        h2_purity_numerator_eecpt, 
        h3_purity_numerator_eecpt, 
        h4_purity_numerator_eecpt, 
        h5_purity_numerator_eecpt, 
        h6_purity_numerator_eecpt, 
        h7_purity_numerator_eecpt, 
        h8_purity_numerator_eecpt, 
        h9_purity_numerator_eecpt,
    };

    TH2D *h0_purity_denominator_eecpt, 
        *h1_purity_denominator_eecpt, 
        *h2_purity_denominator_eecpt, 
        *h3_purity_denominator_eecpt, 
        *h4_purity_denominator_eecpt, 
        *h5_purity_denominator_eecpt, 
        *h6_purity_denominator_eecpt, 
        *h7_purity_denominator_eecpt, 
        *h8_purity_denominator_eecpt, 
        *h9_purity_denominator_eecpt;

    std::vector<TH2D *> histos_purity_denominator_eecpt = {
        h0_purity_denominator_eecpt, 
        h1_purity_denominator_eecpt, 
        h2_purity_denominator_eecpt, 
        h3_purity_denominator_eecpt, 
        h4_purity_denominator_eecpt, 
        h5_purity_denominator_eecpt, 
        h6_purity_denominator_eecpt, 
        h7_purity_denominator_eecpt, 
        h8_purity_denominator_eecpt, 
        h9_purity_denominator_eecpt,
    };

    TH2D *h0_efficiency_numerator_eecpt, 
        *h1_efficiency_numerator_eecpt, 
        *h2_efficiency_numerator_eecpt, 
        *h3_efficiency_numerator_eecpt, 
        *h4_efficiency_numerator_eecpt, 
        *h5_efficiency_numerator_eecpt, 
        *h6_efficiency_numerator_eecpt, 
        *h7_efficiency_numerator_eecpt, 
        *h8_efficiency_numerator_eecpt, 
        *h9_efficiency_numerator_eecpt;

    std::vector<TH2D *> histos_efficiency_numerator_eecpt = {
        h0_efficiency_numerator_eecpt, 
        h1_efficiency_numerator_eecpt, 
        h2_efficiency_numerator_eecpt, 
        h3_efficiency_numerator_eecpt, 
        h4_efficiency_numerator_eecpt, 
        h5_efficiency_numerator_eecpt, 
        h6_efficiency_numerator_eecpt, 
        h7_efficiency_numerator_eecpt, 
        h8_efficiency_numerator_eecpt, 
        h9_efficiency_numerator_eecpt,
    };

    TH2D *h0_efficiency_denominator_eecpt, 
        *h1_efficiency_denominator_eecpt, 
        *h2_efficiency_denominator_eecpt, 
        *h3_efficiency_denominator_eecpt, 
        *h4_efficiency_denominator_eecpt, 
        *h5_efficiency_denominator_eecpt, 
        *h6_efficiency_denominator_eecpt, 
        *h7_efficiency_denominator_eecpt, 
        *h8_efficiency_denominator_eecpt, 
        *h9_efficiency_denominator_eecpt;

    std::vector<TH2D *> histos_efficiency_denominator_eecpt = {
        h0_efficiency_denominator_eecpt, 
        h1_efficiency_denominator_eecpt, 
        h2_efficiency_denominator_eecpt, 
        h3_efficiency_denominator_eecpt, 
        h4_efficiency_denominator_eecpt, 
        h5_efficiency_denominator_eecpt, 
        h6_efficiency_denominator_eecpt, 
        h7_efficiency_denominator_eecpt, 
        h8_efficiency_denominator_eecpt, 
        h9_efficiency_denominator_eecpt,
    };

    RooUnfoldResponse *response0_eecpt, 
        *response1_eecpt,
        *response2_eecpt,
        *response3_eecpt,
        *response4_eecpt,
        *response5_eecpt,
        *response6_eecpt,
        *response7_eecpt,
        *response8_eecpt,
        *response9_eecpt;

    std::vector<RooUnfoldResponse *> responses_eecpt = {
        response0_eecpt, 
        response1_eecpt,
        response2_eecpt,
        response3_eecpt,
        response4_eecpt,
        response5_eecpt,
        response6_eecpt,
        response7_eecpt,
        response8_eecpt,
        response9_eecpt,
    };
    

    // initialize histograms and responses
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_purity_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
        histos_purity_denominator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_purity_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_numerator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_efficiency_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_denominator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_efficiency_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);

        responses_eecpt[i] = new RooUnfoldResponse(histos_purity_denominator_eecpt[i],histos_efficiency_denominator_eecpt[i], "response"+TString(Form("%d",i))+"_eecpt", "response for 2d: eec and jet pt"); 
    }        


    // Loop over tree entries
    Long64_t nentries = tree->GetEntries();
    std::cout << "Entries: " << nentries << std::endl;
    for (Long64_t ient = 0; ient < nentries; ient++) {
        //Print progress
        if (ient%1000000==0) std::cout << "ient=" << ient << std::endl; 
        tree->GetEntry(ient);

   

        Int_t cuts = 1;
        if      (flav == "b1")    cuts = 1;
        else if (flav == "b2")    cuts = 2;
        else if (flav == "nonb")  cuts = 3;
        else if (flav == "all")   cuts = 4;
        else if (flav == "c")     cuts = 5;
        else if (flav == "light") cuts = 6;
        else {
        std::cerr << "ERROR: Unknown label '" << label << "'\n";
        exit(1);
            }

        bool skip = false;

        // Select jet flavour and/or select on the number of b hadrons
        switch(cuts){
            //b-jet with one b hadron
            case 1: 
            if (std::abs(jtHadFlav) < 5) skip = true;
            if (std::abs(jtNbHad) != 1) skip = true;
            break;
            //b-jet with more than 1 b hadron
            case 2:
            if (std::abs(jtHadFlav) < 5) skip = true;
            if (std::abs(jtNbHad) < 2) skip = true;
            break;
            //non-b jets
            case 3:
            if (std::abs(jtHadFlav) == 5) skip = true;
            break;
            //no flavour selection
            case 4:
            skip = false;
            break;
            //c-jets
            case 5:
            if(std::abs(jtHadFlav) != 4) skip = true;
            break;
            //light (non-b non-c) jets
            case 6:
            if(std::abs(jtHadFlav) >= 4) skip = true;
            break;
        }

        if (skip) continue;
    

        if (skipMC(jpt_reco, jpt_gen, pthat)) continue;

        // Random number for jack-knife resampling (Lida)
        double num = distr(generator);

        // --------------------------

        // Check if pass cuts
        bool has_gen_match = (jpt_gen > 0);

        // Fill histograms
        if (!has_gen_match) {   
            // fill fakes
            continue; 
        } 

        
        // Skip jets outside tracker 
        if (std::abs(jt_eta_reco) > 1.6) continue;
        if (std::abs(jt_eta_gen) > 1.6) continue;

        //Select jets passing the reco b-jet tagging 
        if (btag && std::abs(discr) <= 0.99) continue;


        
        // The rest of the histograms don;t include any fakes

        
        //Loop over all dr (no matching needed) to fill the eec histograms
        for (Int_t j = 0; j < ndr_reco; j++){

            Float_t dr_reco_j = dr_reco[j];

            Float_t eec_reco_j = eec_reco[j];



            num = distr(generator);

            
            //checks for underflow/overflow
            if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
            if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;

            bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min);
            bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min);


            // fill eec histograms (both gen and reco dr and eec are filled as the gen one, this will 
            // have to be accounted for in further corrections since we are not unfolding them)
            if (true_pass_cuts_eec) {
                fill_jk_resampling(histos_efficiency_denominator_eecpt, num, dr_reco_j, jpt_gen, weight*eec_reco_j);
                if (num<0.5) h_half0_efficiency_denominator_eecpt->Fill(dr_reco_j, jpt_gen, weight*eec_reco_j);
                else h_half1_efficiency_denominator_eecpt->Fill(dr_reco_j, jpt_gen, weight*eec_reco_j);
            }
            if (reco_pass_cuts_eec) {
                fill_jk_resampling(histos_purity_denominator_eecpt, num, dr_reco_j, jpt_reco, weight*eec_reco_j);
                if (num<0.5) h_half0_purity_denominator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
                else h_half1_purity_denominator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
            }
            if (true_pass_cuts_eec && reco_pass_cuts_eec) {
                fill_jk_resampling(histos_efficiency_numerator_eecpt, num, dr_reco_j, jpt_gen, weight*eec_reco_j);
                fill_jk_resampling(histos_purity_numerator_eecpt, num, dr_reco_j, jpt_reco, weight*eec_reco_j);
                fill_jk_resampling_response(responses_eecpt, num, dr_reco_j, jpt_reco, dr_reco_j, jpt_gen, weight*eec_reco_j);

                if (num<0.5) {
                    h_half0_efficiency_numerator_eecpt->Fill(dr_reco_j, jpt_gen, weight*eec_reco_j);
                    h_half0_purity_numerator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
                    response_half0_eecpt->Fill(dr_reco_j, jpt_reco, dr_reco_j, jpt_gen, weight*eec_reco_j);
                } else {
                    h_half1_efficiency_numerator_eecpt->Fill(dr_reco_j, jpt_gen, weight*eec_reco_j);
                    h_half1_purity_numerator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
                    response_half1_eecpt->Fill(dr_reco_j, jpt_reco, dr_reco_j, jpt_gen, weight*eec_reco_j);
                }
                response_full_eecpt->Fill(dr_reco_j, jpt_reco, dr_reco_j, jpt_gen, weight*eec_reco_j);
            }
        }// pair entry loop
    } // tree entry loop

    // Create purity and efficiency histograms
    TH2D *h0_purity_eecpt, 
        *h1_purity_eecpt,
        *h2_purity_eecpt,
        *h3_purity_eecpt,
        *h4_purity_eecpt,
        *h5_purity_eecpt,
        *h6_purity_eecpt,
        *h7_purity_eecpt,
        *h8_purity_eecpt,
        *h9_purity_eecpt;

    std::vector<TH2D *> histos_purity_eecpt = {
        h0_purity_eecpt,
        h1_purity_eecpt,
        h2_purity_eecpt,
        h3_purity_eecpt,
        h4_purity_eecpt,
        h5_purity_eecpt,
        h6_purity_eecpt,
        h7_purity_eecpt,
        h8_purity_eecpt,
        h9_purity_eecpt,
    };

    TH2D *h0_efficiency_eecpt, 
        *h1_efficiency_eecpt,
        *h2_efficiency_eecpt,
        *h3_efficiency_eecpt,
        *h4_efficiency_eecpt,
        *h5_efficiency_eecpt,
        *h6_efficiency_eecpt,
        *h7_efficiency_eecpt,
        *h8_efficiency_eecpt,
        *h9_efficiency_eecpt;

    std::vector<TH2D *> histos_efficiency_eecpt = {
        h0_efficiency_eecpt,
        h1_efficiency_eecpt,
        h2_efficiency_eecpt,
        h3_efficiency_eecpt,
        h4_efficiency_eecpt,
        h5_efficiency_eecpt,
        h6_efficiency_eecpt,
        h7_efficiency_eecpt,
        h8_efficiency_eecpt,
        h9_efficiency_eecpt,
    };

    
    
    // initialize histograms
    for (int i=0; i<10; i++) {
        histos_purity_eecpt[i] = (TH2D *) histos_purity_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_eecpt");
        histos_purity_eecpt[i]->Divide(histos_purity_numerator_eecpt[i], histos_purity_denominator_eecpt[i], 1., 1., "b");
        histos_efficiency_eecpt[i] = (TH2D *) histos_efficiency_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_eecpt");
        histos_efficiency_eecpt[i]->Divide(histos_efficiency_numerator_eecpt[i], histos_efficiency_denominator_eecpt[i], 1., 1., "b");
    }

    // declare the per half purity + efficiency histograms
    TH2D *h_half0_purity_eecpt = (TH2D *) h_half0_purity_numerator_eecpt->Clone("h_half0_purity_eecpt");
    h_half0_purity_eecpt->Divide(h_half0_purity_numerator_eecpt, h_half0_purity_denominator_eecpt, 1., 1., "b");
    TH2D *h_half0_efficiency_eecpt = (TH2D *) h_half0_efficiency_numerator_eecpt->Clone("h_half0_efficiency_eecpt");
    h_half0_efficiency_eecpt->Divide(h_half0_efficiency_numerator_eecpt, h_half0_efficiency_denominator_eecpt, 1., 1., "b");

    TH2D *h_half1_purity_eecpt = (TH2D *) h_half1_purity_numerator_eecpt->Clone("h_half1_purity_eecpt");
    h_half1_purity_eecpt->Divide(h_half1_purity_numerator_eecpt, h_half1_purity_denominator_eecpt, 1., 1., "b");
    TH2D *h_half1_efficiency_eecpt = (TH2D *) h_half1_efficiency_numerator_eecpt->Clone("h_half1_efficiency_eecpt");
    h_half1_efficiency_eecpt->Divide(h_half1_efficiency_numerator_eecpt, h_half1_efficiency_denominator_eecpt, 1., 1., "b");

    // declare the full purity + efficiency histograms
    TH2D *h_full_purity_numerator_eecpt = (TH2D *) h_half0_purity_numerator_eecpt->Clone("h_full_purity_numerator_eecpt");
    h_full_purity_numerator_eecpt->Add(h_half1_purity_numerator_eecpt);
    TH2D *h_full_purity_denominator_eecpt = (TH2D *) h_half0_purity_denominator_eecpt->Clone("h_full_purity_denominator_eecpt");
    h_full_purity_denominator_eecpt->Add(h_half1_purity_denominator_eecpt);
    TH2D *h_full_purity_eecpt = (TH2D *) h_full_purity_numerator_eecpt->Clone("h_full_purity_eecpt");
    h_full_purity_eecpt->Divide(h_full_purity_numerator_eecpt, h_full_purity_denominator_eecpt, 1., 1., "b");

    TH2D *h_full_efficiency_numerator_eecpt = (TH2D *) h_half0_efficiency_numerator_eecpt->Clone("h_full_efficiency_numerator_eecpt");
    h_full_efficiency_numerator_eecpt->Add(h_half1_efficiency_numerator_eecpt);
    TH2D *h_full_efficiency_denominator_eecpt = (TH2D *) h_half0_efficiency_denominator_eecpt->Clone("h_full_efficiency_denominator_eecpt");
    h_full_efficiency_denominator_eecpt->Add(h_half1_efficiency_denominator_eecpt);
    TH2D *h_full_efficiency_eecpt = (TH2D *) h_full_efficiency_numerator_eecpt->Clone("h_full_efficiency_eecpt");
    h_full_efficiency_eecpt->Divide(h_full_efficiency_numerator_eecpt, h_full_efficiency_denominator_eecpt, 1., 1., "b");

    // Create output file
    std::cout << "Creating file: " << folder+fout_name << std::endl;
    TFile *fout = new TFile(folder+fout_name, "recreate");

    // Write jk resampling histograms + responses (Lida)
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i]->Write();
        histos_purity_denominator_eecpt[i]->Write();
        histos_purity_eecpt[i]->Write();

        histos_efficiency_numerator_eecpt[i]->Write();
        histos_efficiency_denominator_eecpt[i]->Write();
        histos_efficiency_eecpt[i]->Write();

        responses_eecpt[i]->Write();

    }

    // Write per half histograms 
    h_half0_purity_numerator_eecpt->Write();
    h_half0_purity_denominator_eecpt->Write();
    h_half0_purity_eecpt->Write();

    h_half0_efficiency_numerator_eecpt->Write();
    h_half0_efficiency_denominator_eecpt->Write();
    h_half0_efficiency_eecpt->Write();

    response_half0_eecpt->Write();

    h_half1_purity_numerator_eecpt->Write();
    h_half1_purity_denominator_eecpt->Write();
    h_half1_purity_eecpt->Write();

    h_half1_efficiency_numerator_eecpt->Write();
    h_half1_efficiency_denominator_eecpt->Write();
    h_half1_efficiency_eecpt->Write();

    response_half1_eecpt->Write();
    
    // Write full histograms 
    h_full_purity_numerator_eecpt->Write();
    h_full_purity_denominator_eecpt->Write();
    h_full_purity_eecpt->Write();

    h_full_efficiency_numerator_eecpt->Write();
    h_full_efficiency_denominator_eecpt->Write();
    h_full_efficiency_eecpt->Write();

    response_full_eecpt->Write();

    fout->Close();
    delete fout;
}

//Creates a 2D response matrix and purity/efficiency corrections from a tree
void create_response_2D(TString &filename,  TString &dataset, TString &label, TString &folder, bool &btag, Int_t n, Float_t &pT_low, Float_t &pT_high)
{
    TString flav = label;    
    std::cout << "flav:" << flav << std::endl;
   
    //Create the fout name depending on the selection
    TString fout_name = "histos_response_2D_";

    if(!btag) label += "_notag"; 

    fout_name += TString(Form("n%i_", n)) + dataset + "_" + label + ".root";

    std::cout << "fin: " << folder+filename << std::endl;
    TFile *fin = new TFile(folder+filename);

    TString tree_name = "tree";
    
    std::cout << "tree: " << tree_name << std::endl;
    TTree *tree = (TTree *) fin->Get(tree_name);

    // Set tree addresses
    Float_t weight, pthat, jpt_reco, jpt_gen, jt_eta_reco, jt_eta_gen, discr;
    Int_t ndr_reco, ndr_gen, ndr_reco_tot, ndr_gen_tot, jtHadFlav, jtNbHad;
    Float_t dr_reco[4000], dr_gen[4000], eec_reco[4000], eec_gen[4000];


    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("pthat", &pthat);
    tree->SetBranchAddress("ndr_reco", &ndr_reco);
    tree->SetBranchAddress("ndr_gen", &ndr_gen);
    tree->SetBranchAddress("ndr_reco_tot", &ndr_reco_tot);
    tree->SetBranchAddress("ndr_gen_tot", &ndr_gen_tot);
    tree->SetBranchAddress("jpt_reco", &jpt_reco);
    tree->SetBranchAddress("jpt_gen", &jpt_gen);
    tree->SetBranchAddress("dr_reco", &dr_reco);
    tree->SetBranchAddress("dr_gen", &dr_gen);
    tree->SetBranchAddress("eec_reco", &eec_reco);
    tree->SetBranchAddress("eec_gen", &eec_gen);
    tree->SetBranchAddress("jt_eta_reco", &jt_eta_reco);
    tree->SetBranchAddress("jt_eta_gen", &jt_eta_gen);
    tree->SetBranchAddress("discr", &discr);
    tree->SetBranchAddress("jtHadFlav", &jtHadFlav);                                                               
    tree->SetBranchAddress("jtNbHad", &jtNbHad);                                                                                                                                           


    // random number generator for jackknife resampling (Lida)
    const double range_from = 0;
    const double range_to = 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(range_from, range_to);

    
    // Declare histograms
    TH2D *h_half0_purity_numerator_eecpt = new TH2D("h_half0_purity_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half0_purity_denominator_eecpt = new TH2D("h_half0_purity_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half0_efficiency_numerator_eecpt = new TH2D("h_half0_efficiency_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half0_efficiency_denominator_eecpt = new TH2D("h_half0_efficiency_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_purity_numerator_eecpt = new TH2D("h_half1_purity_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_purity_denominator_eecpt = new TH2D("h_half1_purity_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_efficiency_numerator_eecpt = new TH2D("h_half1_efficiency_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_half1_efficiency_denominator_eecpt = new TH2D("h_half1_efficiency_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);

    RooUnfoldResponse *response_half0_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half0", "response for 2d: eec and jet pt"); 
    RooUnfoldResponse *response_half1_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half1", "response for 2d: eec and jet pt"); 
    RooUnfoldResponse *response_full_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_full_eecpt", "response for 2d: eec and jet pt"); 

    // rg jk resampling histogram definiton (Lida)
    TH2D *h0_purity_numerator_eecpt, 
        *h1_purity_numerator_eecpt, 
        *h2_purity_numerator_eecpt, 
        *h3_purity_numerator_eecpt, 
        *h4_purity_numerator_eecpt, 
        *h5_purity_numerator_eecpt, 
        *h6_purity_numerator_eecpt, 
        *h7_purity_numerator_eecpt, 
        *h8_purity_numerator_eecpt, 
        *h9_purity_numerator_eecpt;

    std::vector<TH2D *> histos_purity_numerator_eecpt = {
        h0_purity_numerator_eecpt, 
        h1_purity_numerator_eecpt, 
        h2_purity_numerator_eecpt, 
        h3_purity_numerator_eecpt, 
        h4_purity_numerator_eecpt, 
        h5_purity_numerator_eecpt, 
        h6_purity_numerator_eecpt, 
        h7_purity_numerator_eecpt, 
        h8_purity_numerator_eecpt, 
        h9_purity_numerator_eecpt,
    };

    TH2D *h0_purity_denominator_eecpt, 
        *h1_purity_denominator_eecpt, 
        *h2_purity_denominator_eecpt, 
        *h3_purity_denominator_eecpt, 
        *h4_purity_denominator_eecpt, 
        *h5_purity_denominator_eecpt, 
        *h6_purity_denominator_eecpt, 
        *h7_purity_denominator_eecpt, 
        *h8_purity_denominator_eecpt, 
        *h9_purity_denominator_eecpt;

    std::vector<TH2D *> histos_purity_denominator_eecpt = {
        h0_purity_denominator_eecpt, 
        h1_purity_denominator_eecpt, 
        h2_purity_denominator_eecpt, 
        h3_purity_denominator_eecpt, 
        h4_purity_denominator_eecpt, 
        h5_purity_denominator_eecpt, 
        h6_purity_denominator_eecpt, 
        h7_purity_denominator_eecpt, 
        h8_purity_denominator_eecpt, 
        h9_purity_denominator_eecpt,
    };

    TH2D *h0_efficiency_numerator_eecpt, 
        *h1_efficiency_numerator_eecpt, 
        *h2_efficiency_numerator_eecpt, 
        *h3_efficiency_numerator_eecpt, 
        *h4_efficiency_numerator_eecpt, 
        *h5_efficiency_numerator_eecpt, 
        *h6_efficiency_numerator_eecpt, 
        *h7_efficiency_numerator_eecpt, 
        *h8_efficiency_numerator_eecpt, 
        *h9_efficiency_numerator_eecpt;

    std::vector<TH2D *> histos_efficiency_numerator_eecpt = {
        h0_efficiency_numerator_eecpt, 
        h1_efficiency_numerator_eecpt, 
        h2_efficiency_numerator_eecpt, 
        h3_efficiency_numerator_eecpt, 
        h4_efficiency_numerator_eecpt, 
        h5_efficiency_numerator_eecpt, 
        h6_efficiency_numerator_eecpt, 
        h7_efficiency_numerator_eecpt, 
        h8_efficiency_numerator_eecpt, 
        h9_efficiency_numerator_eecpt,
    };

    TH2D *h0_efficiency_denominator_eecpt, 
        *h1_efficiency_denominator_eecpt, 
        *h2_efficiency_denominator_eecpt, 
        *h3_efficiency_denominator_eecpt, 
        *h4_efficiency_denominator_eecpt, 
        *h5_efficiency_denominator_eecpt, 
        *h6_efficiency_denominator_eecpt, 
        *h7_efficiency_denominator_eecpt, 
        *h8_efficiency_denominator_eecpt, 
        *h9_efficiency_denominator_eecpt;

    std::vector<TH2D *> histos_efficiency_denominator_eecpt = {
        h0_efficiency_denominator_eecpt, 
        h1_efficiency_denominator_eecpt, 
        h2_efficiency_denominator_eecpt, 
        h3_efficiency_denominator_eecpt, 
        h4_efficiency_denominator_eecpt, 
        h5_efficiency_denominator_eecpt, 
        h6_efficiency_denominator_eecpt, 
        h7_efficiency_denominator_eecpt, 
        h8_efficiency_denominator_eecpt, 
        h9_efficiency_denominator_eecpt,
    };

    RooUnfoldResponse *response0_eecpt, 
        *response1_eecpt,
        *response2_eecpt,
        *response3_eecpt,
        *response4_eecpt,
        *response5_eecpt,
        *response6_eecpt,
        *response7_eecpt,
        *response8_eecpt,
        *response9_eecpt;

    std::vector<RooUnfoldResponse *> responses_eecpt = {
        response0_eecpt, 
        response1_eecpt,
        response2_eecpt,
        response3_eecpt,
        response4_eecpt,
        response5_eecpt,
        response6_eecpt,
        response7_eecpt,
        response8_eecpt,
        response9_eecpt,
    };
    

    // initialize histograms and responses
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_purity_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
        histos_purity_denominator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_purity_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_numerator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_efficiency_numerator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_denominator_eecpt[i] = new TH2D("h"+TString(Form("%d",i))+"_efficiency_denominator_eecpt", "x=dr, y=jtpt", dr_bins, dr_binsVector, jtpt_bins, jtpt_binsVector);

        responses_eecpt[i] = new RooUnfoldResponse(histos_purity_denominator_eecpt[i],histos_efficiency_denominator_eecpt[i], "response"+TString(Form("%d",i))+"_eecpt", "response for 2d: eec and jet pt"); 
    }        

    // Loop over tree entries
    Long64_t nentries = tree->GetEntries();
    std::cout << "Entries: " << nentries << std::endl;
    for (Long64_t ient = 0; ient < nentries; ient++) {
        //Print progress
        if (ient%1000000==0) cout << "ient=" << ient << std::endl; 
    	tree->GetEntry(ient);
	Int_t cuts = 1;                                                                                                                                               
        if      (flav == "b1")    cuts = 1;                                                                                                                           
        else if (flav == "b2")    cuts = 2;                                                                                                                           
        else if (flav == "nonb")  cuts = 3;                                                                                                                           
        else if (flav == "all")   cuts = 4;                                                                                                                           
        else if (flav == "c")     cuts = 5;                                                                                                                           
        else if (flav == "light") cuts = 6;                                                                                                                           
        else {                                                                                                                                                      
	  std::cerr << "ERROR: Unknown label '" << label << "'\n";                                                                                                  
	  exit(1);                                                                                                                                                       
	}       

	bool skip = false;                                                                  
    // Select jet flavour and/or select on the number of b hadrons                                                                                                                                                                                                                                          
        switch(cuts){                                                                                                                                                
	  //b-jet with one b hadron                                                                                                                                  
	case 1:                                                                                                                                                               	  
    if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                         	  
    if (std::abs(jtNbHad) != 1) skip = true;                                                                                                                        	  
    break;                                                                                                                                                       	  
    //b-jet with more than 1 b hadron                                                                                                                           
	case 2:                                                                                                                                                              	  
    if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                          	  
    if (std::abs(jtNbHad) < 2) skip = true;                                                                                                                     	  
    break;                                                                                                                                                                   	  
    //non-b jets                                                                                                                                                                                                                                                                                          
	case 3:                                                                                                                                                              	  
    if (std::abs(jtHadFlav) == 5) skip = true;                                                                                                                         	  
    break;                                                                                                                                                             	  
    //no flavour selection                                                                                                                                                                                                                                                                                
	case 4:                                                                                                                                                              	  
    skip = false;                                                                                                                                                    	  
    break;                                                                                                                                                         	  
    //c-jets                                                                                                                                                     
	case 5:                                                                                                                                                              	  
    if(std::abs(jtHadFlav) != 4) skip = true;                                                                                                                 	  
    break;                                                                                                                                                          	  
    //light (non-b non-c) jets                                                                                                                                                                                                                                                                            
	case 6:                                                                                                                                                       	  
    if(std::abs(jtHadFlav) >= 4) skip = true;                                                                                                                         	  
    break;                                                                                                                                                                  }                                                                                                                                                                     if (skip) continue;      
	/////////////////////////////

        if (skipMC(jpt_reco, jpt_gen, pthat)) continue;

        // Random number for jack-knife resampling (Lida)
        double num = distr(generator);

        // --------------------------

        // Check if pass cuts
        bool has_gen_match = (jpt_gen > 0);

        // Fill histograms
        if (!has_gen_match) {   
            // fill fakes
            continue; 
        } 

        
        // Skip jets outside tracker 
        if (std::abs(jt_eta_reco) > 1.6) continue;
        if (std::abs(jt_eta_gen) > 1.6) continue;

        //Select jets passing the reco b-jet tagging 
        if (btag && std::abs(discr) <= 0.99) continue;


        
        // The rest of the histograms don;t include any fakes

        
        //Loop over matched pairs to fill the response matrix 
        for (Int_t j = 0; j < ndr_reco; j++){

            Float_t dr_reco_j = dr_reco[j];
            Float_t dr_gen_j = dr_gen[j];

            Float_t eec_reco_j = eec_reco[j];



            num = distr(generator);

            
            //checks for underflow/overflow
            if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
            if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;

            
            if(dr_gen_j >= dr_max) dr_gen_j = dr_max_fill;
            if(dr_gen_j < dr_min) dr_gen_j = dr_min_fill;

            bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min && dr_reco_j < dr_max && dr_reco_j >= dr_min);
            bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min && dr_gen_j < dr_max && dr_gen_j >= dr_min);

            if (true_pass_cuts_eec && reco_pass_cuts_eec) {
                fill_jk_resampling(histos_efficiency_numerator_eecpt, num, dr_gen_j, jpt_gen, weight*eec_reco_j);
                fill_jk_resampling(histos_purity_numerator_eecpt, num, dr_reco_j, jpt_reco, weight*eec_reco_j);
                fill_jk_resampling_response(responses_eecpt, num, dr_reco_j, jpt_reco, dr_gen_j, jpt_gen, weight*eec_reco_j);

                if (num<0.5) {
                    h_half0_efficiency_numerator_eecpt->Fill(dr_gen_j, jpt_gen, weight*eec_reco_j);
                    h_half0_purity_numerator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
                    response_half0_eecpt->Fill(dr_reco_j, jpt_reco, dr_gen_j, jpt_gen, weight*eec_reco_j);
                } else {
                    h_half1_efficiency_numerator_eecpt->Fill(dr_gen_j, jpt_gen, weight*eec_reco_j);
                    h_half1_purity_numerator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
                    response_half1_eecpt->Fill(dr_reco_j, jpt_reco, dr_gen_j, jpt_gen, weight*eec_reco_j);
                }
                response_full_eecpt->Fill(dr_reco_j, jpt_reco, dr_gen_j, jpt_gen, weight*eec_reco_j);
            }
        }// pair entry loop

        //Loop over all (matched and non-matched) dr to fill the eec purity and efficiency
        //For purity (reco)
        for (Int_t j = 0; j < ndr_reco_tot; j++){

            Float_t dr_reco_j = dr_reco[j];

            Float_t eec_reco_j = eec_reco[j];



            num = distr(generator);

            
            //checks for underflow/overflow
            if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
            if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;


            bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min && dr_reco_j < dr_max && dr_reco_j >= dr_min);


            // fill purity
            if (reco_pass_cuts_eec) {
                fill_jk_resampling(histos_purity_denominator_eecpt, num, dr_reco_j, jpt_reco, weight*eec_reco_j);
                if (num<0.5) h_half0_purity_denominator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
                else h_half1_purity_denominator_eecpt->Fill(dr_reco_j, jpt_reco, weight*eec_reco_j);
            }
        }// pair entry loop

        //For efficiency (gen)
        for (Int_t j = 0; j < ndr_gen_tot; j++){

            Float_t dr_gen_j = dr_gen[j];

            Float_t eec_gen_j = eec_gen[j];



            num = distr(generator);

            
            //checks for underflow/overflow            
            if(dr_gen_j >= dr_max) dr_gen_j = dr_max_fill;
            if(dr_gen_j < dr_min) dr_gen_j = dr_min_fill;

            bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min && dr_gen_j < dr_max && dr_gen_j >= dr_min);



            // fill eec histograms
            if (true_pass_cuts_eec) {
                fill_jk_resampling(histos_efficiency_denominator_eecpt, num, dr_gen_j, jpt_gen, weight*eec_gen_j);
                if (num<0.5) h_half0_efficiency_denominator_eecpt->Fill(dr_gen_j, jpt_gen, weight*eec_gen_j);
                else h_half1_efficiency_denominator_eecpt->Fill(dr_gen_j, jpt_gen, weight*eec_gen_j);
            }
        }// pair entry loop
    } // tree entry loop

    // Create purity and efficiency histograms
    TH2D *h0_purity_eecpt, 
        *h1_purity_eecpt,
        *h2_purity_eecpt,
        *h3_purity_eecpt,
        *h4_purity_eecpt,
        *h5_purity_eecpt,
        *h6_purity_eecpt,
        *h7_purity_eecpt,
        *h8_purity_eecpt,
        *h9_purity_eecpt;

    std::vector<TH2D *> histos_purity_eecpt = {
        h0_purity_eecpt,
        h1_purity_eecpt,
        h2_purity_eecpt,
        h3_purity_eecpt,
        h4_purity_eecpt,
        h5_purity_eecpt,
        h6_purity_eecpt,
        h7_purity_eecpt,
        h8_purity_eecpt,
        h9_purity_eecpt,
    };

    TH2D *h0_efficiency_eecpt, 
        *h1_efficiency_eecpt,
        *h2_efficiency_eecpt,
        *h3_efficiency_eecpt,
        *h4_efficiency_eecpt,
        *h5_efficiency_eecpt,
        *h6_efficiency_eecpt,
        *h7_efficiency_eecpt,
        *h8_efficiency_eecpt,
        *h9_efficiency_eecpt;

    std::vector<TH2D *> histos_efficiency_eecpt = {
        h0_efficiency_eecpt,
        h1_efficiency_eecpt,
        h2_efficiency_eecpt,
        h3_efficiency_eecpt,
        h4_efficiency_eecpt,
        h5_efficiency_eecpt,
        h6_efficiency_eecpt,
        h7_efficiency_eecpt,
        h8_efficiency_eecpt,
        h9_efficiency_eecpt,
    };

    
    
    // initialize histograms
    for (int i=0; i<10; i++) {
        histos_purity_eecpt[i] = (TH2D *) histos_purity_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_eecpt");
        histos_purity_eecpt[i]->Divide(histos_purity_numerator_eecpt[i], histos_purity_denominator_eecpt[i], 1., 1., "b");
        histos_efficiency_eecpt[i] = (TH2D *) histos_efficiency_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_eecpt");
        histos_efficiency_eecpt[i]->Divide(histos_efficiency_numerator_eecpt[i], histos_efficiency_denominator_eecpt[i], 1., 1., "b");
    }

    // declare the per half purity + efficiency histograms
    TH2D *h_half0_purity_eecpt = (TH2D *) h_half0_purity_numerator_eecpt->Clone("h_half0_purity_eecpt");
    h_half0_purity_eecpt->Divide(h_half0_purity_numerator_eecpt, h_half0_purity_denominator_eecpt, 1., 1., "b");
    TH2D *h_half0_efficiency_eecpt = (TH2D *) h_half0_efficiency_numerator_eecpt->Clone("h_half0_efficiency_eecpt");
    h_half0_efficiency_eecpt->Divide(h_half0_efficiency_numerator_eecpt, h_half0_efficiency_denominator_eecpt, 1., 1., "b");

    TH2D *h_half1_purity_eecpt = (TH2D *) h_half1_purity_numerator_eecpt->Clone("h_half1_purity_eecpt");
    h_half1_purity_eecpt->Divide(h_half1_purity_numerator_eecpt, h_half1_purity_denominator_eecpt, 1., 1., "b");
    TH2D *h_half1_efficiency_eecpt = (TH2D *) h_half1_efficiency_numerator_eecpt->Clone("h_half1_efficiency_eecpt");
    h_half1_efficiency_eecpt->Divide(h_half1_efficiency_numerator_eecpt, h_half1_efficiency_denominator_eecpt, 1., 1., "b");

    // declare the full purity + efficiency histograms
    TH2D *h_full_purity_numerator_eecpt = (TH2D *) h_half0_purity_numerator_eecpt->Clone("h_full_purity_numerator_eecpt");
    h_full_purity_numerator_eecpt->Add(h_half1_purity_numerator_eecpt);
    TH2D *h_full_purity_denominator_eecpt = (TH2D *) h_half0_purity_denominator_eecpt->Clone("h_full_purity_denominator_eecpt");
    h_full_purity_denominator_eecpt->Add(h_half1_purity_denominator_eecpt);
    TH2D *h_full_purity_eecpt = (TH2D *) h_full_purity_numerator_eecpt->Clone("h_full_purity_eecpt");
    h_full_purity_eecpt->Divide(h_full_purity_numerator_eecpt, h_full_purity_denominator_eecpt, 1., 1., "b");

    TH2D *h_full_efficiency_numerator_eecpt = (TH2D *) h_half0_efficiency_numerator_eecpt->Clone("h_full_efficiency_numerator_eecpt");
    h_full_efficiency_numerator_eecpt->Add(h_half1_efficiency_numerator_eecpt);
    TH2D *h_full_efficiency_denominator_eecpt = (TH2D *) h_half0_efficiency_denominator_eecpt->Clone("h_full_efficiency_denominator_eecpt");
    h_full_efficiency_denominator_eecpt->Add(h_half1_efficiency_denominator_eecpt);
    TH2D *h_full_efficiency_eecpt = (TH2D *) h_full_efficiency_numerator_eecpt->Clone("h_full_efficiency_eecpt");
    h_full_efficiency_eecpt->Divide(h_full_efficiency_numerator_eecpt, h_full_efficiency_denominator_eecpt, 1., 1., "b");

    // Create output file
    std::cout << "Creating file: " << folder+fout_name << std::endl;
    TFile *fout = new TFile(folder+fout_name, "recreate");

    // Write jk resampling histograms + responses (Lida)
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i]->Write();
        histos_purity_denominator_eecpt[i]->Write();
        histos_purity_eecpt[i]->Write();

        histos_efficiency_numerator_eecpt[i]->Write();
        histos_efficiency_denominator_eecpt[i]->Write();
        histos_efficiency_eecpt[i]->Write();

        responses_eecpt[i]->Write();

    }

    // Write per half histograms 
    h_half0_purity_numerator_eecpt->Write();
    h_half0_purity_denominator_eecpt->Write();
    h_half0_purity_eecpt->Write();

    h_half0_efficiency_numerator_eecpt->Write();
    h_half0_efficiency_denominator_eecpt->Write();
    h_half0_efficiency_eecpt->Write();

    response_half0_eecpt->Write();

    h_half1_purity_numerator_eecpt->Write();
    h_half1_purity_denominator_eecpt->Write();
    h_half1_purity_eecpt->Write();

    h_half1_efficiency_numerator_eecpt->Write();
    h_half1_efficiency_denominator_eecpt->Write();
    h_half1_efficiency_eecpt->Write();

    response_half1_eecpt->Write();
    
    // Write full histograms 
    h_full_purity_numerator_eecpt->Write();
    h_full_purity_denominator_eecpt->Write();
    h_full_purity_eecpt->Write();

    h_full_efficiency_numerator_eecpt->Write();
    h_full_efficiency_denominator_eecpt->Write();
    h_full_efficiency_eecpt->Write();

    response_full_eecpt->Write();

    fout->Close();
    delete fout;
}

//Creates a 3D response matrix and purity/efficiency corrections from a tree
void create_response_3D(TString &filename,  TString &sample, TString &label, TString &folder, bool btag, Int_t &n, Float_t &pT_low, Float_t &pT_high)
{    //Create the fout name depending on the selection

    TString flav = label; 
    std::cout << "flav:" << flav << std::endl;
    TString fout_name = "histos_response_3D_";

    if(!btag) label += "_notag"; 

    fout_name += TString(Form("n%i_", n)) + sample + "_" + label + ".root";


    std::cout << "fin: " << filename << std::endl;
    TFile *fin = new TFile(folder + filename);

    TString tree_name = "tree";
    std::cout << "tree: " << tree_name << std::endl;
    TTree *tree = (TTree *) fin->Get(tree_name);

    // Set tree addresses


    Float_t weight, pthat, jpt_reco, jpt_gen, jt_eta_reco, jt_eta_gen, discr;                                                                                         
    Int_t ndr_reco, ndr_gen, ndr_reco_tot, ndr_gen_tot, jtHadFlav, jtNbHad;                                                                                           
    Float_t dr_reco[4000], dr_gen[4000], eec_reco[4000], eec_gen[4000];                                                                                                
        

    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("pthat", &pthat);
    tree->SetBranchAddress("ndr_reco", &ndr_reco);
    tree->SetBranchAddress("ndr_gen", &ndr_gen);
    tree->SetBranchAddress("ndr_reco_tot", &ndr_reco_tot);
    tree->SetBranchAddress("ndr_gen_tot", &ndr_gen_tot);
    tree->SetBranchAddress("jpt_reco", &jpt_reco);
    tree->SetBranchAddress("jpt_gen", &jpt_gen);
    tree->SetBranchAddress("dr_reco", &dr_reco);
    tree->SetBranchAddress("dr_gen", &dr_gen);
    tree->SetBranchAddress("eec_reco", &eec_reco);
    tree->SetBranchAddress("eec_gen", &eec_gen);
    tree->SetBranchAddress("jt_eta_reco", &jt_eta_reco);
    tree->SetBranchAddress("jt_eta_gen", &jt_eta_gen);
    tree->SetBranchAddress("discr", &discr);
    tree->SetBranchAddress("jtHadFlav", &jtHadFlav);                                                                                         
    tree->SetBranchAddress("jtNbHad", &jtNbHad);                                                                                             
        

    Int_t dr_bins = bins_dr;

    // random number generator for jackknife resampling (Lida)
    const double range_from = 0;
    const double range_to = 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(range_from, range_to);

    //Check percentage of non-matched tracks
    Double_t matched_reco = 0;
    Double_t matched_gen = 0;
    Double_t notmatched_reco = 0;
    Double_t notmatched_gen = 0;

    
    // Declare histograms
    TH3D *h_half0_purity_numerator_eecpt = new TH3D("h_half0_purity_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half0_purity_denominator_eecpt = new TH3D("h_half0_purity_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half0_efficiency_numerator_eecpt = new TH3D("h_half0_efficiency_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half0_efficiency_denominator_eecpt = new TH3D("h_half0_efficiency_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_purity_numerator_eecpt = new TH3D("h_half1_purity_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_purity_denominator_eecpt = new TH3D("h_half1_purity_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_efficiency_numerator_eecpt = new TH3D("h_half1_efficiency_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_efficiency_denominator_eecpt = new TH3D("h_half1_efficiency_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    // histograms only for the non matched pairs
    TH3D *h_notmatched_reco = new TH3D("h_notmatched_reco", "h_notmatched_reco", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_notmatched_gen = new TH3D("h_notmatched_gen", "h_notmatched_gen", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);

    RooUnfoldResponse *response_half0_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half0", "response for 3d: eec and jet pt"); 
    RooUnfoldResponse *response_half1_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half1", "response for 3d: eec and jet pt"); 
    RooUnfoldResponse *response_full_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_full_eecpt", "response for 3d: eec and jet pt"); 


    // rg jk resampling histogram definiton (Lida)
    TH3D *h0_purity_numerator_eecpt, 
        *h1_purity_numerator_eecpt, 
        *h2_purity_numerator_eecpt, 
        *h3_purity_numerator_eecpt, 
        *h4_purity_numerator_eecpt, 
        *h5_purity_numerator_eecpt, 
        *h6_purity_numerator_eecpt, 
        *h7_purity_numerator_eecpt, 
        *h8_purity_numerator_eecpt, 
        *h9_purity_numerator_eecpt;

    std::vector<TH3D *> histos_purity_numerator_eecpt = {
        h0_purity_numerator_eecpt, 
        h1_purity_numerator_eecpt, 
        h2_purity_numerator_eecpt, 
        h3_purity_numerator_eecpt, 
        h4_purity_numerator_eecpt, 
        h5_purity_numerator_eecpt, 
        h6_purity_numerator_eecpt, 
        h7_purity_numerator_eecpt, 
        h8_purity_numerator_eecpt, 
        h9_purity_numerator_eecpt,
    };

    TH3D *h0_purity_denominator_eecpt, 
        *h1_purity_denominator_eecpt, 
        *h2_purity_denominator_eecpt, 
        *h3_purity_denominator_eecpt, 
        *h4_purity_denominator_eecpt, 
        *h5_purity_denominator_eecpt, 
        *h6_purity_denominator_eecpt, 
        *h7_purity_denominator_eecpt, 
        *h8_purity_denominator_eecpt, 
        *h9_purity_denominator_eecpt;

    std::vector<TH3D *> histos_purity_denominator_eecpt = {
        h0_purity_denominator_eecpt, 
        h1_purity_denominator_eecpt, 
        h2_purity_denominator_eecpt, 
        h3_purity_denominator_eecpt, 
        h4_purity_denominator_eecpt, 
        h5_purity_denominator_eecpt, 
        h6_purity_denominator_eecpt, 
        h7_purity_denominator_eecpt, 
        h8_purity_denominator_eecpt, 
        h9_purity_denominator_eecpt,
    };

    TH3D *h0_efficiency_numerator_eecpt, 
        *h1_efficiency_numerator_eecpt, 
        *h2_efficiency_numerator_eecpt, 
        *h3_efficiency_numerator_eecpt, 
        *h4_efficiency_numerator_eecpt, 
        *h5_efficiency_numerator_eecpt, 
        *h6_efficiency_numerator_eecpt, 
        *h7_efficiency_numerator_eecpt, 
        *h8_efficiency_numerator_eecpt, 
        *h9_efficiency_numerator_eecpt;

    std::vector<TH3D *> histos_efficiency_numerator_eecpt = {
        h0_efficiency_numerator_eecpt, 
        h1_efficiency_numerator_eecpt, 
        h2_efficiency_numerator_eecpt, 
        h3_efficiency_numerator_eecpt, 
        h4_efficiency_numerator_eecpt, 
        h5_efficiency_numerator_eecpt, 
        h6_efficiency_numerator_eecpt, 
        h7_efficiency_numerator_eecpt, 
        h8_efficiency_numerator_eecpt, 
        h9_efficiency_numerator_eecpt,
    };

    TH3D *h0_efficiency_denominator_eecpt, 
        *h1_efficiency_denominator_eecpt, 
        *h2_efficiency_denominator_eecpt, 
        *h3_efficiency_denominator_eecpt, 
        *h4_efficiency_denominator_eecpt, 
        *h5_efficiency_denominator_eecpt, 
        *h6_efficiency_denominator_eecpt, 
        *h7_efficiency_denominator_eecpt, 
        *h8_efficiency_denominator_eecpt, 
        *h9_efficiency_denominator_eecpt;

    std::vector<TH3D *> histos_efficiency_denominator_eecpt = {
        h0_efficiency_denominator_eecpt, 
        h1_efficiency_denominator_eecpt, 
        h2_efficiency_denominator_eecpt, 
        h3_efficiency_denominator_eecpt, 
        h4_efficiency_denominator_eecpt, 
        h5_efficiency_denominator_eecpt, 
        h6_efficiency_denominator_eecpt, 
        h7_efficiency_denominator_eecpt, 
        h8_efficiency_denominator_eecpt, 
        h9_efficiency_denominator_eecpt,
    };

    RooUnfoldResponse *response0_eecpt, 
        *response1_eecpt,
        *response2_eecpt,
        *response3_eecpt,
        *response4_eecpt,
        *response5_eecpt,
        *response6_eecpt,
        *response7_eecpt,
        *response8_eecpt,
        *response9_eecpt;

    std::vector<RooUnfoldResponse *> responses_eecpt = {
        response0_eecpt, 
        response1_eecpt,
        response2_eecpt,
        response3_eecpt,
        response4_eecpt,
        response5_eecpt,
        response6_eecpt,
        response7_eecpt,
        response8_eecpt,
        response9_eecpt,
    };
    
    
    // initialize histograms and responses
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_purity_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
        histos_purity_denominator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_purity_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_numerator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_efficiency_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_denominator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_efficiency_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);

        
        responses_eecpt[i] = new RooUnfoldResponse(histos_purity_denominator_eecpt[i],histos_efficiency_denominator_eecpt[i], "response"+TString(Form("%d",i))+"_eecpt", "response for 3d: eec and jet pt"); 
    }        

    // Loop over tree entries
    Long64_t nentries = tree->GetEntries();
    std::cout << "Entries: " << nentries << std::endl;
    for (Long64_t ient = 0; ient < nentries; ient++) {

        //Print progress
        if (ient%1000000==0) cout << "ient=" << ient << std::endl;

        tree->GetEntry(ient);


        Int_t cuts = 1;                                                                                                                                               
            if      (flav == "b1")    cuts = 1;                                                                                                                           
            else if (flav == "b2")    cuts = 2;                                                                                                                           
            else if (flav == "nonb")  cuts = 3;                                                                                                                           
            else if (flav == "all")   cuts = 4;                                                                                                                           
            else if (flav == "c")     cuts = 5;                                                                                                                           
            else if (flav == "light") cuts = 6;                                                                                                                           
            else {                                                                                                                                                      
        std::cerr << "ERROR: Unknown label '" << label << "'\n";                                                                                                  
        exit(1);                                                                                                                                                       
        }       

        bool skip = false;                                                                  
        // Select jet flavour and/or select on the number of b hadrons                                                                                                                                                                                                                                          
            switch(cuts){                                                                                                                                                
        //b-jet with one b hadron                                                                                                                                  
        case 1:                                                                                                                                                       
	  if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                 
	  if (std::abs(jtNbHad) != 1) skip = true;                                                                                                                        	  
        break;                                                                                                                                                       	  
        //b-jet with more than 1 b hadron                                                                                                                           
        case 2:                                                                                                                                                              	  
        if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                          	  
        if (std::abs(jtNbHad) < 2) skip = true;                                                                                                                     	  
        break;                                                                                                                                                                   	  
        //non-b jets                                                                                                                                                                                                                                                                                          
        case 3:                                                                                                                                                              	  
        if (std::abs(jtHadFlav) == 5) skip = true;                                                                                                                         	  
        break;                                                                                                                                                             	  
        //no flavour selection                                                                                                                                                                                                                                                                                
        case 4:                                                                                                                                                              	  
        skip = false;                                                                                                                                                    	  
        break;                                                                                                                                                         	  
        //c-jets                                                                                                                                                     
        case 5:                                                                                                                                                              	  
        if(std::abs(jtHadFlav) != 4) skip = true;                                                                                                                 	  
        break;                                                                                                                                                          	  
        //light (non-b non-c) jets                                                                                                                                                                                                                                                                            
        case 6:                                                                                                                                                       	  
        if(std::abs(jtHadFlav) >= 4) skip = true;                                                                                                                         	  
        break;                                                                                                                                                                  }                                                                                                                                                                     if (skip) continue;      
        /////////////////////////////





        if (skipMC(jpt_reco, jpt_gen, pthat)) continue;

        // Random number for jack-knife resampling (Lida)
        double num = distr(generator);



        // Check if jet has a match at gen
        bool has_gen_match = (jpt_gen > 0);

        // Fill histograms only if the jet as a match at gen
        if (!has_gen_match) {   
            // fill fakes
            continue; 
        } 

        // Skip jets outside tracker 
        if (std::abs(jt_eta_reco) > 1.6) continue;
        if (std::abs(jt_eta_gen) > 1.6) continue;

        //Select jets passing the reco b-jet tagging
        if (btag && std::abs(discr) <= 0.99) continue;


        // The rest of the histograms don't include any fakes


        //Loop over matched dr to fill the eec histograms

        //Debug
        if(ndr_gen != ndr_reco){
            std::cout << "Different dr entries for reco and gen" << std::endl;
        }

        for (Int_t j = 0; j < ndr_reco; j++){

            Float_t dr_reco_j = dr_reco[j];
            Float_t dr_gen_j = dr_gen[j];

            Float_t eec_reco_j = eec_reco[j];
            Float_t eec_gen_j = eec_gen[j];



            num = distr(generator);

            
            //checks for underflow/overflow
            if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
            if(dr_gen_j >= dr_max) dr_gen_j = dr_max_fill;
            if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;
            if(dr_gen_j < dr_min) dr_gen_j = dr_min_fill;
            if(eec_reco_j >= eec_max) eec_reco_j = eec_max_fill;
            if(eec_gen_j >= eec_max) eec_gen_j = eec_max_fill;

            //Entries are/aren't passing cuts at gen or reco level
            bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min && dr_reco_j < dr_max && dr_reco_j >= dr_min && eec_reco_j < eec_max && eec_reco_j >= eec_min);
            bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min && dr_gen_j < dr_max && dr_gen_j >= dr_min && eec_gen_j < eec_max && eec_gen_j >= eec_min);


            // fill eec histograms
            if (true_pass_cuts_eec && reco_pass_cuts_eec) {

                //count matched pairs 
                matched_gen += 1;
                matched_reco += 1;

                //Lida
                fill_jk_resampling(histos_efficiency_numerator_eecpt, num, dr_gen_j, eec_gen_j, jpt_gen, weight);
                fill_jk_resampling(histos_purity_numerator_eecpt, num, dr_reco_j, eec_reco_j, jpt_reco, weight);
                fill_jk_resampling_response(responses_eecpt, num, dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);

                if (num<0.5) {
                    h_half0_efficiency_numerator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                    h_half0_purity_numerator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                    response_half0_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
                } else {
                    h_half1_efficiency_numerator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                    h_half1_purity_numerator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                    response_half1_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
                }
                response_full_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
            }
        }//pair entry loop

        //Loop over all (matched/non matched) dr to fill the eec purity and efficiency

        //For purity
        for (Int_t j = 0; j < ndr_reco_tot; j++){

            Float_t dr_reco_j = dr_reco[j];

            Float_t eec_reco_j = eec_reco[j];



            num = distr(generator);

            
            //checks for underflow/overflow
            if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
            if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;
            if(eec_reco_j >= eec_max) eec_reco_j = eec_max_fill;

            //Check if cuts at reco are passed
            bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min && dr_reco_j < dr_max && dr_reco_j >= dr_min && eec_reco_j < eec_max && eec_reco_j >= eec_min);
            

            // fill purity
            
            if (reco_pass_cuts_eec) {
                fill_jk_resampling(histos_purity_denominator_eecpt, num, dr_reco_j, eec_reco_j, jpt_reco, weight);
                if (num<0.5) h_half0_purity_denominator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                else h_half1_purity_denominator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                if(j >= ndr_reco){
                    notmatched_reco += 1;
                    h_notmatched_reco->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                }
            }

        }// pair entry loop

        //For efficiency
        for (Int_t j = 0; j < ndr_gen_tot; j++){

            Float_t dr_gen_j = dr_gen[j];

            Float_t eec_gen_j = eec_gen[j];



            num = distr(generator);

            
            //checks for underflow/overflow            
            if(dr_gen_j >= dr_max) dr_gen_j = dr_max_fill;
            if(dr_gen_j < dr_min) dr_gen_j = dr_min_fill;
            if(eec_gen_j >= eec_max) eec_gen_j = eec_max_fill;

            //Checks if cuts are passed at gen
            bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min && dr_gen_j < dr_max && dr_gen_j >= dr_min && eec_gen_j < eec_max && eec_gen_j >= eec_min);


            // fill efficiency

            if (true_pass_cuts_eec) {
                fill_jk_resampling(histos_efficiency_denominator_eecpt, num, dr_gen_j, eec_gen_j, jpt_gen, weight);
                if (num<0.5) h_half0_efficiency_denominator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                else h_half1_efficiency_denominator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                if(j >= ndr_gen){
                    notmatched_gen += 1;
                    h_notmatched_gen->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                }
            }

        }// pair entry loop
    } // tree entry loop

    // Create purity and efficiency histograms
    TH3D *h0_purity_eecpt, 
        *h1_purity_eecpt,
        *h2_purity_eecpt,
        *h3_purity_eecpt,
        *h4_purity_eecpt,
        *h5_purity_eecpt,
        *h6_purity_eecpt,
        *h7_purity_eecpt,
        *h8_purity_eecpt,
        *h9_purity_eecpt;

    std::vector<TH3D *> histos_purity_eecpt = {
        h0_purity_eecpt,
        h1_purity_eecpt,
        h2_purity_eecpt,
        h3_purity_eecpt,
        h4_purity_eecpt,
        h5_purity_eecpt,
        h6_purity_eecpt,
        h7_purity_eecpt,
        h8_purity_eecpt,
        h9_purity_eecpt,
    };

    TH3D *h0_efficiency_eecpt, 
        *h1_efficiency_eecpt,
        *h2_efficiency_eecpt,
        *h3_efficiency_eecpt,
        *h4_efficiency_eecpt,
        *h5_efficiency_eecpt,
        *h6_efficiency_eecpt,
        *h7_efficiency_eecpt,
        *h8_efficiency_eecpt,
        *h9_efficiency_eecpt;

    std::vector<TH3D *> histos_efficiency_eecpt = {
        h0_efficiency_eecpt,
        h1_efficiency_eecpt,
        h2_efficiency_eecpt,
        h3_efficiency_eecpt,
        h4_efficiency_eecpt,
        h5_efficiency_eecpt,
        h6_efficiency_eecpt,
        h7_efficiency_eecpt,
        h8_efficiency_eecpt,
        h9_efficiency_eecpt,
    };

    
    
    // initialize histograms
    for (int i=0; i<10; i++) {
        histos_purity_eecpt[i] = (TH3D *) histos_purity_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_eecpt");
        histos_purity_eecpt[i]->Divide(histos_purity_numerator_eecpt[i], histos_purity_denominator_eecpt[i], 1., 1., "b");
        histos_efficiency_eecpt[i] = (TH3D *) histos_efficiency_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_eecpt");
        histos_efficiency_eecpt[i]->Divide(histos_efficiency_numerator_eecpt[i], histos_efficiency_denominator_eecpt[i], 1., 1., "b");
    }

    // declare the per half purity + efficiency histograms
    TH3D *h_half0_purity_eecpt = (TH3D *) h_half0_purity_numerator_eecpt->Clone("h_half0_purity_eecpt");
    h_half0_purity_eecpt->Divide(h_half0_purity_numerator_eecpt, h_half0_purity_denominator_eecpt, 1., 1., "b");
    TH3D *h_half0_efficiency_eecpt = (TH3D *) h_half0_efficiency_numerator_eecpt->Clone("h_half0_efficiency_eecpt");
    h_half0_efficiency_eecpt->Divide(h_half0_efficiency_numerator_eecpt, h_half0_efficiency_denominator_eecpt, 1., 1., "b");

    TH3D *h_half1_purity_eecpt = (TH3D *) h_half1_purity_numerator_eecpt->Clone("h_half1_purity_eecpt");
    h_half1_purity_eecpt->Divide(h_half1_purity_numerator_eecpt, h_half1_purity_denominator_eecpt, 1., 1., "b");
    TH3D *h_half1_efficiency_eecpt = (TH3D *) h_half1_efficiency_numerator_eecpt->Clone("h_half1_efficiency_eecpt");
    h_half1_efficiency_eecpt->Divide(h_half1_efficiency_numerator_eecpt, h_half1_efficiency_denominator_eecpt, 1., 1., "b");

    // declare the full purity + efficiency histograms
    TH3D *h_full_purity_numerator_eecpt = (TH3D *) h_half0_purity_numerator_eecpt->Clone("h_full_purity_numerator_eecpt");
    h_full_purity_numerator_eecpt->Add(h_half1_purity_numerator_eecpt);
    TH3D *h_full_purity_denominator_eecpt = (TH3D *) h_half0_purity_denominator_eecpt->Clone("h_full_purity_denominator_eecpt");
    h_full_purity_denominator_eecpt->Add(h_half1_purity_denominator_eecpt);
    TH3D *h_full_purity_eecpt = (TH3D *) h_full_purity_numerator_eecpt->Clone("h_full_purity_eecpt");
    h_full_purity_eecpt->Divide(h_full_purity_numerator_eecpt, h_full_purity_denominator_eecpt, 1., 1., "b");

    TH3D *h_full_efficiency_numerator_eecpt = (TH3D *) h_half0_efficiency_numerator_eecpt->Clone("h_full_efficiency_numerator_eecpt");
    h_full_efficiency_numerator_eecpt->Add(h_half1_efficiency_numerator_eecpt);
    TH3D *h_full_efficiency_denominator_eecpt = (TH3D *) h_half0_efficiency_denominator_eecpt->Clone("h_full_efficiency_denominator_eecpt");
    h_full_efficiency_denominator_eecpt->Add(h_half1_efficiency_denominator_eecpt);
    TH3D *h_full_efficiency_eecpt = (TH3D *) h_full_efficiency_numerator_eecpt->Clone("h_full_efficiency_eecpt");
    h_full_efficiency_eecpt->Divide(h_full_efficiency_numerator_eecpt, h_full_efficiency_denominator_eecpt, 1., 1., "b");

    // Create output file
    std::cout << "Creating file: " << folder+fout_name << std::endl;
    TFile *fout = new TFile(folder+fout_name, "recreate");

    // Write jk resampling histograms + responses (Lida)
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i]->Write();
        histos_purity_denominator_eecpt[i]->Write();
        histos_purity_eecpt[i]->Write();

        histos_efficiency_numerator_eecpt[i]->Write();
        histos_efficiency_denominator_eecpt[i]->Write();
        histos_efficiency_eecpt[i]->Write();

        responses_eecpt[i]->Write();

    }

    // Write per half histograms 
    h_half0_purity_numerator_eecpt->Write();
    h_half0_purity_denominator_eecpt->Write();
    h_half0_purity_eecpt->Write();

    h_half0_efficiency_numerator_eecpt->Write();
    h_half0_efficiency_denominator_eecpt->Write();
    h_half0_efficiency_eecpt->Write();

    response_half0_eecpt->Write();

    h_half1_purity_numerator_eecpt->Write();
    h_half1_purity_denominator_eecpt->Write();
    h_half1_purity_eecpt->Write();

    h_half1_efficiency_numerator_eecpt->Write();
    h_half1_efficiency_denominator_eecpt->Write();
    h_half1_efficiency_eecpt->Write();

    response_half1_eecpt->Write();

    
    // Write full histograms 
    h_full_purity_numerator_eecpt->Write();
    h_full_purity_denominator_eecpt->Write();
    h_full_purity_eecpt->Write();

    h_full_efficiency_numerator_eecpt->Write();
    h_full_efficiency_denominator_eecpt->Write();
    h_full_efficiency_eecpt->Write();

    response_full_eecpt->Write();

    h_notmatched_reco->Write();
    h_notmatched_gen->Write();

    fout->Close();
    delete fout;

    //Check pair matching efficiency
    std::cout << "Reco pairs not matched = " << notmatched_reco/(notmatched_reco+matched_reco)*100 << " percent" << std::endl;
    std::cout << "Gen pairs not matched = " << notmatched_gen/(notmatched_gen+matched_gen)*100 << " percent" << std::endl;
}

//Creates a 3D response matrix and purity/efficiency corrections from a tree for the inclusive sample
//(or generally for any sample where the tree is split into different files)
void create_response_3D_inclusive(TString &sample, TString &label, TString &folder, bool btag, Int_t &n, Float_t &pT_low, Float_t &pT_high)
{  

    //Create the fout name depending on the selection
    TString fout_name = "histos_response_3D_";

    if(!btag) label += "_notag"; 

    fout_name += TString(Form("n%i_", n)) + sample + "_" + label + ".root";

    std::cout<< "HOLA" << std::endl;
    
    //See max and min of the energy weights (for check)
    Float_t minimum_entry = 0;
    Float_t maximum_entry = 0;

    // random number generator for jackknife resampling
    const double range_from = 0;
    const double range_to = 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(range_from, range_to);

    
    //Check percentage of non-matched tracks
    Double_t matched_reco = 0;
    Double_t matched_gen = 0;
    Double_t notmatched_reco = 0;
    Double_t notmatched_gen = 0;

    
    // Declare histograms
    TH3D *h_half0_purity_numerator_eecpt = new TH3D("h_half0_purity_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half0_purity_denominator_eecpt = new TH3D("h_half0_purity_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half0_efficiency_numerator_eecpt = new TH3D("h_half0_efficiency_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half0_efficiency_denominator_eecpt = new TH3D("h_half0_efficiency_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_purity_numerator_eecpt = new TH3D("h_half1_purity_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_purity_denominator_eecpt = new TH3D("h_half1_purity_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_efficiency_numerator_eecpt = new TH3D("h_half1_efficiency_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_half1_efficiency_denominator_eecpt = new TH3D("h_half1_efficiency_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    // histograms only for the non matched tracks
    TH3D *h_notmatched_reco = new TH3D("h_notmatched_reco", "h_notmatched_reco", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
    TH3D *h_notmatched_gen = new TH3D("h_notmatched_gen", "h_notmatched_gen", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);

    RooUnfoldResponse *response_half0_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half0", "response for 3d: eec and jet pt"); 
    RooUnfoldResponse *response_half1_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_eecpt_half1", "response for 3d: eec and jet pt"); 
    RooUnfoldResponse *response_full_eecpt = new RooUnfoldResponse(h_half0_purity_denominator_eecpt, h_half0_efficiency_denominator_eecpt, "response_full_eecpt", "response for 3d: eec and jet pt"); 


    // rg jk resampling histogram definiton (Lida)
    TH3D *h0_purity_numerator_eecpt, 
        *h1_purity_numerator_eecpt, 
        *h2_purity_numerator_eecpt, 
        *h3_purity_numerator_eecpt, 
        *h4_purity_numerator_eecpt, 
        *h5_purity_numerator_eecpt, 
        *h6_purity_numerator_eecpt, 
        *h7_purity_numerator_eecpt, 
        *h8_purity_numerator_eecpt, 
        *h9_purity_numerator_eecpt;

    std::vector<TH3D *> histos_purity_numerator_eecpt = {
        h0_purity_numerator_eecpt, 
        h1_purity_numerator_eecpt, 
        h2_purity_numerator_eecpt, 
        h3_purity_numerator_eecpt, 
        h4_purity_numerator_eecpt, 
        h5_purity_numerator_eecpt, 
        h6_purity_numerator_eecpt, 
        h7_purity_numerator_eecpt, 
        h8_purity_numerator_eecpt, 
        h9_purity_numerator_eecpt,
    };

    TH3D *h0_purity_denominator_eecpt, 
        *h1_purity_denominator_eecpt, 
        *h2_purity_denominator_eecpt, 
        *h3_purity_denominator_eecpt, 
        *h4_purity_denominator_eecpt, 
        *h5_purity_denominator_eecpt, 
        *h6_purity_denominator_eecpt, 
        *h7_purity_denominator_eecpt, 
        *h8_purity_denominator_eecpt, 
        *h9_purity_denominator_eecpt;

    std::vector<TH3D *> histos_purity_denominator_eecpt = {
        h0_purity_denominator_eecpt, 
        h1_purity_denominator_eecpt, 
        h2_purity_denominator_eecpt, 
        h3_purity_denominator_eecpt, 
        h4_purity_denominator_eecpt, 
        h5_purity_denominator_eecpt, 
        h6_purity_denominator_eecpt, 
        h7_purity_denominator_eecpt, 
        h8_purity_denominator_eecpt, 
        h9_purity_denominator_eecpt,
    };

    TH3D *h0_efficiency_numerator_eecpt, 
        *h1_efficiency_numerator_eecpt, 
        *h2_efficiency_numerator_eecpt, 
        *h3_efficiency_numerator_eecpt, 
        *h4_efficiency_numerator_eecpt, 
        *h5_efficiency_numerator_eecpt, 
        *h6_efficiency_numerator_eecpt, 
        *h7_efficiency_numerator_eecpt, 
        *h8_efficiency_numerator_eecpt, 
        *h9_efficiency_numerator_eecpt;

    std::vector<TH3D *> histos_efficiency_numerator_eecpt = {
        h0_efficiency_numerator_eecpt, 
        h1_efficiency_numerator_eecpt, 
        h2_efficiency_numerator_eecpt, 
        h3_efficiency_numerator_eecpt, 
        h4_efficiency_numerator_eecpt, 
        h5_efficiency_numerator_eecpt, 
        h6_efficiency_numerator_eecpt, 
        h7_efficiency_numerator_eecpt, 
        h8_efficiency_numerator_eecpt, 
        h9_efficiency_numerator_eecpt,
    };

    TH3D *h0_efficiency_denominator_eecpt, 
        *h1_efficiency_denominator_eecpt, 
        *h2_efficiency_denominator_eecpt, 
        *h3_efficiency_denominator_eecpt, 
        *h4_efficiency_denominator_eecpt, 
        *h5_efficiency_denominator_eecpt, 
        *h6_efficiency_denominator_eecpt, 
        *h7_efficiency_denominator_eecpt, 
        *h8_efficiency_denominator_eecpt, 
        *h9_efficiency_denominator_eecpt;

    std::vector<TH3D *> histos_efficiency_denominator_eecpt = {
        h0_efficiency_denominator_eecpt, 
        h1_efficiency_denominator_eecpt, 
        h2_efficiency_denominator_eecpt, 
        h3_efficiency_denominator_eecpt, 
        h4_efficiency_denominator_eecpt, 
        h5_efficiency_denominator_eecpt, 
        h6_efficiency_denominator_eecpt, 
        h7_efficiency_denominator_eecpt, 
        h8_efficiency_denominator_eecpt, 
        h9_efficiency_denominator_eecpt,
    };

    RooUnfoldResponse *response0_eecpt, 
        *response1_eecpt,
        *response2_eecpt,
        *response3_eecpt,
        *response4_eecpt,
        *response5_eecpt,
        *response6_eecpt,
        *response7_eecpt,
        *response8_eecpt,
        *response9_eecpt;

    std::vector<RooUnfoldResponse *> responses_eecpt = {
        response0_eecpt, 
        response1_eecpt,
        response2_eecpt,
        response3_eecpt,
        response4_eecpt,
        response5_eecpt,
        response6_eecpt,
        response7_eecpt,
        response8_eecpt,
        response9_eecpt,
    };
    
    
    // initialize histograms and responses
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_purity_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
        histos_purity_denominator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_purity_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_numerator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_efficiency_numerator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);
        histos_efficiency_denominator_eecpt[i] = new TH3D("h"+TString(Form("%d",i))+"_efficiency_denominator_eecpt", "x=dr, y=eec, z=jtpt", dr_bins, dr_binsVector, eec_bins, eec_binsVector, jtpt_bins, jtpt_binsVector);

        
        responses_eecpt[i] = new RooUnfoldResponse(histos_purity_denominator_eecpt[i],histos_efficiency_denominator_eecpt[i], "response"+TString(Form("%d",i))+"_eecpt", "response for 3d: eec and jet pt"); 
    }        


    // Loop over tree files, the tree is split in 30 files because it is too big
    for(Int_t i = 0; i <= 0; i++){
        //Create the fout name depending on the selection
        TString fin_name =  sample;//"trees_nocuts_matched_noaggr_";

        //fin_name += TString(Form("n%i_",n)) + label + "_" + sample + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high)));

        //if(i != 0) fin_name += Form("_%i.root", i);
        //else fin_name += ".root";

        std::cout << "fin: " << folder+fin_name << std::endl;
        TFile *fin = new TFile(folder+fin_name);
    
        TString tree_name = "tree";
        std::cout << "tree: " << tree_name << std::endl;
        TTree *tree = (TTree *) fin->Get(tree_name);

        // Set tree addresses
        
        
        Int_t ndr_reco_tot;
        Int_t ndr_gen_tot;
        
	Int_t ndr_reco, ndr_gen, jtHadFlav, jtNbHad;                                                                                             
	Float_t weight, pthat, jpt_reco, jpt_gen, mb_reco, mb_gen, jt_eta_reco, jt_eta_gen, discr;                                               
	Float_t dr_reco[6000], dr_gen[6000], eec_reco[6000], eec_gen[6000];                                                                      
        
    
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("pthat", &pthat);
        tree->SetBranchAddress("ndr_reco", &ndr_reco);
        tree->SetBranchAddress("ndr_gen", &ndr_gen);
        tree->SetBranchAddress("ndr_reco_tot", &ndr_reco_tot);
        tree->SetBranchAddress("ndr_gen_tot", &ndr_gen_tot);
        tree->SetBranchAddress("jpt_reco", &jpt_reco);
        tree->SetBranchAddress("jpt_gen", &jpt_gen);
        tree->SetBranchAddress("dr_reco", &dr_reco);
        tree->SetBranchAddress("dr_gen", &dr_gen);
        tree->SetBranchAddress("eec_reco", &eec_reco);
        tree->SetBranchAddress("eec_gen", &eec_gen);
        tree->SetBranchAddress("jt_eta_reco", &jt_eta_reco);
        tree->SetBranchAddress("jt_eta_gen", &jt_eta_gen);
        tree->SetBranchAddress("discr", &discr);
	    tree->SetBranchAddress("jtHadFlav", &jtHadFlav);                                                                                         
	    tree->SetBranchAddress("jtNbHad", &jtNbHad);                                                                                             
         


        Long64_t nentries = tree->GetEntries();
        std::cout << "Entries: " << nentries << " for i = "<< i << std::endl;
        for (Long64_t ient = beg_event; ient < end_event; ient++) {
            if (ient%1000000==0) cout << "ient=" << ient << std::endl; 
            tree->GetEntry(ient);

//NO flav selection

            if (skipMC(jpt_reco, jpt_gen, pthat)) continue;
    

            // Random number for jack-knife resampling (Lida)
            double num = distr(generator);
    
    
            // Check if pass cuts
            bool has_gen_match = (jpt_gen > 0);
    
            // Fill histograms if the jet has a match at gen
            if (!has_gen_match) {   
                // fill fakes
                continue; 
            } 
    
            // Skip jets outside tracker 
            if (std::abs(jt_eta_reco) > 1.6) continue;
            if (std::abs(jt_eta_gen) > 1.6) continue;
    
            //Select jets passing the reco b-jet tagging 
            if (btag && std::abs(discr) <= 0.99) continue;
    
    
            // The rest of the histograms don;t include any fakes
    
    
            //Loop over dr to fill the eec histograms

            //debug
            if(ndr_gen != ndr_reco){
                std::cout << "Different dr entries for reco and gen" << std::endl;
            }
    
            //Loop over matched pairs
            for (Int_t j = 0; j < ndr_reco; j++){
    
                Float_t dr_reco_j = dr_reco[j];
                Float_t dr_gen_j = dr_gen[j];
    
                Float_t eec_reco_j = eec_reco[j];
                Float_t eec_gen_j = eec_gen[j];
    
    
    
                num = distr(generator);
    
                
                //checks for underflow/overflow
                if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
                if(dr_gen_j >= dr_max) dr_gen_j = dr_max_fill;
                if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;
                if(dr_gen_j < dr_min) dr_gen_j = dr_min_fill;
                if(eec_reco_j >= eec_max) eec_reco_j = eec_max_fill;
                if(eec_gen_j >= eec_max) eec_gen_j = eec_max_fill;
    
                bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min && dr_reco_j < dr_max && dr_reco_j >= dr_min && eec_reco_j < eec_max && eec_reco_j >= eec_min);
                bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min && dr_gen_j < dr_max && dr_gen_j >= dr_min && eec_gen_j < eec_max && eec_gen_j >= eec_min);
    
    
                // fill eec histograms
                if (true_pass_cuts_eec && reco_pass_cuts_eec) {
                    matched_gen += 1;
                    matched_reco += 1;
                    
                    fill_jk_resampling(histos_efficiency_numerator_eecpt, num, dr_gen_j, eec_gen_j, jpt_gen, weight);
                    fill_jk_resampling(histos_purity_numerator_eecpt, num, dr_reco_j, eec_reco_j, jpt_reco, weight);
                    fill_jk_resampling_response(responses_eecpt, num, dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
    
                    if (num<0.5) {
                        h_half0_efficiency_numerator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                        h_half0_purity_numerator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                        response_half0_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
                    } else {
                        h_half1_efficiency_numerator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                        h_half1_purity_numerator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                        response_half1_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
                    }
                    response_full_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, dr_gen_j, eec_gen_j, jpt_gen, weight);
                }
            }// pair entry loop
    
            //Loop over dr to fill the eec purity and efficiency

            //Loop over all pairs (matched and non-matched)
            for (Int_t j = 0; j < ndr_reco_tot; j++){
    
                Float_t dr_reco_j = dr_reco[j];
    
                Float_t eec_reco_j = eec_reco[j];
    
                //Find max and mix for check
                if(dr_reco_j < minimum_entry) minimum_entry = dr_reco_j;
                if(eec_reco_j > maximum_entry) maximum_entry = eec_reco_j;
    
    
                num = distr(generator);
    
                
                //checks for underflow/overflow
                if(dr_reco_j >= dr_max) dr_reco_j = dr_max_fill;
                if(dr_reco_j < dr_min) dr_reco_j = dr_min_fill;
                if(eec_reco_j >= eec_max) eec_reco_j = eec_max_fill;
    
                bool reco_pass_cuts_eec = (jpt_reco < jtpt_max && jpt_reco >= jtpt_min && dr_reco_j < dr_max && dr_reco_j >= dr_min && eec_reco_j < eec_max && eec_reco_j >= eec_min);
                
    
                // fill purity
                
                if (reco_pass_cuts_eec) {
                    fill_jk_resampling(histos_purity_denominator_eecpt, num, dr_reco_j, eec_reco_j, jpt_reco, weight);
                    if (num<0.5) h_half0_purity_denominator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                    else h_half1_purity_denominator_eecpt->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                    if(j >= ndr_reco){
                        notmatched_reco += 1;
                        h_notmatched_reco->Fill(dr_reco_j, eec_reco_j, jpt_reco, weight);
                    }
                }
            }// pair entry loop
    
            for (Int_t j = 0; j < ndr_gen_tot; j++){
    
                Float_t dr_gen_j = dr_gen[j];
    
                Float_t eec_gen_j = eec_gen[j];
    
    
    
                num = distr(generator);
    
                
                //checks for underflow/overflow            
                if(dr_gen_j >= dr_max) dr_gen_j = dr_max_fill;
                if(dr_gen_j < dr_min) dr_gen_j = dr_min_fill;
                if(eec_gen_j >= eec_max) eec_gen_j = eec_max_fill;
    
                bool true_pass_cuts_eec = (jpt_gen < jtpt_max && jpt_gen >= jtpt_min && dr_gen_j < dr_max && dr_gen_j >= dr_min && eec_gen_j < eec_max && eec_gen_j >= eec_min);
    
    
    
                // fill eec histograms
                if (true_pass_cuts_eec) {
                    fill_jk_resampling(histos_efficiency_denominator_eecpt, num, dr_gen_j, eec_gen_j, jpt_gen, weight);
                    if (num<0.5) h_half0_efficiency_denominator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                    else h_half1_efficiency_denominator_eecpt->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                    if(j >= ndr_gen){
                        notmatched_gen += 1;
                        h_notmatched_gen->Fill(dr_gen_j, eec_gen_j, jpt_gen, weight);
                    }
                }
            }// pair entry loop
        } // tree entry loop
    }

    // Create purity and efficiency histograms
    TH3D *h0_purity_eecpt, 
        *h1_purity_eecpt,
        *h2_purity_eecpt,
        *h3_purity_eecpt,
        *h4_purity_eecpt,
        *h5_purity_eecpt,
        *h6_purity_eecpt,
        *h7_purity_eecpt,
        *h8_purity_eecpt,
        *h9_purity_eecpt;

    std::vector<TH3D *> histos_purity_eecpt = {
        h0_purity_eecpt,
        h1_purity_eecpt,
        h2_purity_eecpt,
        h3_purity_eecpt,
        h4_purity_eecpt,
        h5_purity_eecpt,
        h6_purity_eecpt,
        h7_purity_eecpt,
        h8_purity_eecpt,
        h9_purity_eecpt,
    };

    TH3D *h0_efficiency_eecpt, 
        *h1_efficiency_eecpt,
        *h2_efficiency_eecpt,
        *h3_efficiency_eecpt,
        *h4_efficiency_eecpt,
        *h5_efficiency_eecpt,
        *h6_efficiency_eecpt,
        *h7_efficiency_eecpt,
        *h8_efficiency_eecpt,
        *h9_efficiency_eecpt;

    std::vector<TH3D *> histos_efficiency_eecpt = {
        h0_efficiency_eecpt,
        h1_efficiency_eecpt,
        h2_efficiency_eecpt,
        h3_efficiency_eecpt,
        h4_efficiency_eecpt,
        h5_efficiency_eecpt,
        h6_efficiency_eecpt,
        h7_efficiency_eecpt,
        h8_efficiency_eecpt,
        h9_efficiency_eecpt,
    };

    
    
    // initialize histograms
    for (int i=0; i<10; i++) {
        histos_purity_eecpt[i] = (TH3D *) histos_purity_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_eecpt");
        histos_purity_eecpt[i]->Divide(histos_purity_numerator_eecpt[i], histos_purity_denominator_eecpt[i], 1., 1., "b");
        histos_efficiency_eecpt[i] = (TH3D *) histos_efficiency_numerator_eecpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_eecpt");
        histos_efficiency_eecpt[i]->Divide(histos_efficiency_numerator_eecpt[i], histos_efficiency_denominator_eecpt[i], 1., 1., "b");
    }

    // declare the per half purity + efficiency histograms
    TH3D *h_half0_purity_eecpt = (TH3D *) h_half0_purity_numerator_eecpt->Clone("h_half0_purity_eecpt");
    h_half0_purity_eecpt->Divide(h_half0_purity_numerator_eecpt, h_half0_purity_denominator_eecpt, 1., 1., "b");
    TH3D *h_half0_efficiency_eecpt = (TH3D *) h_half0_efficiency_numerator_eecpt->Clone("h_half0_efficiency_eecpt");
    h_half0_efficiency_eecpt->Divide(h_half0_efficiency_numerator_eecpt, h_half0_efficiency_denominator_eecpt, 1., 1., "b");

    TH3D *h_half1_purity_eecpt = (TH3D *) h_half1_purity_numerator_eecpt->Clone("h_half1_purity_eecpt");
    h_half1_purity_eecpt->Divide(h_half1_purity_numerator_eecpt, h_half1_purity_denominator_eecpt, 1., 1., "b");
    TH3D *h_half1_efficiency_eecpt = (TH3D *) h_half1_efficiency_numerator_eecpt->Clone("h_half1_efficiency_eecpt");
    h_half1_efficiency_eecpt->Divide(h_half1_efficiency_numerator_eecpt, h_half1_efficiency_denominator_eecpt, 1., 1., "b");

    // declare the full purity + efficiency histograms
    TH3D *h_full_purity_numerator_eecpt = (TH3D *) h_half0_purity_numerator_eecpt->Clone("h_full_purity_numerator_eecpt");
    h_full_purity_numerator_eecpt->Add(h_half1_purity_numerator_eecpt);
    TH3D *h_full_purity_denominator_eecpt = (TH3D *) h_half0_purity_denominator_eecpt->Clone("h_full_purity_denominator_eecpt");
    h_full_purity_denominator_eecpt->Add(h_half1_purity_denominator_eecpt);
    TH3D *h_full_purity_eecpt = (TH3D *) h_full_purity_numerator_eecpt->Clone("h_full_purity_eecpt");
    h_full_purity_eecpt->Divide(h_full_purity_numerator_eecpt, h_full_purity_denominator_eecpt, 1., 1., "b");

    TH3D *h_full_efficiency_numerator_eecpt = (TH3D *) h_half0_efficiency_numerator_eecpt->Clone("h_full_efficiency_numerator_eecpt");
    h_full_efficiency_numerator_eecpt->Add(h_half1_efficiency_numerator_eecpt);
    TH3D *h_full_efficiency_denominator_eecpt = (TH3D *) h_half0_efficiency_denominator_eecpt->Clone("h_full_efficiency_denominator_eecpt");
    h_full_efficiency_denominator_eecpt->Add(h_half1_efficiency_denominator_eecpt);
    TH3D *h_full_efficiency_eecpt = (TH3D *) h_full_efficiency_numerator_eecpt->Clone("h_full_efficiency_eecpt");
    h_full_efficiency_eecpt->Divide(h_full_efficiency_numerator_eecpt, h_full_efficiency_denominator_eecpt, 1., 1., "b");

    // Create output file
    std::cout << "Creating file: " << folder + fout_name << std::endl;
    TFile *fout = new TFile(folder + fout_name, "recreate");

    // Write jk resampling histograms + responses (Lida)
    for (int i=0; i<10; i++) {
        histos_purity_numerator_eecpt[i]->Write();
        histos_purity_denominator_eecpt[i]->Write();
        histos_purity_eecpt[i]->Write();

        histos_efficiency_numerator_eecpt[i]->Write();
        histos_efficiency_denominator_eecpt[i]->Write();
        histos_efficiency_eecpt[i]->Write();

        responses_eecpt[i]->Write();

    }

    // Write per half histograms 
    h_half0_purity_numerator_eecpt->Write();
    h_half0_purity_denominator_eecpt->Write();
    h_half0_purity_eecpt->Write();

    h_half0_efficiency_numerator_eecpt->Write();
    h_half0_efficiency_denominator_eecpt->Write();
    h_half0_efficiency_eecpt->Write();

    response_half0_eecpt->Write();

    h_half1_purity_numerator_eecpt->Write();
    h_half1_purity_denominator_eecpt->Write();
    h_half1_purity_eecpt->Write();

    h_half1_efficiency_numerator_eecpt->Write();
    h_half1_efficiency_denominator_eecpt->Write();
    h_half1_efficiency_eecpt->Write();

    response_half1_eecpt->Write();

    
    // Write full histograms 
    h_full_purity_numerator_eecpt->Write();
    h_full_purity_denominator_eecpt->Write();
    h_full_purity_eecpt->Write();

    h_full_efficiency_numerator_eecpt->Write();
    h_full_efficiency_denominator_eecpt->Write();
    h_full_efficiency_eecpt->Write();

    response_full_eecpt->Write();

    h_notmatched_reco->Write();
    h_notmatched_gen->Write();

    std::cout << "min dr = " << minimum_entry << ", max eec = " << maximum_entry << std::endl;

    fout->Close();
    delete fout;

    std::cout << "Reco pairs not matched = " << notmatched_reco/(notmatched_reco+matched_reco)*100 << " percent" << std::endl;
    std::cout << "Gen pairs not matched = " << notmatched_gen/(notmatched_gen+matched_gen)*100 << " percent" << std::endl;

}

//________________________________________________________________________________________________________
//_______________________Check dr and eec bin migrations__________________________________________________
//________________________________________________________________________________________________________

//Draws the migration (so you can just plot it without rerunning the histogram filling)
void draw_migration(RooUnfoldResponse* &response, TString &observable, TString &filename,  TString &sample, TString &label, TString &folder){
   
    //Define the canvas
    TCanvas *c = new TCanvas("c", " ",500,500,904,804);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetRightMargin(0.2);

    //Get response matrix as a 2D histogram
    TMatrixD response_matrix = response->Mresponse();
    TH2D *h = new TH2D(response_matrix);

    //Draw histograms
    h->SetStats(0);
    if(observable == "dr") h->SetTitle("\\mbox{Effect of bin-to-bin migration for the }\\Delta\\mbox{r distribution}");
    else h->SetTitle("Effect of bin-to-bin migration for the " + observable + " distribution");
    h->GetZaxis()->SetTitle("Migration probability");
    h->GetZaxis()->SetTitleOffset(1.5);
    if(observable == "dr"){
        h->GetXaxis()->SetTitle("\\Delta\\mbox{r at detector level}");
        h->GetYaxis()->SetTitle("\\Delta\\mbox{r at particle level}");
    }
    else{
        h->GetXaxis()->SetTitle(observable + " at detector level");
        h->GetYaxis()->SetTitle(observable + " at particle level");
    }
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    gStyle->SetPaintTextFormat("4.2f");
    h->Draw("colz text");

    //Save plot as pdf
    c->Print(folder + observable +"_migration_effect_" + sample + "_" + label + ".pdf");
    //c->Print(folder + observable +"_migration_effect_" + sample + "_" + label + ".cpp");
}

//Fills the dr and eec bin migration histograms (already normalized as "response matrices")
void get_eec_dr_migration(TString &filename, TString &sample, TString &label, TString &folder, bool &btag){
    TString flav = label; 
    std::cout << "flav:" << flav << std::endl;

    
    TString fin_name = filename; 

    std::cout << "fin: " << fin_name << std::endl;
    TFile *fin = new TFile(fin_name);

    TString tree_name = "tree";

    std::cout << "tree: " << tree_name << std::endl;
    TTree *tree = (TTree *) fin->Get(tree_name);

    // Set tree addresses
    Int_t evt_nr;
    Double_t weight;
    Double_t pthat;
    Int_t ndr_reco;
    Int_t ndr_gen;
    Int_t ntrk_reco;
    Int_t ntrk_gen;
    Int_t passcuts_reco;
    Int_t passcuts_gen;
    Int_t njet_reco;
    Int_t njet_gen;
    Int_t jt_index_reco;
    Int_t jt_index_gen;
    Double_t jpt_reco;
    Double_t jpt_gen;
    Float_t mb_reco;
    Float_t mb_gen;
    Float_t dr_reco[4000];
    Float_t dr_gen[4000];
    Float_t eec_reco[4000];
    Float_t eec_gen[4000];
    Double_t jt_eta_reco;
    Double_t jt_eta_gen;
    Double_t discr;
    Int_t jtHadFlav;
    Int_t jtNbHad;


    //tree->SetBranchAddress("evt_nr", &evt_nr);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("pthat", &pthat);
    tree->SetBranchAddress("ndr_reco", &ndr_reco);
    tree->SetBranchAddress("ndr_gen", &ndr_gen);
    tree->SetBranchAddress("ntrk_reco", &ntrk_reco);
    tree->SetBranchAddress("ntrk_gen", &ntrk_gen);
    tree->SetBranchAddress("passcuts_reco", &passcuts_reco);
    tree->SetBranchAddress("passcuts_gen", &passcuts_gen);
    tree->SetBranchAddress("njet_reco", &njet_reco);
    tree->SetBranchAddress("njet_gen", &njet_gen);
    tree->SetBranchAddress("jt_index_reco", &jt_index_reco);
    tree->SetBranchAddress("jt_index_gen", &jt_index_gen);
    tree->SetBranchAddress("jpt_reco", &jpt_reco);
    tree->SetBranchAddress("jpt_gen", &jpt_gen);
    tree->SetBranchAddress("mB_reco", &mb_reco);
    tree->SetBranchAddress("mB_gen", &mb_gen);
    tree->SetBranchAddress("dr_reco", &dr_reco);
    tree->SetBranchAddress("dr_gen", &dr_gen);
    tree->SetBranchAddress("eec_reco", &eec_reco);
    tree->SetBranchAddress("eec_gen", &eec_gen);
    tree->SetBranchAddress("jt_eta_reco", &jt_eta_reco);
    tree->SetBranchAddress("jt_eta_gen", &jt_eta_gen);
    tree->SetBranchAddress("discr", &discr);
    tree->SetBranchAddress("jtHadFlav", &jtHadFlav);                                                                                         
	tree->SetBranchAddress("jtNbHad", &jtNbHad);       

    //Define dr and eec histograms and response matrices
    TH1D *h_dr_migration = new TH1D("h_dr_migration", "h_dr_migration", dr_bins, dr_binsVector);
    RooUnfoldResponse *response_dr = new RooUnfoldResponse(h_dr_migration, h_dr_migration, "response_dr", "response for 1d: dr"); 

    TH1D *h_eec_migration = new TH1D("h_eec_migration", "h_eec_migration", eec_bins, eec_binsVector);
    RooUnfoldResponse *response_eec = new RooUnfoldResponse(h_eec_migration, h_eec_migration, "response_eec", "response for 1d: eec"); 

    // Loop over tree entries
    Long64_t nentries = tree->GetEntries();
    std::cout << "Entries: " << nentries << std::endl;
    for (Long64_t ient = 0; ient < nentries; ient++) {
        //Print progress
        if (ient%1000000==0) cout << "ient=" << ient << std::endl; 

        tree->GetEntry(ient);
        
        Int_t cuts = 1;                                                                                                                                               
        if      (flav == "b1")    cuts = 1;                                                                                                                           
        else if (flav == "b2")    cuts = 2;                                                                                                                           
        else if (flav == "nonb")  cuts = 3;                                                                                                                           
        else if (flav == "all")   cuts = 4;                                                                                                                           
        else if (flav == "c")     cuts = 5;                                                                                                                           
        else if (flav == "light") cuts = 6;                                                                                                                           
        else {                                                                                                                                                      
        std::cerr << "ERROR: Unknown label '" << label << "'\n";                                                                                                  
        exit(1);                                                                                                                                                       
        }       

        bool skip = false;                                                                  
        // Select jet flavour and/or select on the number of b hadrons                                                                                                                                                                                                                                          
        switch(cuts){                                                                                                                                                
        //b-jet with one b hadron                                                                                                                                  
        case 1:                                                                                                                                                               	  
        if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                         	  
        if (std::abs(jtNbHad) != 1) skip = true;                                                                                                                        	  
        break;                                                                                                                                                       	  
        //b-jet with more than 1 b hadron                                                                                                                           
        case 2:                                                                                                                                                              	  
        if (std::abs(jtHadFlav) < 5) skip = true;                                                                                                                          	  
        if (std::abs(jtNbHad) < 2) skip = true;                                                                                                                     	  
        break;                                                                                                                                                                   	  
        //non-b jets                                                                                                                                                                                                                                                                                          
        case 3:                                                                                                                                                              	  
        if (std::abs(jtHadFlav) == 5) skip = true;                                                                                                                         	  
        break;                                                                                                                                                             	  
        //no flavour selection                                                                                                                                                                                                                                                                                
        case 4:                                                                                                                                                              	  
        skip = false;                                                                                                                                                    	  
        break;                                                                                                                                                         	  
        //c-jets                                                                                                                                                     
        case 5:                                                                                                                                                              	  
        if(std::abs(jtHadFlav) != 4) skip = true;                                                                                                                 	  
        break;                                                                                                                                                          	  
        //light (non-b non-c) jets                                                                                                                                                                                                                                                                            
        case 6:                                                                                                                                                       	  
        if(std::abs(jtHadFlav) >= 4) skip = true;                                                                                                                         	  
        break;                                                                                                                                                                  
        }                                                                                                                                                                     
        if (skip) continue;      
        /////////////////////////////

        if (skipMC(jpt_reco, jpt_gen, pthat)) continue;

        // Check if pass cuts
        bool has_gen_match = (jpt_gen > 0);

        // Fill histograms
        if (!has_gen_match) {   
            // fill fakes
            continue; 
        } 
        
        // Skip jets outside tracker 
        if (std::abs(jt_eta_reco) > 1.6) continue;
        if (std::abs(jt_eta_gen) > 1.6) continue;

        //Select jets passing the reco b-jet tagging (to see the bias that the tagging introduces on the measurements)
        if (btag && std::abs(discr) <= 0.99) continue;


        
        // The rest of the histograms don;t include any fakes

        //Debug
        if(ndr_gen != ndr_reco) std::cout << "!! Different ndr_reco and ndr_gen" << std::endl;
        
        //Loop over matched dr to fill the dr and eec histograms
        for (Int_t j = 0; j < ndr_reco; j++){
            Float_t dr_reco_j = dr_reco[j];
            Float_t dr_gen_j = dr_gen[j];
            
            Float_t eec_reco_j = eec_reco[j];
            Float_t eec_gen_j = eec_gen[j];

            response_dr->Fill(dr_reco_j, dr_gen_j, weight);
            response_eec->Fill(eec_reco_j, eec_gen_j, weight);
        }
    }

    //Draw dr migration
    TString observable = "dr";
    draw_migration(response_dr, observable, filename, sample, label, folder);

    //Save dr migration histograms
    TString fout_dr_name = "dr_migration_hist_" + sample + "_" + label + ".root";
    TFile *fout_dr = new TFile(folder + fout_dr_name, "recreate");
    response_dr->Write();
    fout_dr->Close();

    //Draw eec histogram
    observable = "eec";
    draw_migration(response_eec, observable, filename, sample, label, folder);

    //Save eec migration histogram
    TString fout_eec_name = "eec_migration_hist_" + sample + "_" + label + ".root";
    TFile *fout_eec = new TFile(folder+fout_eec_name, "recreate");
    response_eec->Write();
    fout_eec->Close();

}

void create_response(Int_t &beg_event, Int_t &end_event){
  TString datasets{"dijet"}; //{"bjet"};//, "dijet"};
    
  TString folder = "/data_CMS/cms/zaidan/eec_trees/fran_bins/";               //$mydata/eec_trees/merged_trees/";
  Float_t pT_low = 80;
  Float_t pT_high = 140;
  TString pT_selection = "80_140";
  bool btag =false;
  Int_t n=1;

  //Create labels
  std::vector<TString> labels_vec{"b1"};//, "moreb", "other","mc"};

  //    for(Int_t i = 0; i < datasets.size(); i++){
  //        for(Int_t j = 0; j < labels_vec.size(); j++){
  TString filename = "merged_trees_nocuts_matched_noaggr_n1_" + datasets + "_inclusive_"+ pT_selection + ".root";
  std::cout << "Processing file: " << filename << std::endl;
  //	        create_response_1D(filename,  datasets.at(i), labels_vec.at(j), folder, btag, n, pT_low, pT_high);
  //create_response_2D(filename,  datasets.at(i), labels_vec.at(j), folder, btag, n, pT_low, pT_high);
  //create_response_3D(filename,  datasets.at(i), labels_vec.at(j), folder, btag, n, pT_low, pT_high);
  create_response_3D_inclusive(filename,  datasets, folder, btag, n, pT_low, pT_high, beg_event, end_event);
    
  //get_eec_dr_migration(filename, datasets.at(i), labels_vec.at(j), folder, btag);
}

//  }
//}
