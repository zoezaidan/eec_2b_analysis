// Define binnings

// 
//dr
const Int_t dr_binsVectorSize = 21;
Int_t bins_dr = dr_binsVectorSize - 1;
Int_t dr_bins =  dr_binsVectorSize - 1;
// Equally sized bins in logarithmic scale... Zoe binning instead of Frac.
Double_t dr_binsVector[dr_binsVectorSize]= 
{0.002, 
0.00272884, 
0.00372329, 
0.00508014, 
0.00693145, 
0.00945742, 
0.0129039, 
0.0176064, 
0.0240225, 
0.0327768, 
0.0447214, 
0.0610188, 
0.0832553, 
0.113595, 
0.154992, 
0.211474, 
0.28854, 
0.39369, 
0.537159, 
0.732911, 
1};

Double_t dr_min = dr_binsVector[0];
Double_t dr_max = dr_binsVector[bins_dr];

//Values for histogram filling
// Franc. first bin was 0.15 to 0.4
Double_t dr_shiftbin = 0.0001; //  0.0001
Double_t dr_min_fill = 0.0021; //  0.016 = -0.15 + 0.01 
Double_t dr_max_fill = 0.99;// 0.39 = 0.4- 0.01

//mB
const Int_t mb_binsVectorSize = 13;
Int_t mb_bins = mb_binsVectorSize - 1;
Double_t mb_binsVector[mb_binsVectorSize] = {
    0., 
    1.,
    2.,
    3.,
    4.,
    5.,
    6.,
    7.,
    8.,
    9.,
    10.,
    11.,
    12.,
};
Double_t mb_min = mb_binsVector[0];
Double_t mb_max = mb_binsVector[mb_bins];
Int_t bins_mb = 12;
Double_t mb_max_fill = 11.9;


//EEC
const Int_t eec_binsVectorSize = 12; // 6
Int_t bins_eec = eec_binsVectorSize - 1;
Int_t eec_bins = bins_eec;
Double_t eec_binsVector[eec_binsVectorSize] = {
    0., 
    100.,
    200.,
    300.,
    400.,
    500.,
    600.,
    700.,
    800.,
    900.,
    1000.,
    1e+08
};
Double_t eec_max = eec_binsVector[bins_eec];
Double_t eec_min = eec_binsVector[0];
Double_t eec_step = (eec_max-eec_min)/bins_eec;
Double_t eec_max_fill = 1e+08 -1; //499 

//Jet Pt
const Int_t jtpt_binsVectorSize = 4;
Int_t jtpt_bins = jtpt_binsVectorSize - 1;
Double_t jtpt_binsVector[jtpt_binsVectorSize] = {
    80., 
    100., 
    120.,
    140.
};
Double_t jtpt_min = jtpt_binsVector[0];
Double_t jtpt_max = jtpt_binsVector[jtpt_bins];


Int_t bins_pt = 3;
Double_t pt_min = 80;
Double_t pt_max = 140;

//Get the dimension label
Int_t mb_dim = 0;
Int_t dr_dim = 1;
Int_t eec_dim = 2;
Int_t pt_dim = 3;

//Recover the weighted dr distribution from a 2D distribution where dr and eec axes are separated
void recover_eec_distr(TH1D* &h_1D, TH2D* &h, Double_t last_eec = 1000){

    //save the eec along the bins (evaluate at the center of the bin)
    Float_t eec_step = (eec_max-eec_min)/eec_bins;
    Float_t eec = 0;
    Int_t eec_bins_tot = eec_bins;
    
    //save bin entries for Fill()
    Float_t mB, dr, x;

    //For each dr fill in each eec bin the corresponding bin content times the eec weight (at the center of the bin)
    for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){

        dr = h->GetXaxis()->GetBinLowEdge(ibin_dr)+dr_shiftbin; //shift by a small amount to make sure the correct bin is filled
        Float_t dr_error = 0;

        for(Int_t ibin_eec = 1; ibin_eec <= eec_bins_tot; ibin_eec++){

            eec = h->GetYaxis()->GetBinLowEdge(ibin_eec)+eec_step/2;
            //If last eec bin fill with a higher weight
            if(eec_bins_tot == ibin_eec) eec = last_eec;
            
            //Save the bin errors (added in quadrature)
            Float_t dr_err_temp = h->GetBinError(ibin_dr, ibin_eec)*eec;
            dr_error += pow(dr_err_temp, 2);

            h_1D->Fill(dr, h->GetBinContent(ibin_dr, ibin_eec)*eec);

        }

        //Approximate the bin error        
        h_1D->SetBinError(ibin_dr, std::sqrt(dr_error));

    }

}


