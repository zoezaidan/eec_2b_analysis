// Define binnings

//dr
const Int_t dr_binsVectorSize = 21;
Int_t bins_dr = dr_binsVectorSize - 1;
Int_t dr_bins = bins_dr;
Double_t dr_binsVector[dr_binsVectorSize];
void build_dr_bins()
{   
 const Double_t start = 0.002;
 const Double_t end   = 1;

    Double_t log_start = std::log(start);
    Double_t log_end   = std::log(end);
    Double_t step = (log_end - log_start) / (dr_binsVectorSize - 1);

    for (Int_t i = 0; i < dr_binsVectorSize; ++i)
    {
        dr_binsVector[i] = std::exp(log_start + step * i);
    }
}
bool _dr_bins_init = (build_dr_bins(), true);
Double_t dr_min = dr_binsVector[0];
Double_t dr_max = dr_binsVector[bins_dr];

//Values for histogram filling
Double_t dr_shiftbin = 0.00001;
Double_t dr_min_fill = 0.0021;
Double_t dr_max_fill = 0.99;

//mB
const Int_t mb_binsVectorSize = 9;
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
};
Double_t mb_min = mb_binsVector[0];
Double_t mb_max = mb_binsVector[mb_bins];
Int_t bins_mb = 7;
Double_t mb_max_fill = 6.9;

//EEC
const Int_t eec_binsVectorSize = 6;
Int_t bins_eec = eec_binsVectorSize - 1;
Int_t eec_bins = bins_eec;
Double_t eec_binsVector[eec_binsVectorSize] = {
    0., 
    100.,
    200.,
    300.,
    400.,
    500.,
    /*    600.,
    700.,
    800.,
    900.,
    1000., */
};
Double_t eec_max = eec_binsVector[bins_eec];
Double_t eec_min = eec_binsVector[0];
Double_t eec_step = (eec_max-eec_min)/bins_eec;
Double_t eec_max_fill = 499;

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

        //dr = h->GetXaxis()->GetBinCenter(ibin_dr);
        //Float_t dr_error = 0;   


/*

        for (Int_t ibin_eec = 1; ibin_eec <= eec_bins_tot; ibin_eec++) {

        eec = h->GetYaxis()->GetBinCenter(ibin_eec);
        if (ibin_eec == eec_bins_tot) eec = last_eec;

        double val = h->GetBinContent(ibin_dr, ibin_eec);
        h_1D->Fill(dr, val * eec);
        double err = h->GetBinError(ibin_dr, ibin_eec);
        if (!std::isnan(err) && !std::isinf(err))
            dr_error += pow(err * eec, 2);
    }

        */ 

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

