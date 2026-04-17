// needed headers
#include "TH2D.h"
#include "TAxis.h"


// Define binnings

// 
//dr
const Int_t dr_binsVectorSize = 16;
Int_t bins_dr = dr_binsVectorSize - 1;
Int_t dr_bins =  dr_binsVectorSize - 1;
// Equally sized bins in logarithmic scale... Zoe binning instead of Frac.
Double_t dr_binsVector[dr_binsVectorSize]= 
{
 0.002000000
, 0.002856000
, 0.004079000
, 0.005825000
, 0.008320000
, 0.011883000
, 0.016953000
, 0.024155000
, 0.034399000
, 0.049013000
, 0.069840000
, 0.099570000
, 0.141980000
, 0.202700000
, 0.289100000
, 0.400000000
};

Double_t dr_min = dr_binsVector[0];
Double_t dr_max = dr_binsVector[bins_dr];

//Values for histogram filling
// Franc. first bin was 0.15 to 0.4
Double_t dr_shiftbin = 0.0001; //  0.0001
Double_t dr_min_fill = 0.0021; //  0.016 = -0.15 + 0.01 
Double_t dr_max_fill = 0.39;// 0.39 = 0.4- 0.01

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



