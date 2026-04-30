// Organizational steps 

// --- Put all the template fit results in one directory, TemplateFitResult
// Under it: have the directories of the nominal abd varaitions, the root files that is input to to other drawings and calcuations.
// + the dir of FitResult_summary_S_B plots 
// + the directory of EEC_plots 

//--- Make Directory for Syetamtic uncertaintiy result only: Plots + its root file

//-- Change the name of the summary histos_templatefit_test.root (very important) --> To something means it is summary from Templatefit to raw EEC of signal. 



// --- How to run? 
Simply go to the path of file: template_fit.cpp 
and run it! 
it runs: 
template fit --> then draw Fit result S/B fractions --> Use them to Draw raw EEC (signal).
Then, it computes the syst. uncert. due to sticking the 0B to the total bkg template fit.

