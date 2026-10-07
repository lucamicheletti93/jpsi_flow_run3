void run_flow_fitter()
{
    gROOT->ProcessLine(".L FlowAnalysis_Helper.cxx+");
    gROOT->ProcessLine(".L FlowAnalysis_Fitting.cxx+");
    gROOT->ProcessLine(".L flow_fitter.C+");
    gROOT->ProcessLine("flow_fitter(10,0,4,1,0,0)");
}