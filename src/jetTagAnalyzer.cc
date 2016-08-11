//#include "jetTagAnalyzer.h"
#include "json/json.h"
#include "json/json-forwards.h"
//#include "jetSelector.h"
//#include "globalJetProbabilities.h"
#include "jetProbabilityMasterComputer.h"
#include <iostream>
#include <fstream>
#include <string.h>  
#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "rootlogon.C"

int main(int argc, char* argv[]) {

  defaultStyle style;
  style.setStyle();
  // names of the json files to be parsed
  std::string sample_json;
  std::string selection_json;
  std::string run_config_json;
  std::string globalProb_json;
  std::string ip_smear_json;
  // name of the json containing the global probabilities
  std::string globalProbabilities;
  // name of the tree to use in each sample
  std::string treeName;

  TH1D pileupWeightHist;
  TH1D pileupHist;
  TH1D mcPileupHist;

  // debug status
  int debug = 0;

  // indices corresponding to how the sample is divided
  int nDivisions = 1; // the number of total divisions of the sample
  int probIndex = -1; // the index of events (nEvn % nDivisions) for events used to derive the probaiblities
  int valiIndex = -1; // the index of events (nEvn % nDivisions) for events used to validate probabilities 
  bool runChop = false;
  std::string runChopSuffix = "";
  std::string runChopProbSuffix = "";

  // give default values for the configurations if not specified
  // when we want to use specific probabilities
  if (argc == 2) {
    sample_json	    = "sample.json";
    selection_json  = "selection.json";
    run_config_json = "run_config.json";
    ip_smear_json   = "ip_smear.json";
    globalProb_json = argv[1];
  }
  // when we want to do cross validation the probaiblities method
  else if(argc == 4) {
    sample_json	    = "sample.json";
    selection_json  = "selection.json";
    run_config_json = "run_config.json";
    ip_smear_json   = "ip_smear.json";
    globalProb_json = "none";

    // declare the chopping
    runChop	    = true;
    nDivisions	    = std::atoi(argv[1]);
    probIndex	    = std::atoi(argv[2]);
    valiIndex	    = std::atoi(argv[3]);

    runChopSuffix +=  "_"+std::to_string(nDivisions)+"div_";
    runChopSuffix += std::to_string(probIndex)+"prob_";
    runChopSuffix += std::to_string(valiIndex)+"vali";

    runChopProbSuffix +=  "_"+std::to_string(nDivisions)+"div_";
    runChopProbSuffix += std::to_string(probIndex)+"prob";

    std::cout << "-----------------------------" << std::endl;
    std::cout << " Running Cross Validation of Sample. " << std::endl;
    std::cout << " Training Sample for Probaiblities: " << probIndex << " of " << nDivisions << std::endl;
    std::cout << " Validiation Sample Index: " << valiIndex << " of " << nDivisions << std::endl;
    std::cout << "-----------------------------" << std::endl;

    if (probIndex < 0 || valiIndex < 0 || 
	probIndex > nDivisions - 1 ||
	valiIndex > nDivisions - 1 || 
	valiIndex == probIndex) {
      std::cout << "[WARNING] Possible Bad Indicies for Cross Validition" << std::endl;
      std::cout << "[WARNING]\t nDivisions" << nDivisions << std::endl;
      std::cout << "[WARNING]\t prob Index" << probIndex << std::endl;
      std::cout << "[WARNING]\t vali Index" << valiIndex << std::endl;
      //exit(1);      
    }
  }
  // when we want to do cross validation the probaiblities method
  else {
    sample_json	    = "sample.json";
    selection_json  = "selection.json";
    run_config_json = "run_config.json";
    ip_smear_json   = "ip_smear.json";
    globalProb_json = "none";
  }
  // else { //otherwise parse from the arguments
  //   sample_json	    = argv[1];
  //   selection_json  = argv[2];
  //   run_config_json = argv[3];
  //   globalProb_json = argv[4];
  // }

  // check if global probabilities were provided


  bool	probProvided	 = (globalProb_json != "none");

  if(debug > -1){
    std::cout << "sample      json name: " << sample_json << std::endl;
    std::cout << "selection   json name: " << selection_json << std::endl;
    std::cout << "run config  json name: " << run_config_json << std::endl;
    std::cout << "globalProb  json name: " << globalProb_json << std::endl;
  }

  // json to pull from
  Json::Value	sample_root;
  Json::Value	selection_root;
  Json::Value	run_config_root;
  Json::Value   globalProb_root;
  Json::Value   ip_smear_root;


  // file to input into the json reader
  std::ifstream sample_doc(sample_json, std::ifstream::binary);
  std::ifstream selection_doc(selection_json, std::ifstream::binary);
  std::ifstream run_config_doc(run_config_json, std::ifstream::binary);
  std::ifstream globalProb_doc(globalProb_json, std::ifstream::binary);
  std::ifstream ip_smear_doc(ip_smear_json, std::ifstream::binary);

  // pointer to the globalProbabilities
  globalJetProbabilities * globalJetProbToApply = NULL; 

  // input the files into the JSON values
  if(debug > -1) std::cout << "Reading in JSONs..." << std::endl;
  if(debug > -1) std::cout << "\t..." << sample_json << std::endl;
  sample_doc     >> sample_root;
  if(debug > -1) std::cout << "\t..." << selection_json << std::endl;
  selection_doc  >> selection_root;
  if(debug > -1) std::cout << "\t..." << run_config_json << std::endl;
  run_config_doc >> run_config_root;
  if(debug > -1) std::cout << "\t..." << ip_smear_json << std::endl;
  ip_smear_doc >> ip_smear_root;

  // read parameters from the run config
  if(debug > -1) std::cout << "Building Output File..." << std::endl;
  std::string	outputFileName		    = run_config_root.get("outputFileName","test").asString();
  std::string	pileupFileName		    = run_config_root.get("puProfile","").asString();
  std::string	mcPileupFileName	    = run_config_root.get("mcPuProfile","").asString();
  bool	        applyPuWeight		    = run_config_root.get("applyPUWeight",false).asBool();
  std::string	outputDir		    = run_config_root.get("outputDir","/tmp/hardenbr/").asString();
  debug					    = run_config_root.get("debug",0).asInt();
  long int	beginEvent		    = run_config_root.get("beginEvent",0).asInt();
  long int	endEvent		    = run_config_root.get("endEvent",-1).asInt();
  //limit on the number of ntags to compute
  int		maxJetTags		    = run_config_root.get("maxJetTags",3).asInt(); 
  bool		runSignalContam		    = run_config_root["signalContam"].get("run",false).asBool(); 
  bool		writeTree		    = run_config_root.get("writeTree",false).asBool(); 
  bool		onlyProbs		    = run_config_root.get("onlyProbs",false).asBool(); 
  bool		removeSignalRegionFromProbs = run_config_root.get("removeSignalRegionFromProbs",false).asBool(); 

  std::string	hypotheticalChopProb = outputDir+"/"+"prob" + runChopProbSuffix + ".json";
  if(runChop) {
    std::cout << "Checking for already calculated probabilities for cross validation: " << hypotheticalChopProb << std::endl;
  }

  std::ifstream hypotheticalChopProb_doc(hypotheticalChopProb.c_str());
  bool		chopProbExists	     = hypotheticalChopProb_doc.good() && runChop;

  if(debug > -1) std::cout << "\t..." << outputFileName << std::endl;
  // read in the json if the config is provided 

  if(debug > -1) std::cout << "Checking for provided global probabilities configuration..." << std::endl;

  if(probProvided) {
    globalProb_doc >> globalProb_root;
    std::cout << "Global Jet Probabilities provided. Parsing JSON... " << std::endl;
  }
  else if(chopProbExists) {
    hypotheticalChopProb_doc >> globalProb_root;
    std::cout << "Probabilities have already been calculated for cross validation, loading.... " << std::endl;
    probProvided = true;
  }
  else {
    std::cout << "No global jet probabilities provided. Will compute per sample... " << std::endl;
  }

  std::cout << "Run Signal Contamination..." << runSignalContam << std::endl;

  if(debug > -1) std::cout << "Building Jet Selection..." << std::endl;

  // if we are doing the cross validation we need to modify the selection json

  // build the jet selector
  jetSelector jetSel(selection_root, runChop, nDivisions, probIndex, valiIndex, debug);

  // parse the name of the tree to use for the globalProbabilities
  treeName = selection_root.get("tree","jets").asString();  

  // access the sample parameterization
  const Json::Value samples = sample_root["samples"];
  std::cout << "Number of samples: " << samples.size() << std::endl;

  // Get the eventweight for the pile-up. This only needs to be be done once
  float puWeight = 1;
  if(pileupFileName != "") {
    if(debug > -1) std::cout << "[PU WEIGHTING] Getting Pileup Distribution....." << std::endl;
    TFile puFile(pileupFileName.c_str(),"READ");
    pileupHist  = *(TH1D*)puFile.Get("pileup");
    pileupHist.Scale(1.0 / pileupHist.Integral());
    std::cout << "[DATA PU PROFILE] " << std::endl;
    pileupHist.Print();

    TFile mcPuFile(mcPileupFileName.c_str(),"READ");
    mcPileupHist = *(TH1D*)mcPuFile.Get("pileup");
    mcPileupHist.Scale(1.0 / mcPileupHist.Integral());
    std::cout << "[MC PU PROFILE] " << std::endl;
    mcPileupHist.Print();

    pileupWeightHist = *(TH1D*)pileupHist.Clone("weights");

    TH1F clone("pileup_Rebin","pileup_rebin",50,0,50);
    //    clone.Add(&mcPileupHist);

    for(int ii=0; ii < 50; ++ii) {
      int bin = mcPileupHist.GetBin(ii);
      clone.SetBinContent(bin, mcPileupHist.GetBinContent(bin));
    }

    pileupWeightHist.Divide(&clone);
     
    if(debug > 3) std::cout << "[PU WEIGHTING] Finished Pileup Distribution....." << std::endl;
  }

  // output file containing combined information of each sample
  std::string rootFileName = outputDir+"/"+outputFileName;
  if(runChop) {
    rootFileName+=runChopSuffix;
  }
  rootFileName += ".root";

  TFile outputFile((rootFileName).c_str(), "RECREATE");
  outputFile.cd();

  if(debug > -1) std::cout << "Begin Sample Loop..." << std::endl;
  // loop over each sample and build the global jet probabilities 
  for(int ss = 0; ss < samples.size(); ++ss ) {
    bool runSample  = samples[ss].get("runSample", false).asBool();
    if(!runSample) continue;

    // extract sample information from json
    std::string label = samples[ss].get("label", "NO_LABEL").asString();
    std::string path  = samples[ss].get("path", "NO_PATH_PROVIDED").asString();
    std::string stack = samples[ss].get("stack", "NO_STACK_PROVIDED").asString();
    bool	isMC  = samples[ss].get("isMC", true).asBool();
    bool	isSig = samples[ss].get("isSig", false).asBool();
    double	xsec  = samples[ss].get("xsec", 1).asFloat();
    float       x_limit_label = 0;
    float       y_limit_label = 0;

    // running pdf totals
    float	    pdfWeightTotal   = 0;
    float	    pdfWeightTotalUp = 0;
    float	    pdfWeightTotalDn = 0;

    // parse the values to be used for the limit calculation (ONLY IF SIGNAL)
    if(isSig) {
      x_limit_label = samples[ss].get("x_limit_label", -999).asFloat();      
      y_limit_label = samples[ss].get("y_limit_label", -999).asFloat();      
    }

    // build  names for tagging histograms  based on the label
    std::string nTagHistTrueName	 = label + "_nTagTrue";
    std::string nTagHistPredName	 = label + "_nTagPred";
    std::string nTagHistTruePDFName  	 = label + "_nTagTruePDF";
    std::string nTagHistTrueHLTName	 = label + "_nTagTrueHLT";
    std::string nTagHistTrueTrackingName = label + "_nTagTrueTracking";
    

    // initialize histograms
    TH1D nTagHistTrue(nTagHistTrueName.c_str(), "N Tags Obs.", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistPred(nTagHistPredName.c_str(), "N Tags Exp.", maxJetTags, 1-.5, maxJetTags+1-.5);
    // pdf RMS Weight applied
    TH1D nTagHistTruePDFUp((nTagHistTruePDFName+"Up").c_str(), "N Tags PDF Up", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistTruePDFDn((nTagHistTruePDFName+"Dn").c_str(), "N Tags PDF Dn", maxJetTags, 1-.5, maxJetTags+1-.5);
    // prediction fake rate statistical error propogated 
    TH1D nTagHistPredErrUp((nTagHistPredName+"ErrUp").c_str(), 
			   "Fake Rate Stat.", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistPredErrDn((nTagHistPredName+"ErrDn").c_str(), 
			   "Fake Rate Stat.", maxJetTags, 1-.5, maxJetTags+1-.5);
    // Systematics for the hlt simulation
    // up
    TH1D nTagHistTrueHLT((nTagHistTrueHLTName+"Nominal").c_str(), 
			   "HLT Simul..", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistTrueHLTSYSUp((nTagHistTrueHLTName+"ErrUp").c_str(), 
			   "HLT Simul..", maxJetTags, 1-.5, maxJetTags+1-.5);
    // dn
    TH1D nTagHistTrueHLTSYSDn((nTagHistTrueHLTName+"ErrDn").c_str(), 
			   "HLT Simul.", maxJetTags, 1-.5, maxJetTags+1-.5);
    // Systematics for the offline tracking variable modeling
    // up
    TH1D nTagHistTrueTrackingSYSUp((nTagHistTrueTrackingName+"ErrUp").c_str(), 
			   "Tag Var Modeling", maxJetTags, 1-.5, maxJetTags+1-.5);
    // dn
    TH1D nTagHistTrueTrackingSYSDn((nTagHistTrueTrackingName+"ErrDn").c_str(), 
				   "Tag Var Modeling.", maxJetTags, 1-.5, maxJetTags+1-.5);
    
    if(debug > -1) {
      std::cout << "label: " << label << std::endl;
      std::cout << "\t path: " << path << std::endl;
      std::cout << "\t stack: " << stack << std::endl;
      std::cout << "\t isMC: " << isMC << std::endl;
      std::cout << "\t isSig: " << isSig << std::endl;
      std::cout << "\t xsec: " << xsec << std::endl;
    }

    ///////
    // STEP 1: GENERATE THE JET PROBABILITIES IF NECESSARY
    ///////

    if(debug > 3) std::cout << "Opening file to analyze....." << std::endl;
    TFile thisFile(path.c_str(),"READ");
    thisFile.cd();
    if(debug > 3) std::cout << "Getting the tree of the file....." << std::endl;
    TTree * tree       = (TTree*)(thisFile.Get(treeName.c_str()));    
    TTree * caloHTTree = (TTree*)(thisFile.Get(treeName.c_str()));    

    std::string contamPath;

    double  total_processed	  = tree->GetEntries();
    double  filtered         	  = tree->GetEntries();

    // get the eventWeight related to the processed piece of the sample
    if(isMC && !isSig) {
      std::cout << "Getting the filter count tree ....." << std::endl;
      TTree * filterTree = (TTree*)(thisFile.Get("filtCount"));    
      filterTree->Draw("nEventsTotal>>event_total_hist","nEventsTotal","goff");
      filterTree->Draw("nEventsFiltered>>event_filt_hist","nEventsFiltered","goff");
      TH1F* totalEventHist    = (TH1F*)gDirectory->Get("event_total_hist");
      TH1F* filteredEventHist = (TH1F*)gDirectory->Get("event_filt_hist");
      total_processed	      = totalEventHist->Integral();
      filtered		      = filteredEventHist->Integral();
      std::cout << "Total MC events processed: " << total_processed <<  std::endl;
      std::cout << "Total MC events post-filter: " << filtered <<  std::endl;
      std::cout << "Corresponding Event Weight:: " << xsec / total_processed <<  std::endl;
    }    

    // First get some global information about the sample
    int	    nPassEventSelection	  = 0;
    int	    nPassTriggerSelection = 0;
    double  eventWeight		  = xsec / total_processed;
    double  totalWeight = 1;    
    // is globalProbabilities were not computed, produce them for each sample
    if(!probProvided) {
      if(debug > -1) std::cout << "Probabilities not provided.." << std::endl;
      if(debug > -1) std::cout << "Constructing localProbabilities..." << std::endl;
      // process the sample by constructing the global probabilities
      globalJetProbabilities * localSampleProb = new globalJetProbabilities(label, stack, isMC, isSig, eventWeight, xsec, tree, jetSel, debug);     

      // subtract the signal region
      if(removeSignalRegionFromProbs) localSampleProb->removeSignalRegion(tree, jetSel, false, 1);
      
      if(runSignalContam) {

	// parse the parameters
	float	    norm      = run_config_root["signalContam"].get("norm",0).asFloat();
	std::string contLabel = run_config_root["signalContam"].get("label","NO_LABEL").asString();

	std::cout << "Running Signal Contamination with sample : "  << contLabel << " using normalization " << norm << std::endl;

	// find the sample
	bool foundContam = false;
	for(int sss = 0; sss < samples.size(); ++sss ) {	    
	  // extract sample informatino from json
	  std::string label = samples[sss].get("label", "NO_LABEL").asString();	  
	  if (label != contLabel) continue;
	  foundContam = true;	  

	  // parse the path
	  contamPath  = samples[sss].get("path", "NO_PATH_PROVIDED").asString();
	  
	  TFile   contamFile(contamPath.c_str(), "READ");
	  TTree * contamTree = (TTree*)contamFile.Get(treeName.c_str());
	  
	  std::cout << "\n\n Adding signal contamination to probabilities from path: " << contamPath.c_str() << std::endl;
	  // add the signal contamination to the local probabilities
	  localSampleProb->addSignalContamination(contamTree, jetSel, norm);	  

	  // and then remove the piece contained in the signal region
	  std::cout << "\n\n removing signal region from probabilities " << std::endl;
	  if(removeSignalRegionFromProbs) localSampleProb->removeSignalRegion(contamTree, jetSel, true, norm);
	} // end loop over signal samples for contamination	
	
	// make sure we found
	if(!foundContam) {
	  std::cout << "[ERROR] did not find sample for signal contamaination....exiting" << std::endl;
	  exit(1);
	}	
	
      } // end additional signal contamination
      
      // parse the json assembled
      std::cout << "\n\n getting sample json to output probabilities" << std::endl;
      Json::Value  sampleJson = localSampleProb->getProbabilitiesJSON();

      // write the json to the stream in the output directory
      std::string jsonOutputName = outputDir + "/prob";

      if(runChop) { 
	jsonOutputName += runChopProbSuffix;
      }
      else {
	jsonOutputName +="_"+label;
      }
      jsonOutputName += ".json";

      std::cout << "\n\n Writing JSON output probabilities to file: " << jsonOutputName << std::endl;
      std::ofstream json_stream;
      json_stream.open(jsonOutputName);      
      // use the style writter
      Json::StyledWriter styledWriter;
      json_stream << styledWriter.write(sampleJson);
      json_stream.close();

      // write the histograms for the file
      if(debug > -1) std::cout << "Writing probability hists..." << std::endl;
      outputFile.cd();      
      localSampleProb->getTaggedHist().Write();
      localSampleProb->getAllHist().Write();
      localSampleProb->getRatioGraph().Write();
      localSampleProb->getCentralEffHist().Write();
      localSampleProb->getCentralEffHistErrUp().Write();
      localSampleProb->getCentralEffHistErrDn().Write();
      pileupWeightHist.Write();
      // set the probabilities to apply to this samples local probabilities
      globalJetProbToApply = localSampleProb; 
    } // end building global probabilities
    // probaiblities were provided, or previously calculated
    else {
      std::cout << "....Loading JSON Probabilities......." << std::endl;
      globalJetProbabilities * loadedProbabilities = new globalJetProbabilities(label, stack, isMC, isSig, eventWeight, xsec, globalProb_root, debug);
      globalJetProbToApply = loadedProbabilities;

      outputFile.cd();
      // write the efficiency histograms
      loadedProbabilities->getCentralEffHist().Write();
      loadedProbabilities->getCentralEffHistErrUp().Write();
      loadedProbabilities->getCentralEffHistErrDn().Write();
      // write the histograms going into the efficiency
      loadedProbabilities->getAllHist().Write();
      loadedProbabilities->getTaggedHist().Write();      

    } // end loading global probaiblities from json


    // check the status of all objects in the global probabilities to by applied
    globalJetProbToApply->printHistStatus();

    ///////////////////
    // STEP 2: APPLY THE GLOBAL / LOCAL PROBABILITIES TO THE SAMPLE 
    // EVENT BY EVENT
    ///////////////////

    if(debug > -1) std::cout << "------ BEGINNING STEP TWO ----- " << std::endl;
    // make one file per sample for the individual tree

    std::string outputFileName = outputDir + "/tree_" + label;
    if(runChop) outputFileName += "_" + runChopSuffix;
    outputFileName += ".root";
    if(debug > -1) std::cout << " Building sample output file:  " << outputFileName << std::endl;
    TFile sampleOutputFile(outputFileName.c_str(), "RECREATE");    

    // create the tree that will contain the variables from
    // a shallow copy of the original tree
    TTree* jetVariableTree = jetSel.shallowCopyTree(tree);

    // arrays related to the probability tagging    
    int nCat		      = 11;
    int nTagWeight[11]	      = {0,0,0,0,0,0,0,0,0,0};
    double probNTags[11]      = {0,0,0,0,0,0,0,0,0,0};
    double probNTagsErrUp[11] = {0,0,0,0,0,0,0,0,0,0};
    double probNTagsErrDn[11] = {0,0,0,0,0,0,0,0,0,0};
    int   nJetTagArray[11]    = {0,1,2,3,4,5,6,7,8,9,10};
    // displaced
    int nRegDispTracks[11]	 = {0,0,0,0,0,0,0,0,0,0};
    int nRegDispTracksUp[11]     = {0,0,0,0,0,0,0,0,0,0};
    int nRegDispTracksDn[11]     = {0,0,0,0,0,0,0,0,0,0};
    // prompt tracks
    int nRegPromTracks[11]      = {0,0,0,0,0,0,0,0,0,0};
    int nRegPromTracksUp[11]     = {0,0,0,0,0,0,0,0,0,0};
    int nRegPromTracksDn[11]     = {0,0,0,0,0,0,0,0,0,0};
    // HLT status
    // inclusive
    int isInclusive[11];
    int isInclusiveUp[11];
    int isInclusiveDn[11];
    // displaced track
    int isDisplacedTrack[11];
    int isDisplacedTrackUp[11];
    int isDisplacedTrackDn[11];
   
    // jet indexed
    double probJTag[100];
    int   isTagged[100];
    int   isTaggedUp[100];
    int   isTaggedDn[100];

    // event indexed / flat number
    int		nTagged		      = 0;
    int		nTaggedUp	      = 0;
    int		nTaggedDn	      = 0;
    long int	evNum		      = 0;
    float       caloHT		      = 0;
    float       subLeadingJetPt	      = 0;
    float       leadingJetPt	      = 0;
    // ONLINE TRACKING SYSTEMATICS FOR SIGNAL
    // inclusive canddiates
    int		nInclusiveJets	      = 0;
    int		nInclusiveJetsUp      = 0;
    int		nInclusiveJetsDn      = 0;
    // displaced track path candidates
    int		nDisplacedTrackJets   = 0;
    int		nDisplacedTrackJetsUp = 0;
    int		nDisplacedTrackJetsDn = 0;


    // set the branches for the probabilities
    jetVariableTree->Branch("evNum", &evNum, "evNum/I");

    jetVariableTree->Branch("nCat", &nCat, "nCat/I");
    jetVariableTree->Branch("index", &nJetTagArray, "index[nCat]/I");
    jetVariableTree->Branch("nTagWeight", &nTagWeight, "nTagWeight[nCat]/I");
    jetVariableTree->Branch("probNTags", &probNTags, "probNTags[nCat]/D");
    jetVariableTree->Branch("probNTagsErrUp", &probNTagsErrUp, "probNTagsErrUp[nCat]/D");
    jetVariableTree->Branch("probNTagsErrDn", &probNTagsErrDn, "probNTagsErrDn[nCat]/D");
    // 
    // jet indexed
    //
    jetVariableTree->Branch("probJTag", &probJTag, "probJTag[nCaloJets]/D");    
    jetVariableTree->Branch("isTagged", &isTagged, "isTagged[nCaloJets]/I");    
    jetVariableTree->Branch("isTaggedUp", &isTaggedUp, "isTaggedUp[nCaloJets]/I");    
    jetVariableTree->Branch("isTaggedDn", &isTaggedDn, "isTaggedDn[nCaloJets]/I");    
    // inclusive
    jetVariableTree->Branch("isInclusive", &isInclusive, "isInclusive[nCaloJets]/I");    
    jetVariableTree->Branch("isInclusiveUp", &isInclusiveUp, "isInclusiveUp[nCaloJets]/I");    
    jetVariableTree->Branch("isInclusiveDn", &isInclusiveDn, "isInclusiveDn[nCaloJets]/I");    
    // displaced
    jetVariableTree->Branch("isDisplacedTrack", &isDisplacedTrack, "isDisplacedTrack[nCaloJets]/I");    
    jetVariableTree->Branch("isDisplacedTrackUp", &isDisplacedTrackUp, "isDisplacedTrackUp[nCaloJets]/I");    
    jetVariableTree->Branch("isDisplacedTrackDn", &isDisplacedTrackDn, "isDisplacedTrackDn[nCaloJets]/I");    

    // track counts generated from the jet PDFS
    jetVariableTree->Branch("nRegPromTracks", &nRegPromTracks, "nRegPromTracks[nCaloJets]/I");    
    jetVariableTree->Branch("nRegPromTracksUp", &nRegPromTracksUp, "nRegPromTracksUp[nCaloJets]/I");    
    jetVariableTree->Branch("nRegPromTracksDn", &nRegPromTracksDn, "nRegPromTracksDn[nCaloJets]/I");    
    jetVariableTree->Branch("nRegDispTracks", &nRegDispTracks, "nRegDispTracks[nCaloJets]/I");    
    jetVariableTree->Branch("nRegDispTracksUp", &nRegDispTracksUp, "nRegDispTracksUp[nCaloJets]/I");    
    jetVariableTree->Branch("nRegDispTracksDn", &nRegDispTracksDn, "nRegDispTracksDn[nCaloJets]/I");        

    // 
    // event indexed / flat number
    //
    jetVariableTree->Branch("nTagged", &nTagged, "nTagged/I");    
    jetVariableTree->Branch("nTaggedUp", &nTaggedUp, "nTaggedUp/I");    
    jetVariableTree->Branch("nTaggedDn", &nTaggedDn, "nTaggedDn/I");    
    // online ip requirements
    // inclusive
    jetVariableTree->Branch("nInclusiveJets", &nInclusiveJets, "nInclusiveJets/I");    
    jetVariableTree->Branch("nInclusiveJetsUp", &nInclusiveJetsUp, "nInclusiveJetsUp/I");    
    jetVariableTree->Branch("nInclusiveJetsDn", &nInclusiveJetsDn, "nInclusiveJetsDn/I");    
    // displaced
    jetVariableTree->Branch("nDisplacedTrackJets", &nDisplacedTrackJets, "nDisplacedTrackJets/I");    
    jetVariableTree->Branch("nDisplacedTrackJetsUp", &nDisplacedTrackJetsUp, "nDisplacedTrackJetsUp/I");    
    jetVariableTree->Branch("nDisplacedTrackJetsDn", &nDisplacedTrackJetsDn, "nDisplacedTrackJetsDn/I");    
    // event level variables
    jetVariableTree->Branch("caloHT", &caloHT, "caloHT/F");
    jetVariableTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
    jetVariableTree->Branch("subLeadingJetPt", &subLeadingJetPt, "subLeadingJetPt/F");
        
    if(debug > -1) std::cout << "Initializing master probability calculator  " << std::endl;
    // construct the master jet probability calculator to do the n tag calculations
    if(globalJetProbToApply == NULL) std::cout << "global Jet Probabilities not chosen....NULL POINTER" << std::endl;
    jetProbabilityMasterComputer masterJetCalc(globalJetProbToApply, tree, debug);

    // loop event by event
    if(debug > -1) std::cout << " Beginning Event Loop  " << std::endl;
    int nEvents = tree->GetEntries();
    // fill the number of events passing the trigger

    nPassTriggerSelection = tree->GetEntries((globalJetProbToApply->triggerCutOnlyString).c_str());

    std::cout << "\n TOTAL NUMBER OF EVENTS TO PROCESS: " << nEvents << std::endl;

    // add variables for determining the event rate
    jetVariableTree->Branch("eventWeight", &eventWeight, "eventWeight/D");
    jetVariableTree->Branch("puWeight", &puWeight, "puWeight/F");
    jetVariableTree->Branch("xsec", &xsec, "xsec/D");
    jetVariableTree->Branch("eventsAnal", &total_processed, "eventsAnal/D");
    jetVariableTree->Branch("eventsInFile", &nEvents, "eventsInFile/I");
    jetVariableTree->Branch("eventsPassingSel", &nPassEventSelection, "eventsPassingSel/I");      
    jetVariableTree->Branch("eventsPassingTrigger", &nPassTriggerSelection, "eventsPassingTrigger/I");      

    if(endEvent > 0 && endEvent+1 < nEvents) nEvents = endEvent+1;
    if(onlyProbs) nEvents = 0;
    if(debug > 5) std::cout << "Pre-Event Loop beginEvent: " << beginEvent << " endEvent " << endEvent << std::endl;

    for(long int event = 0; event < nEvents; ++event) {
      // start at the first event designated
      if(debug > 5) std::cout << "begin event loop -- event #" << event << " " << std::endl;
      if(event < beginEvent) continue; 
      tree->GetEntry(event);
      if(event % 20000 == 0) std::cout << "Processing Event # --- "  << event << std::endl;
      if(event % 1000 == 0 && isSig) std::cout << "Processing Signal Event # --- "  << event << std::endl;
      // check the index matches for the validation sample and that the kinematic / trigger requirements are satisfied
      int   evNum		= tree->GetLeaf("evNum")->GetValue(0);
      bool  passValidationIndex	= (evNum % nDivisions == valiIndex) || !runChop;    // is the correct validation sample or no chop

      // speed up when we are running pulls
      // no need to calculate the number of tags if we are performing a chop with nDiv > 2
      // because we are not subtracting the signal region from the probabilities
      if(nDivisions > 2 && runChop && !passValidationIndex) continue;

      // Calculate the number of tagged jets
      if(debug > 2) std::cout << "Getting the nTagged Jets Vector.." << std::endl;
      const std::vector<bool> taggedVector   = jetSel.getJetTaggedVector(tree, event, false, false);      
      std::vector<bool> taggedVectorUp;// = taggedVector;
      std::vector<bool> taggedVectorDn;// = taggedVector;
      
      // only run the smearing if we are running signal
      if(isSig) {
	std::vector<bool> taggedVectorUp_temp =  jetSel.getJetTaggedVector(tree, event, true, true);
	std::vector<bool> taggedVectorDn_temp =  jetSel.getJetTaggedVector(tree, event, true, false);

	for(int vvv = 0; vvv < int(taggedVector.size()); ++vvv) {
	  bool valUp =  taggedVectorUp_temp[vvv];
	  bool valDn =  taggedVectorDn_temp[vvv];
	  taggedVectorUp.push_back(valUp);
	  taggedVectorDn.push_back(valDn);
	}

      }
      else {
	for(int vvv = 0; vvv < int(taggedVector.size()); ++vvv) {
	  bool val =  taggedVector[vvv];
	  taggedVectorUp.push_back(val);
	  taggedVectorDn.push_back(val);
	}
      }

      // variable for tracking the number of tags in the event
      // track systematically 
      nTagged	= 0;		// the total number of jets tagged per event      
      nTaggedUp = 0;
      nTaggedDn = 0;

      if(debug > 2) std::cout << "Filling output Tree Branch.. with tagvector size:" << taggedVector.size() << std::endl;
      if(debug > 5) std::cout << "[jetTagAnalyzer] vector of tags";
      for(int jj = 0; jj < taggedVector.size(); ++jj) {
	if(debug > 5) std::cout << taggedVector[jj];	
	// check for each vector being tagged
	if(taggedVector[jj]) nTagged++; 
	if(taggedVectorUp[jj]) nTaggedUp++;
	if(taggedVectorDn[jj]) nTaggedDn++;

	isTagged[jj]   = taggedVector[jj] ? 1 : 0;
	isTaggedUp[jj] = taggedVectorUp[jj] ? 1 : 0;
	isTaggedDn[jj] = taggedVectorDn[jj] ? 1 : 0;
      }
      
      // is in the signal region for <= 2 divisions
      bool  passSignalTagRegion		= (nTagged >= 2 && nDivisions <= 2);  
      // passes the trigger and kinematic requirements
      bool  eventPassKinematicSelection = jetSel.doesEventPassSelection(tree, event);
      // two ways to get into the validation sample (signal region or correct validation index)
      bool  eventPassSelection		= eventPassKinematicSelection && (passValidationIndex || passSignalTagRegion); 
 
      // when we use 1 or 2 divisions, keep the full signal region       
      if(eventPassSelection) nPassEventSelection++;	

      if(debug > 2) std::cout << "Event passes event selection?? "  << eventPassSelection << std::endl;
      evNum = event;

      if(debug > 9) std::cout << "\t\t\tProcessing Event # --- "  << event << std::endl;
      // get the calo HT calculation 
      caloHTTree->GetEntry(event);
      if(debug > 9) std::cout << "Retrieved Entry" << std::endl;
      TLeaf *   caloHTLeaf	    = caloHTTree->GetLeaf("eventCaloHT");
      int       nHT		    = caloHTLeaf->GetNdata();
      TLeaf *   leadingJetPtLeaf    = caloHTTree->GetLeaf("eventCaloHT");
      int       nPt		    = leadingJetPtLeaf->GetNdata();
      TLeaf *   subLeadingJetPtLeaf = caloHTTree->GetLeaf("eventCaloHT");
      int       nSubPt		    = subLeadingJetPtLeaf->GetNdata();
      if(debug > 9) std::cout << "Retrieved Leafs..." << std::endl;
      if(nHT != 1) {
	std::cout << "caloHT leaf does not contain a single entry" << std::endl;	
	exit(1);
      }
      else {
	caloHT		= caloHTLeaf->GetValue(0);
	leadingJetPt	= 0; //leadingJetPtLeaf->GetValue(0);
	subLeadingJetPt = 0; //subLeadingJetPtLeaf->GetValue(0);
      }

      // get the puWeight
      if(isMC && pileupFileName != "") {
	if(debug > 9) std::cout << "Retrieving PU INFO..." << std::endl;
	TLeaf *	genPULeaf = tree->GetLeaf("genPU");
	float	truePU	  = genPULeaf->GetValue(0);           
	puWeight	  = pileupWeightHist.GetBinContent(pileupWeightHist.GetBin(truePU));
	if(applyPuWeight) totalWeight *= puWeight;
	if(debug > 0) std::cout << "Found PU weight..." << puWeight << std::endl;
      }      

      // PDF SMEARING OF ONLINE TRACKING FOR AOD SIGNAL
      // only build online tracking quantities for signal
      // 
      if(isSig) {

	if (debug > 2) std::cout << "building up online tracking from JSON" << std::endl;
	// reset all the counters
	nInclusiveJets	      = 0;
	nInclusiveJetsUp      = 0;
	nInclusiveJetsDn      = 0;
	// displaced track path candidates
	nDisplacedTrackJets   = 0;
	nDisplacedTrackJetsUp = 0;
	nDisplacedTrackJetsDn = 0;	

	if (debug > 2) std::cout << "throwing toys" << std::endl;
	// throw the toys using the json pdf's
	std::vector<std::tuple<int,int,int> > promOnline = 
	  jetSel.buildOnlineTrackingFromJSON(tree, event, int(y_limit_label), ip_smear_root, true);
	std::vector<std::tuple<int,int,int> > dispOnline = 
	  jetSel.buildOnlineTrackingFromJSON(tree, event, int(y_limit_label), ip_smear_root, false);      

	// // keep track of the candidate status
	// std::vector<bool> isPromVec, isPromUpVec, isPromDnVec;
	// std::vector<bool> isDispVec, isDispUpVec, isDispDnVec;
	// std::vector<bool> isPromAndDispVec;//#(50,false);
	// std::vector<bool> isPromAndDispUpVec;//(50,false);
	// std::vector<bool> isPromAndDispDnVec;//(50,false);
	
	int nJets =  tree->GetLeaf("nCaloJets")->GetValue(0);
	if (nJets > 13) {
	  //std::cout << "TOO MANY JETS CONTINUE" << std::endl;
	  continue; 
	}
	if (debug > 2) std::cout << "loop over jets for online tracking" << std::endl;	
	// check the smearing is consitent in size
	if(promOnline.size() != dispOnline.size()) {
	  std::cout << "pdf JSONS do not match in size" << std::endl;
	}
	else { // determine if the jets are prompt or displaced
	  for(int pdfJet = 0; pdfJet < nJets; ++pdfJet) {
	    if (debug > 2) std::cout << "prase results from jet selector" << std::endl;	
	    // parse the results from the tuples calculated by the jet selector
	    // prom
	    int nProm	= std::get<0>(promOnline[pdfJet]);
	    int nPromUp = std::get<1>(promOnline[pdfJet]);
	    int nPromDn = std::get<2>(promOnline[pdfJet]);
	    // disp
	    int nDisp	= std::get<0>(dispOnline[pdfJet]);
	    int nDispUp = std::get<1>(dispOnline[pdfJet]);
	    int nDispDn = std::get<2>(dispOnline[pdfJet]);	    

	    // prom track multiplicities
	    bool    isProm	    = nProm <= 2; 
	    bool    isPromUp	    = nPromUp <= 2; 
	    bool    isPromDn	    = nPromDn <= 2; 
	    // disp track multiplicities
	    bool    isDisp	    = nDisp >= 1; 
	    bool    isDispUp	    = nDispUp >=1; 
	    bool    isDispDn	    = nDispDn >=1; 
	    // prom and disp track multiplicities
	    bool    isPromAndDisp   = nDisp && nProm;
	    bool    isPromAndDispUp = nDispUp && nPromUp;
	    bool    isPromAndDispDn = nDispDn && nPromDn;

	    if (debug > 2) std::cout << "fill jet index candidates" << std::endl;	
	    // fill the jet index variables by category
	    isInclusive[pdfJet]	       = isProm ? 1 : 0;
	    isInclusiveUp[pdfJet]      = isPromUp ? 1 : 0;
	    isInclusiveDn[pdfJet]      = isPromDn ? 1 : 0;
	    isDisplacedTrack[pdfJet]   = isPromAndDisp ? 1 : 0;
	    isDisplacedTrackUp[pdfJet] = isPromAndDispUp ? 1 : 0;
	    isDisplacedTrackDn[pdfJet] = isPromAndDispDn ? 1 : 0;

	    if (debug > 2) std::cout << "set regional track quantities" << std::endl;	
	    // set the regional quantities
	    nRegPromTracks[pdfJet]   = nProm;
	    nRegPromTracksUp[pdfJet] = nPromUp;
	    nRegPromTracksDn[pdfJet] = nPromDn;
	    nRegDispTracks[pdfJet]   = nDisp;
	    nRegDispTracksUp[pdfJet] = nDispUp;
	    nRegDispTracksDn[pdfJet] = nDispDn;

	    // // for the inclusive trigger candidates
	    // isPromVec.push_back(isProm);
	    // isPromUpVec.push_back(isPromUp);
	    // isPromDnVec.push_back(isPromDn);
	    // // not used, except for the displace dtrack trigger candidates
	    // isDispVec.push_back(isDisp);
	    // isDispUpVec.push_back(isDispUp);
	    // isDispDnVec.push_back(isDispDn);
	    // // for the displaced track trigger
	    // isPromAndDispVec.push_back(isPromAndDisp);
	    // if (debug > 3) std::cout << " isPromAndDispVec size " << isPromAndDispVec.size() << " pushing back " << isPromAndDisp << std::endl; 
	    // isPromAndDispUpVec.push_back(isPromAndDispUp);
	    // isPromAndDispDnVec.push_back(isPromAndDispDn);	    

	    // isPromAndDispVec.resize(nJets);
	    // isPromAndDispUpVec.resize(nJets);
	    // isPromAndDispDnVec.resize(nJets);

	    if (debug > 2) std::cout << "count the candidates" << std::endl;	
	    // inclusive jets
	    if(isProm) nInclusiveJets++;
	    if(isPromUp) nInclusiveJetsUp++;
	    if(isPromDn) nInclusiveJetsDn++;
	    // displaced and prompt
	    if(isPromAndDisp) nDisplacedTrackJets++;
	    if(isPromAndDispUp) nDisplacedTrackJetsUp++;
	    if(isPromAndDispDn) nDisplacedTrackJetsDn++;	    
	  }
	  if (debug > 2) std::cout << "leaving jet loop scope" << std::endl;	
	}
	if (debug > 2) std::cout << "leaving else scope" << std::endl;	
      }
      if (debug > 2) std::cout << "leaving isSig scope" << std::endl;	
	// if (debug > 2) std::cout << "building candidates based on pseduo online tracking" << std::endl;
	// // determine how many prompt and displaced (HLT) jets there are in an event
	// for(int cand = 0; cand < nJets; ++cand) {
	//   // prompt
	//   if (debug > 2) std::cout << "inclusive candidates " << std::endl;
	//   if(isPromVec[cand]) nInclusiveJets++;
	//   if(isPromUpVec[cand]) nInclusiveJetsUp++;
	//   if(isPromDnVec[cand]) nInclusiveJetsDn++;
	//   // displaced and prompt
	//   if (debug > 2) std::cout << "displacedTrack candidates " << std::endl;
	//   if (debug > 2) std::cout << " isPromAndDispVec size" << isPromAndDispVec.size() << 
	// 		   "isPromAndDispVecUp" << isPromAndDispUpVec.size() <<  "isPromAndDispVecDn" << isPromAndDispDnVec.size() << std::endl;
	//   if (debug > 2) std::cout << "counting  candidates " << std::endl;
	//   if(isPromAndDispVec[cand]) nDisplacedTrackJets++;
	//   if(isPromAndDispUpVec[cand]) nDisplacedTrackJetsUp++;
	//   if(isPromAndDispDnVec[cand]) nDisplacedTrackJetsDn++;
	//   if (debug > 2) std::cout << "candidates complete " << std::endl;
	//} // end loop over vectors of jet canddiates
	//} // end if statement for doing tracking PDF smearing

      if (debug > 2) std::cout << "Checking probability vector Njets  " << std::endl;
      // get the probability vector for n tagged scenarios
      // n jet tag probabilities
      const std::vector<double>    nTagProbVector   = masterJetCalc.getNJetTaggedVector(event, maxJetTags);
      // single jet probabilities
      if (debug > 2) std::cout << "Checking probability vector Njets complete  " << std::endl;
      if (debug > 2) std::cout << "Checking probability vector Jets  " << std::endl;
      const std::vector<double>    jetProbabilityVector	    = 	masterJetCalc.getJetProbabilityVector(event);
      // errors +/-
      const std::vector<std::pair<double,double>> nTagProbErrVector = 	masterJetCalc.getNJetErrorVector(event, maxJetTags);
      if (debug > 2) std::cout << "Checking probability vector Jets complete  " << std::endl;
      if(debug > 2) std::cout << "Getting Vector Size for event: "  << event << std::endl;
      int nJets = tree->GetLeaf("nCaloJets")->GetValue(0);//nTagProbVector.size();
      float pdfWeight	     = 0;

      if(isSig) {
	pdfWeight = fabs(tree->GetLeaf("genWeightsRMS")->GetValue(0)) / fabs(tree->GetLeaf("genWeight")->GetValue(0)) ;//nTagProbVector.size();
	pdfWeightTotal += 1;
	pdfWeightTotalUp += (1 + pdfWeight);
	pdfWeightTotalDn += (1 - pdfWeight);
      }

      // fill the array with zeros
      for(int jj = 0; jj <= maxJetTags; ++jj) probNTags[jj] = 0;
      
      // look at up to maxJetTags
      for(int jj = 0; (jj <= nJets) && (jj <= maxJetTags); ++jj ) {
	if(debug > 2) std::cout << "Checking for jets jj "  << jj <<  " nJets " << nJets << std::endl;
	double	prob   = nTagProbVector[jj];
	double	probUp = nTagProbErrVector[jj].first;
	double	probDn = nTagProbErrVector[jj].second;

	probNTags[jj]	   = (prob >= 0 && prob <= 1) ? prob * totalWeight : 0;
	probNTagsErrUp[jj] = (probUp >= 0 && probUp <= 1) ? probUp * totalWeight : 0;
	probNTagsErrDn[jj] = (probDn >= 0 && probDn <= 1) ? probDn * totalWeight : 0;	

	if(debug> 5)  std::cout << "probTagsErrUp[jj] jj=" << jj << " probUp=" << probUp << std::endl;
	if(debug> 5)  std::cout << "probTagsErrDn[jj] jj=" << jj << " probUp=" << probUp << std::endl;
	   
	double weightErrUp = (probNTags[jj] + probNTagsErrUp[jj]);
	double weightErrDn = probNTags[jj] - probNTagsErrDn[jj];

	if (weightErrDn < 0) { 
	  weightErrDn = 0;
	}
	if (weightErrUp > 1) {
	  weightErrUp = 1;
	}

	// fill the prediciton histogram
	if(debug> 2) {
	  std::cout << "Fill histograms with probabilities p = " << 
	    prob << " pUp =  " << probUp << " pDn" << probDn << 
	    " weight up " << weightErrUp << " weight dn " << weightErrDn << 
	    " jj " << jj << std::endl;
	}

	// check that the weight make sense
	if(probNTags[jj] > weightErrUp || probNTags[jj] < weightErrDn || weightErrUp < 0 || weightErrDn < 0) {
	  std::cout << "\nWARNING: NON-SENSICAL ERRORS FOR PROBABILITIES" << std::endl;
	  std::cout << "probNTags" << probNTags[jj] << std::endl;
	  std::cout << "weigherrUp" << weightErrUp << std::endl;
	  std::cout << "weigherrDn" << weightErrDn << std::endl;
	}

	// fill histograms and add a factor 2 if we are estimating the background using the 2 sample devision
	if(eventPassSelection) {
	  if(debug > 2) std::cout << "Event passes selectio nso Fill histograms with probabilities " << std::endl;
	  if(jj >= 2 && nDivisions == 2 && runChop && removeSignalRegionFromProbs) {
	    nTagHistPred.Fill(jj, probNTags[jj]*2);
	    nTagHistPredErrUp.Fill(jj, weightErrUp*2);
	    nTagHistPredErrDn.Fill(jj, weightErrDn*2);
	  }
	  else {
	    nTagHistPred.Fill(jj, probNTags[jj]);
	    nTagHistPredErrUp.Fill(jj, weightErrUp);
	    nTagHistPredErrDn.Fill(jj, weightErrDn);
	  }
	} // close if the event passes the event selection

	// fill the systematics histograms for the signal		
	if(debug > 3) { 
	  std::cout << " histogram status " << std::endl;
	  nTagHistPred.Print();	
	  std::cout << "\tintegral =  " << nTagHistPred.Integral() << std::endl;
	  nTagHistPredErrUp.Print();
	  std::cout << "\tintegral =  " << nTagHistPredErrUp.Integral() << std::endl;
	  nTagHistPredErrDn.Print();
	  std::cout << "\tintegral =  " << nTagHistPredErrDn.Integral() << std::endl;
	}

	if(debug > 5 && nJets > 0) {
	  std::cout << "probJTag[jj]: " << probJTag[jj] << std::endl;
	  std::cout << "jetProbabilityVector[jj]: " << jetProbabilityVector[jj] << std::endl;
	}
	
	// get the jet probability
	if(nJets > 0){
	  probJTag[jj] = jetProbabilityVector[jj];
	}
	if(debug > 3) 	  std::cout << "------Jet loop interation complete---- " << std::endl;
      }

      // set the weight vector to the number of true tags 
      nTagWeight[nTagged] = 1;
      // fill the histogram
      if(eventPassSelection) { 
	nTagHistTrue.Fill(nTagged, 1);	
	nTagHistTruePDFUp.Fill(nTagged, (1 + pdfWeight));
	nTagHistTruePDFDn.Fill(nTagged, (1 - pdfWeight));
      }

      if(isSig) {
	// counts from the smeared tagging systematic
	if(eventPassSelection) {
	  nTagHistTrueTrackingSYSUp.Fill(nTaggedUp, 1);
	  nTagHistTrueTrackingSYSDn.Fill(nTaggedDn, 1);
	}

	// build the trigger matching conditions explicitly 
	// when 2dip and 2dip sig are varied up the triggers are more likely to fire
	// nominal values
	bool passInclusive  =  caloHT > 650 && nInclusiveJets >= 2;
	bool passDisplaced  =  caloHT > 450 && nDisplacedTrackJets >= 2;
	// inclusive varied up/dn
	bool passInclusiveUp = caloHT > 650 && nInclusiveJetsUp >= 2; // && nInclusiveJets < 2;
	bool passInclusiveDn = caloHT > 650 && nInclusiveJetsDn >= 2; // && nInclusiveJets < 2;
	// displaced varied up/dn
	bool passDisplacedUp = caloHT > 450 && nDisplacedTrackJetsUp >= 2;// && nDisplacedTrackJetsUp < 2;
	bool passDisplacedDn = caloHT > 450 && nDisplacedTrackJetsDn >= 2;// && nDisplacedTrackJetsDn >= 2;

	// build an OR of the trigger matching conditions
	bool eventPassHLT   = passInclusive || passDisplaced;
	bool eventPassHLTUp = passInclusiveUp || passDisplacedUp;
	bool eventPassHLTDn = passInclusiveDn || passDisplacedDn;
	
	// fill the histograms
	if(eventPassHLT) nTagHistTrueHLT.Fill(nTagged, 1);
	if(eventPassHLTUp) nTagHistTrueHLTSYSUp.Fill(nTagged, 1);
	if(eventPassHLTDn) nTagHistTrueHLTSYSDn.Fill(nTagged, 1);	    
	//nTagHistTrue.Fill(nTagged, totalWeight);
      }

      if(debug > 1) std::cout << " -- Number of jets tagged...." << nTagged << std::endl;

      // fill the Tree
      if(debug > 2) std::cout << "Filling the Tree.." << std::endl;
      if(writeTree && eventPassSelection) jetVariableTree->Fill();

      // reset the weight back to zero
      nTagWeight[nTagged] = 0;
      if(debug > 2) std::cout << "leaving event loop scope.." << std::endl;
    }
    // //////////////
    // END loop over events in a single sample
    // //////////////

    if(debug > 2) std::cout << "checking for signal contamin.." << std::endl;
    // make copies this far that have no inclusion of the contamination histograms
    TH1D * nTagHistPred_noContam = (TH1D*)nTagHistPred.Clone((std::string(nTagHistPred.GetName()) + "_noContam").c_str());
    TH1D * nTagHistTrue_noContam = (TH1D*)nTagHistPred.Clone((std::string(nTagHistTrue.GetName()) + "_noContam").c_str());    

    // make a copy of the histogram to only keep track of the histograms with contamation
    TH1D * nTagHistTrue_contam = (TH1D*)nTagHistPred.Clone((std::string(nTagHistTrue.GetName()) + "_contam").c_str());    
    
    //nTagHistTrue_contam->SetTitle("Signal Inj.");
    // clear out the information
    nTagHistTrue_contam->Reset();

    // now add the histograms for the contamination correctly scaled
    if(runSignalContam) {
      if(debug > 2) std::cout << "[contamination] building global probabilities for contamination" << std::endl;
      
      // get the contamation sample tree
      std::cout << "[contamination] Looking for contamination file with path " << contamPath << std::endl;
      TFile	contamFile(contamPath.c_str(), "READ");
      TTree*	contamTree = (TTree*)contamFile.Get(treeName.c_str());
      //      contamTree->Print();
      TTree*	contamTreeHT = (TTree*)contamFile.Get("eventInfo");

      // build the master jet probability comupter using hte global probaiblities and the contamination
      jetProbabilityMasterComputer contamJetCalc(globalJetProbToApply, contamTree, debug);
      
      // weight the events by the normalization determined by the run configuration
      if(debug > 2) std::cout << "[contamination] Parsing scale factors for normalization " <<std::endl;
      float norm    = run_config_root["signalContam"].get("norm",0).asFloat();
      if(debug > 2) std::cout << "[contamination] Parsing number of events in contamination tree " << std::endl;
      int   nEvents = contamTree->GetEntries();
      float scale_factor = norm / float(nEvents);

      if(debug > -1) std::cout << "[contamination]Processing Contamination to add " <<std::endl;
      for(long int event = 0; event < nEvents; ++event) {
	if(debug > 2) std::cout << "[contamination] Checking event selection... "  <<  std::endl;	

	// get the probabilities of n-tags for the signal
	std::vector<double> nTagProbVector = contamJetCalc.getNJetTaggedVector(event, maxJetTags);
	// get the actual number of tags per event
	std::vector<bool>   taggedVector   = jetSel.getJetTaggedVector(contamTree, event, false, false);      
	int		    nJets	   = contamTree->GetLeaf("nCaloJets")->GetValue(0);	

	int nTagged   = 0;

	// loop over the various possibilities of jet tags up to max tags
	for(int jj = 0; (jj <= nJets) && (jj <= maxJetTags); ++jj ) {
	  if(debug > 2) std::cout << "[CONTAMINATION] Checking for jets jj "  << jj <<  " nJets " << nJets << std::endl;
	  if (nJets == 0) break;
	  // contamination only
	  if(taggedVector[jj]) nTagged++;
	} // loop over the possible number of jets in the event

	// check the index matches for the validation sample and that the kinematic / trigger requirements are satisfied
	bool	passValidationIndex	     = (int(event) % nDivisions == valiIndex) || !runChop;   // is the correct validation sample or not running chop
	// is in the signal region for	    <= 2 divisions
	bool	passSignalTagRegion	     = (nTagged >= 2 && nDivisions <= 2);  
	// passes the trigger and kinematic requirements
	bool	eventPassKinematicSelection  = jetSel.doesEventPassSelection(tree, event);
	// two ways to get into the validation sample (signal region or correct validation index)
	bool	eventPassSelection	     = eventPassKinematicSelection && (passValidationIndex || passSignalTagRegion);        
	if(!eventPassSelection) continue;
	
	// fill the number taged for the (samp + contam) and contam only
	nTagHistTrue.Fill(nTagged, scale_factor);
	nTagHistTrue_contam->Fill(nTagged, scale_factor);    

      } // loop over events in the contamination sample       
    } // if statement closing if we are running the contamination study

    // --------------------------
    // write a json result containing the n-tags and systematic variations 
    // based on the histograms
    // --------------------------

    Json::Value resultJSON;
    // build the array
    Json::Value nTagTrue(Json::arrayValue);
    Json::Value nTagPred(Json::arrayValue);
    // statistical uncertainty in fake rate
    Json::Value nTagFakeRateUp(Json::arrayValue);
    Json::Value nTagFakeRateDn(Json::arrayValue);
    // pdf weights
    Json::Value nTagTruePDFUp(Json::arrayValue);
    Json::Value nTagTruePDFDn(Json::arrayValue);
    Json::Value nTagTruePDFRelUp(Json::arrayValue);
    Json::Value nTagTruePDFRelDn(Json::arrayValue);
    // online tracking modeling
    Json::Value nTagTrueHLT(Json::arrayValue);
    Json::Value nTagTrueHLTSYSUp(Json::arrayValue);
    Json::Value nTagTrueHLTSYSDn(Json::arrayValue);
    Json::Value nTagTrueHLTSYSRelUp(Json::arrayValue);
    Json::Value nTagTrueHLTSYSRelDn(Json::arrayValue);
    // tag variable modeling
    Json::Value nTagTrueTrackingSYSUp(Json::arrayValue);
    Json::Value nTagTrueTrackingSYSDn(Json::arrayValue);
    Json::Value nTagTrueTrackingSYSRelUp(Json::arrayValue);
    Json::Value nTagTrueTrackingSYSRelDn(Json::arrayValue);

    // set the array values
    for(int ntag = 1; ntag <= maxJetTags; ++ntag) {
      // true and predicted
      nTagTrue.append(nTagHistTrue.GetBinContent(ntag));
      nTagPred.append(nTagHistPred.GetBinContent(ntag));
      // fake rate 
      nTagFakeRateUp.append(nTagHistPredErrUp.GetBinContent(ntag));
      nTagFakeRateDn.append(nTagHistPredErrDn.GetBinContent(ntag));
      // signal systematics 
      // hlt
      float val0  = nTagHistTrueHLT.GetBinContent(ntag);
      float valUp = nTagHistTrueHLTSYSUp.GetBinContent(ntag);
      float valDn = nTagHistTrueHLTSYSDn.GetBinContent(ntag);
      // absolute values for the HLT Systematic
      nTagTrueHLT.append(val0);
      nTagTrueHLTSYSUp.append(valUp);
      nTagTrueHLTSYSDn.append(valDn);
      // relative systematic bin by bin 
      nTagTrueHLTSYSRelUp.append(val0 > 2 ? fabs(val0-valUp) / val0: -1);
      nTagTrueHLTSYSRelDn.append(val0 > 2 ? fabs(val0-valDn) / val0 : -1);
      // tracking
      val0  = nTagHistTrue.GetBinContent(ntag);
      valUp = nTagHistTrueTrackingSYSUp.GetBinContent(ntag);
      valDn = nTagHistTrueTrackingSYSDn.GetBinContent(ntag);
      // absolute values
      nTagTrueTrackingSYSUp.append(valUp);
      nTagTrueTrackingSYSDn.append(valDn);
      // relative systematic to the nominal value
      nTagTrueTrackingSYSRelUp.append(val0 > 2 ? fabs(val0 - valUp) / val0 : -1);
      nTagTrueTrackingSYSRelDn.append(val0 > 2 ? fabs(val0 - valDn) / val0 : -1);
      // pdf systematic
      val0  = nTagHistTrue.GetBinContent(ntag);
      valUp = nTagHistTruePDFUp.GetBinContent(ntag) * pdfWeightTotal / pdfWeightTotalUp;

      valDn = nTagHistTruePDFDn.GetBinContent(ntag) * pdfWeightTotal / pdfWeightTotalDn;
      nTagTruePDFUp.append(valUp);
      nTagTruePDFDn.append(valDn);
      nTagTruePDFRelUp.append(val0 > 2 ? fabs(val0 - valUp) / val0: -1);
      nTagTruePDFRelDn.append(val0 > 2 ? fabs(val0 - valDn) / val0: -1);      
    }

    // set the values inside the result JSON
    // flags realted to the type of file being analyzed
    resultJSON["label"]				      = label;
    resultJSON["includesContamination"]		      = runSignalContam ? "True" : "False";
    resultJSON["isSignal"]			      = isSig ? "True" : "False";
    resultJSON["x_limit_label"]			      = x_limit_label;
    resultJSON["y_limit_label"]			      = y_limit_label;
    // statistics about events analysed
    resultJSON["nEventsAnalyzed"]		      = total_processed;
    resultJSON["nEventsInFile"]			      = nEvents;
    resultJSON["nEventsPassSelection"]		      = nPassEventSelection;
    resultJSON["nEventsPassTriggerSelection"]	      = nPassTriggerSelection;
    resultJSON["nTagTrue"]			      = nTagTrue;
    resultJSON["nTagTrueHLT"]			      = nTagTrueHLT;
    resultJSON["nTagPred"]			      = nTagPred;
    resultJSON["pdfWeightTotal"]		      = pdfWeightTotal;
    resultJSON["pdfWeightTotalUp"]		      = pdfWeightTotalUp;
    resultJSON["pdfWeightTotalDn"]		      = pdfWeightTotalDn;
    // systematics up and down
    resultJSON["systematics"]["nTagFakeRateUp"]	      = nTagFakeRateUp;
    resultJSON["systematics"]["nTagFakeRateDn"]	      = nTagFakeRateDn;
    // online tracking
    resultJSON["systematics"]["nTagHLTSYSUp"]	      = nTagTrueHLTSYSUp;
    resultJSON["systematics"]["nTagHLTSYSDn"]	      = nTagTrueHLTSYSDn;
    resultJSON["systematics"]["nTagHLTSYSRelUp"]      = nTagTrueHLTSYSRelUp;
    resultJSON["systematics"]["nTagHLTSYSRelDn"]      = nTagTrueHLTSYSRelDn;
    // offline tagging variable modeling
    resultJSON["systematics"]["nTagTrackingSYSUp"]    = nTagTrueTrackingSYSUp;
    resultJSON["systematics"]["nTagTrackingSYSDn"]    = nTagTrueTrackingSYSDn;
    resultJSON["systematics"]["nTagTrackingSYSRelUp"] = nTagTrueTrackingSYSRelUp;
    resultJSON["systematics"]["nTagTrackingSYSRelDn"] = nTagTrueTrackingSYSRelDn;
    resultJSON["systematics"]["nTagTruePDFUp"]        = nTagTruePDFUp;
    resultJSON["systematics"]["nTagTruePDFDn"]        = nTagTruePDFDn;
    resultJSON["systematics"]["nTagTruePDFRelUp"]     = nTagTruePDFRelUp;
    resultJSON["systematics"]["nTagTruePDFRelDn"]     = nTagTruePDFRelDn;

    // keep some information about the probabilities JSON used for this
    resultJSON["probabilitiesProvided?"]	= probProvided ? "True" : "False";
    if(probProvided) resultJSON["providedPath"] = globalProb_json;
    resultJSON["jetTagString"]			= globalJetProbToApply->jetCutString;
    resultJSON["eventCutString"]		= globalJetProbToApply->eventCutString;
    resultJSON["triggerCutOnlyString"]		= globalJetProbToApply->triggerCutOnlyString;
    resultJSON["baselineJetCutString"]		= globalJetProbToApply->baselineJetCutString;
    resultJSON["fakeRateBinningVar"]		= globalJetProbToApply->binningVar;
    resultJSON["xsec"]                          = xsec;

    // write teh result JSO using the styled writter
    std::string jsonResultOutputName = outputDir+"/"+"result_"+label;    
    if(runChop) { 
      jsonResultOutputName+=runChopSuffix;
    }      
    jsonResultOutputName += ".json";


    std::ofstream jsonResult_stream;
    jsonResult_stream.open(jsonResultOutputName);      
      // use the style writter
    Json::StyledWriter styledWriter;
    jsonResult_stream << styledWriter.write(resultJSON);
    jsonResult_stream.close();

    // --------------------------
    // make cosmetic changes to the histograms and build a canvas with all
    // contamination and esimtations  
    // --------------------------
    sampleOutputFile.cd();    

    // build the expectation comparison
    TCanvas canvas("canvas", "prediction vs number of tags per event", 800, 800);

    // lines attributes
    nTagHistPred.SetLineWidth(2);
    nTagHistTrue.SetLineWidth(2);
    nTagHistTrue_contam->SetLineWidth(2);
    nTagHistPredErrUp.SetLineStyle(9);
    nTagHistPredErrDn.SetLineStyle(9);
    nTagHistTrue_contam->SetLineStyle(9);
    nTagHistPred.SetLineColor(kBlue);
    nTagHistTrue.SetLineColor(kBlack);
    nTagHistTrue_contam->SetLineColor(kRed);

    // marker  attributes
    nTagHistPred.SetMarkerStyle(25);
    nTagHistTrue.SetMarkerStyle(20);
    nTagHistTrue.SetMarkerSize(1.5);
    nTagHistTrue_contam->SetMarkerStyle(20);
    nTagHistPred.SetMarkerColor(kBlue);
    nTagHistTrue.SetMarkerColor(kBlack);
    nTagHistTrue_contam->SetMarkerColor(kRed);

    // axes
    nTagHistPredErrUp.GetXaxis()->SetTitle("N Jets Tagged");
    nTagHistPredErrUp.GetYaxis()->SetTitle("N Events");

    // draw onto the canvas
    nTagHistPredErrUp.Draw();
    nTagHistPredErrDn.Draw("same");
    nTagHistPred.Draw("psame");
    nTagHistTrue.Draw("epsame");
    if(runSignalContam) nTagHistTrue_contam->Draw("histsame");

    TLegend *  leg = canvas.BuildLegend();
    leg->SetFillColor(0);      
    leg->SetLineColor(0);      
    
    // write the canvas
    canvas.Write();

    // write the histograms
    nTagHistTrue.Write();
    nTagHistPred.Write();
    nTagHistPredErrUp.Write();
    nTagHistPredErrDn.Write();



    // write the contamination related histograms
    nTagHistPred_noContam->Write();
    nTagHistTrue_noContam->Write();
    nTagHistTrue_contam->Write();

    // write the hlt tracking related systematics histograms
    nTagHistTrueHLT.Write();
    nTagHistTrueHLTSYSUp.Write();
    nTagHistTrueHLTSYSDn.Write();
    nTagHistTrueTrackingSYSUp.Write();
    nTagHistTrueTrackingSYSDn.Write();    
    nTagHistTruePDFUp.Write();
    nTagHistTruePDFDn.Write();

    // write out the tree and close 
    if(writeTree) jetVariableTree->Write();
    sampleOutputFile.Close();    
  } // loop over samples contained in the the sample json  

  if(debug > -1) std::cout << "Closing outputfile..." << std::endl;
  outputFile.Close();  
} // end of main method

