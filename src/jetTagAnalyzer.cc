

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
    globalProb_json = argv[1];
  }
  // when we want to do cross validation the probaiblities method
  else if(argc == 4) {
    sample_json	    = "sample.json";
    selection_json  = "selection.json";
    run_config_json = "run_config.json";
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

  // file to input into the json reader
  std::ifstream sample_doc(sample_json, std::ifstream::binary);
  std::ifstream selection_doc(selection_json, std::ifstream::binary);
  std::ifstream run_config_doc(run_config_json, std::ifstream::binary);
  std::ifstream globalProb_doc(globalProb_json, std::ifstream::binary);

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

    // parse the values to be used for the limit calculation (ONLY IF SIGNAL)
    if(isSig) {
      x_limit_label = samples[ss].get("x_limit_label", -999).asFloat();      
      y_limit_label = samples[ss].get("y_limit_label", -999).asFloat();      
    }

    // build  names for tagging histograms  based on the label
    std::string nTagHistTrueName = label + "_nTagTrue";
    std::string nTagHistPredName = label + "_nTagPred";

    // initialize histograms
    TH1D nTagHistTrue(nTagHistTrueName.c_str(), "N Tags Obs.", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistPred(nTagHistPredName.c_str(), "N Tags Exp.", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistPredErrUp((nTagHistPredName+"ErrUp").c_str(), 
			   "Fake Rate Stat.", maxJetTags, 1-.5, maxJetTags+1-.5);
    TH1D nTagHistPredErrDn((nTagHistPredName+"ErrDn").c_str(), 
			   "Fake Rate Stat.", maxJetTags, 1-.5, maxJetTags+1-.5);
    
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
    TTree * caloHTTree = (TTree*)(thisFile.Get("eventInfo"));    

    std::string contamPath;

    double  total_processed	  = tree->GetEntries();
    double  filtered         	  = tree->GetEntries();

    // get the eventWeight related to the processed piece of the sample
    if(isMC && !isSig) {
      if(debug > 3) std::cout << "Getting the filter count tree ....." << std::endl;
      TTree * filterTree = (TTree*)(thisFile.Get("filtCount"));    
      filterTree->Draw("nEventsTotal>>event_total_hist","nEventsTotal","goff");
      filterTree->Draw("nEventsFiltered>>event_filt_hist","nEventsFiltered","goff");
      TH1F*   totalEventHist	  = (TH1F*)gDirectory->Get("event_total_hist");
      TH1F*   filteredEventHist	  = (TH1F*)gDirectory->Get("event_filt_hist");
      total_processed	  = totalEventHist->Integral();
      filtered         	  = filteredEventHist->Integral();
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


    //get the eventweight for the pilup
    // if(isMC && pileupFileName != "") {
    //   if(debug > 3) std::cout << "Building the pileup Hist for this sample ....." << std::endl;
    //   TH1D* pileupDistMC = (TH1D*)pileupHist.Clone();
    //   std::string histName = "puMC"+std::to_string(ss);
    //   pileupDistMC->SetName(histName.c_str());
    //   pileupDistMC->SetTitle(histName.c_str());
    //   pileupDistMC->Print();
    //   pileupDistMC->Reset();
    //   //tree->Draw(("genPU>>"+histName).c_str(),"","goff");
    //   //pileupDistMC->Print();

    //   if(debug > 3) std::cout << "[PU WEIGHTING] Scaling Pileup Distribution....." << std::endl;
    //   // normalize the two distributions
    //   //pileupDistMC->Scale(1.0 / pileupDistMC->Integral());
      
    //   // get the ratio
    //   if(debug > 3) std::cout << "[PU WEIGHTING] Dividing PU Distributions ....." << std::endl;
    //   pileupWeightHist = (TH1D*)pileupHist.Clone();
    //   pileupWeightHist->Divide(pileupDistMC);
    //   if(debug > 3) std::cout << "[PU WEIGHTING] After Clone and Division ....." << std::endl;
    //   pileupWeightHist->Print();
    // }

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
    // jet indexed
    double probJTag[100];
    int   isTagged[100];
    // event indexed / flat number
    int		nTagged		= 0;
    long int	evNum		= 0;
    float       caloHT		= 0;
    float       subLeadingJetPt = 0;
    float       leadingJetPt	= 0;

    // set the branches for the probabilities
    jetVariableTree->Branch("evNum", &evNum, "evNum/I");

    jetVariableTree->Branch("nCat", &nCat, "nCat/I");
    jetVariableTree->Branch("index", &nJetTagArray, "index[nCat]/I");
    jetVariableTree->Branch("nTagWeight", &nTagWeight, "nTagWeight[nCat]/I");
    jetVariableTree->Branch("probNTags", &probNTags, "probNTags[nCat]/D");
    jetVariableTree->Branch("probNTagsErrUp", &probNTagsErrUp, "probNTagsErrUp[nCat]/D");
    jetVariableTree->Branch("probNTagsErrDn", &probNTagsErrDn, "probNTagsErrDn[nCat]/D");
    // jet indexed
    jetVariableTree->Branch("probJTag", &probJTag, "probJTag[nCaloJets]/D");    
    jetVariableTree->Branch("isTagged", &isTagged, "isTagged[nCaloJets]/I");    
    // event index / flat number
    jetVariableTree->Branch("nTagged", &nTagged, "nTagged/I");    
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

      // Calculate the number of tagged jets
      if(debug > 2) std::cout << "Getting the nTagged Jets Vector.." << std::endl;
      std::vector<bool> taggedVector = jetSel.getJetTaggedVector(tree, event);      
      nTagged = 0; // the total number of jets tagged per event      
      if(debug > 2) std::cout << "Filling output Tree Branch.. with tagvector size:" << taggedVector.size() << std::endl;
      if(debug > 5) std::cout << "[jetTagAnalyzer] vector of tags";
      for(int jj = 0; jj < taggedVector.size(); ++jj) {
	if(debug > 5) std::cout << taggedVector[jj];	
	// check for each vector being tagged
	if(taggedVector[jj]) nTagged++; 
	isTagged[jj] = taggedVector[jj] ? 1 : 0;
      }
      
      // check the index matches for the validation sample and that the kinematic / trigger requirements are satisfied
      bool  passValidationIndex		= (int(event) % nDivisions == valiIndex); // is the correct validation sample
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
      if(event % 5000 == 0) std::cout << "Processing Event # --- "  << event << std::endl;
      if(debug > 9) std::cout << "\t\t\tProcessing Event # --- "  << event << std::endl;
      // get the calo HT calculation 
      caloHTTree->GetEntry(event);
      if(debug > 9) std::cout << "Retrieved Entry" << std::endl;
      TLeaf *   caloHTLeaf	    = caloHTTree->GetLeaf("eventCaloHT");
      int       nHT		    = caloHTLeaf->GetNdata();
      TLeaf *   leadingJetPtLeaf    = caloHTTree->GetLeaf("caloLeadingJetPT");
      int       nPt		    = leadingJetPtLeaf->GetNdata();
      TLeaf *   subLeadingJetPtLeaf = caloHTTree->GetLeaf("caloSubLeadingJetPT");
      int       nSubPt		    = subLeadingJetPtLeaf->GetNdata();
      if(debug > 9) std::cout << "Retrieved Leafs..." << std::endl;
      if(nHT != 1) {
	std::cout << "caloHT leaf does not contain a single entry" << std::endl;	
	exit(1);
      }
      else {
	caloHT		= caloHTLeaf->GetValue(0);
	leadingJetPt	= leadingJetPtLeaf->GetValue(0);
	subLeadingJetPt = subLeadingJetPtLeaf->GetValue(0);
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
      
      // get the probability vector for n tagged scenarios
      // n jet tag probabilities
      std::vector<double>    nTagProbVector		    = 
	masterJetCalc.getNJetTaggedVector(event, maxJetTags);
      // single jet probabilities
      std::vector<double>    jetProbabilityVector	    = 
	masterJetCalc.getJetProbabilityVector(event);
      // errors +/-
      std::vector<std::pair<double,double>> nTagProbErrVector = 
	masterJetCalc.getNJetErrorVector(event, maxJetTags);

      if(debug > 2) std::cout << "Getting Vector Size for event: "  << event << std::endl;
      int nJets = tree->GetLeaf("nCaloJets")->GetValue(0);//nTagProbVector.size();

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

	if(debug> 5 && jj < 3)  std::cout << "probTagsErrUp[jj] jj=" << jj << " probUp=" << probUp << std::endl;
	if(debug> 5 && jj < 3)  std::cout << "probTagsErrDn[jj] jj=" << jj << " probUp=" << probUp << std::endl;
	   
	double weightErrUp = (probNTags[jj] + probNTagsErrUp[jj]);
	double weightErrDn = probNTags[jj] - probNTagsErrDn[jj];

	if (weightErrDn < 0) { 
	  //std::cout << "[WARNING!!!] weight with errors down smaller than 0???" << probNTags[jj] << " errorDn " << probNTagsErrDn[jj] << std::endl;
	  weightErrDn = 0;
	}
	if (weightErrUp > 1) {
	  //std::cout << "[WARNING!!!] weight with errors up greater than 1???" << probNTags[jj] << " errorDn " << probNTagsErrUp[jj] << std::endl;
	  weightErrUp = 1;
	}

	// fill the prediciton histogram
	if(debug> 2) {
	  std::cout << "Fill histograms with probabilities p = " << 
	    prob << " pUp =  " << probUp << " pDn" << probDn << 
	    " weight up " << weightErrUp << " weight dn " << weightErrDn << 
	    " jj " << jj << std::endl;
	}

	// fill histograms
	if(eventPassSelection) {
	  if(debug > 2) std::cout << "Event passes selectio nso Fill histograms with probabilities " << std::endl;
	  nTagHistPred.Fill(jj, probNTags[jj]);
	  nTagHistPredErrUp.Fill(jj, weightErrUp);
	  nTagHistPredErrDn.Fill(jj, weightErrDn);
	}
	
	if(debug > 3) { 
	  std::cout << " histogram status " << std::endl;
	  nTagHistPred.Print();	
	  std::cout << "\tintegral =  " << nTagHistPred.Integral() << std::endl;
	  nTagHistPredErrUp.Print();
	  std::cout << "\tintegral =  " << nTagHistPredErrUp.Integral() << std::endl;
	  nTagHistPredErrDn.Print();
	  std::cout << "\tintegral =  " << nTagHistPredErrDn.Integral() << std::endl;
	}

	if(debug > 5) {
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
      if(eventPassSelection) nTagHistTrue.Fill(nTagged, totalWeight);
      if(debug > 1) std::cout << " -- Number of jets tagged...." << nTagged << std::endl;

      // fill the Tree
      if(debug > 2) std::cout << "Filling the Tree.." << std::endl;
      if(writeTree && eventPassSelection) jetVariableTree->Fill();

      // reset the weight back to zero
      nTagWeight[nTagged] = 0;
    } 
    // //////////////
    // END loop over events in a single sample
    // //////////////

    // make copies this far that have no inclusion of the contamination histograms
    TH1D * nTagHistPred_noContam = (TH1D*)nTagHistPred.Clone((std::string(nTagHistPred.GetName()) + "_noContam").c_str());
    TH1D * nTagHistTrue_noContam = (TH1D*)nTagHistPred.Clone((std::string(nTagHistTrue.GetName()) + "_noContam").c_str());    

    // make a copy of the histogram to only keep track of the histograms with contamation
    TH1D * nTagHistTrue_contam = (TH1D*)nTagHistPred.Clone((std::string(nTagHistTrue.GetName()) + "_contam").c_str());    
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
	std::vector<bool>   taggedVector   = jetSel.getJetTaggedVector(contamTree, event);      
	int		    nJets	   = contamTree->GetLeaf("nCaloJets")->GetValue(0);	

	int nTagged = 0;
	// loop over the various possibilities of jet tags up to max tags
	for(int jj = 0; (jj <= nJets) && (jj <= maxJetTags); ++jj ) {
	  if(debug > 2) std::cout << "[CONTAMINATION] Checking for jets jj "  << jj <<  " nJets " << nJets << std::endl;
	  if (nJets == 0) break;
	  // sample + contamination 
	  nTagHistPred.Fill(jj, nTagProbVector[jj] * scale_factor);	  
	  // contamination only
	  if(taggedVector[jj]) nTagged++;
	} // loop over the possible number of jets in the event

	// check the index matches for the validation sample and that the kinematic / trigger requirements are satisfied
	bool	passValidationIndex	     = (int(event) % nDivisions == valiIndex);   // is the correct validation sample
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
    Json::Value nTagFakeRateUp(Json::arrayValue);
    Json::Value nTagFakeRateDn(Json::arrayValue);

    // set the array values
    for(int ntag = 1; ntag < maxJetTags; ++ntag) {
      nTagTrue.append(nTagHistTrue.GetBinContent(ntag));
      nTagPred.append(nTagHistPred.GetBinContent(ntag));
      nTagFakeRateUp.append(nTagHistPredErrUp.GetBinContent(ntag));
      nTagFakeRateDn.append(nTagHistPredErrDn.GetBinContent(ntag));
    }

    // set the values inside the result JSON
    // flags realted to the type of file being analyzed
    resultJSON["label"]				= label;
    resultJSON["includesContamination"]         = runSignalContam ? "True" : "False";
    resultJSON["isSignal"]                      = isSig ? "True" : "False";
    resultJSON["x_limit_label"]                 = x_limit_label;
    resultJSON["y_limit_label"]                 = y_limit_label;
    // statistics about events analysed
    resultJSON["nEventsAnalyzed"]		= total_processed;
    resultJSON["nEventsInFile"]		        = nEvents;
    resultJSON["nEventsPassSelection"]		= nPassEventSelection;
    resultJSON["nEventsPassTriggerSelection"]   = nPassTriggerSelection;
    resultJSON["nTagTrue"]			= nTagTrue;
    resultJSON["nTagPred"]			= nTagPred;
    resultJSON["systematics"]["nTagFakeRateUp"] = nTagFakeRateUp;
    resultJSON["systematics"]["nTagFakeRateDn"] = nTagFakeRateDn;
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

    // write out the tree and close 
    if(writeTree) jetVariableTree->Write();
    sampleOutputFile.Close();    
  } // loop over samples contained in the the sample json  

  if(debug > -1) std::cout << "Closing outputfile..." << std::endl;
  outputFile.Close();  
} // end of main method


// void processSample(const Json::Value sample) {
  

// }

// void parseSampleJSON() {
// }

// void parseJetSelection() {
// }


// void jetTagAnalyzer::setupOutput() {
// 

// void jetTagAnalyzer::setupProbabilities() {
// }
