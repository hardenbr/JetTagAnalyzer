#include "globalJetProbabilities.h"

globalJetProbabilities::globalJetProbabilities(const std::string& label_, 
					       const std::string& stack,
					       const bool & isMC_ , 
					       const bool & isSig_ , 
					       const double& evWeight_, 
					       const double& xsec_, 
					       Json::Value probabilities,
					       const std::string& divideString_,
					       //const jetSelector& jetSel_,				      
					       const int & debug_) : 
  label(label_),
  isMC(isMC_), 
  isSig(isSig_), 
  evWeight(evWeight_),
  xsec(xsec_),
  divideString(divideString_),
  debug(debug_) {
  //  jetSel(jetSel_){

  std::cout  << "[globalJetProbabilities] USING DIVIDE STRING: " << divideString << std::endl;
  // names for hte histograms generated by the label
  taggedJetHistName = "tagged_" + label;
  allJetHistName    = "allJets_" + label;
  effHistName	    = "effHist_" + label;
  effHistNameUp	    = "effHist_" + label + "Up";
  effHistNameDn	    = "effHist_" + label + "Dn";
  

  if(debug > 2) std::cout << "[globalJetProbabilities] Getting first category" << std::endl;

  // TODO!!!!: generalize to multiple cateogires
  std::string	catName		= "cat0";

  // get the first category
  Json::Value	cat1		  = probabilities[catName];
  // get the arrays for each of the histograms
  Json::Value	binning_json	  = cat1["binning"];
  Json::Value	allHist_json	  = cat1["all"];
  Json::Value	taggedHist_json	  = cat1["tagged"];
  Json::Value	effHist_json	  = cat1["eff"];
  Json::Value	effHistErrUp_json = cat1["effErrUp"];
  Json::Value	effHistErrDn_json = cat1["effErrDn"];

  if(debug > 2) std::cout << "[globalJetProbabilities] Building Binning : ";

  // make a vector of the binning
  std::vector<double> binningVec;
  nBins = int(binning_json.size()) - 1 ;
  // build vectors from the arrays
  for(int bin = 0; bin < int(binning_json.size()); ++bin) {
    double binVal = binning_json[bin].asDouble();    
    if(debug > 2) std::cout << binVal << " ";
    binningVec.push_back(binVal);
  }

  if(debug > 2) std::cout << "\n[globalJetProbabilities] Building Temp Histograms" << std::endl;
  // build the jet histograms
  TH1D allJetHist_temp(allJetHistName.c_str(), "", nBins, &(binningVec[0]));
  TH1D taggedJetHist_temp(taggedJetHistName.c_str(), "", nBins, &(binningVec[0]));
  TH1D ratioHistEff_temp(effHistName.c_str(), "", nBins, &(binningVec[0]));
  TH1D ratioHistEffErrUp_temp(effHistNameUp.c_str(), "", nBins, &(binningVec[0]));
  TH1D ratioHistEffErrDn_temp(effHistNameDn.c_str(), "", nBins, &(binningVec[0]));

  if(debug > 2) std::cout << "[globalJetProbabilities] setting bin content of hists" << std::endl;
  // set the histogram bin content from the json
  for(int edge = 0; edge < nBins; ++edge) {
    int bin = edge + 1;
    allJetHist_temp.SetBinContent(bin, allHist_json[edge].asDouble());
    taggedJetHist_temp.SetBinContent(bin, taggedHist_json[edge].asDouble());
    // the efficiency histogram and the associated errors
    ratioHistEff_temp.SetBinContent(bin, effHist_json[edge].asDouble());
    ratioHistEffErrUp_temp.SetBinContent(bin, effHistErrUp_json[edge].asDouble());
    ratioHistEffErrDn_temp.SetBinContent(bin, effHistErrDn_json[edge].asDouble());
  }

  if(debug > 2) std::cout << "[globalJetProbabilities] Setting Pointers" << std::endl;

  // set the local hists to the temp hists
  taggedJetHist	    = taggedJetHist_temp;
  allJetHist	    = allJetHist_temp;
  ratioHistEff	    = ratioHistEff_temp;
  ratioHistEffErrUp = ratioHistEffErrUp_temp;
  ratioHistEffErrDn = ratioHistEffErrDn_temp;


  // calculate the fake rate
  nJetsTotal	   = allJetHist.Integral();
  nJetsTaggedTotal = taggedJetHist.Integral();
  fakeRate	   = float(taggedJetHist.Integral()) / float(allJetHist.Integral());

  if(debug > 2) std::cout << "[globalJetProbabilities] Parsing cfg strings from config  " << std::endl;
  // set the string variables  from json
  binningVar	       = probabilities[catName].get("binningVar","1").asString();  
  categoryVar	       = probabilities[catName].get("categoryVar","caloJetPt").asString();  
  jetCutString	       = probabilities[catName].get("jetTagString","(1)").asString();
  eventCutString       = probabilities[catName].get("eventTagString","(1)").asString();      
  baselineJetCutString = probabilities[catName].get("baselineJetCutString","(1)").asString();      
  triggerCutOnlyString = probabilities[catName].get("triggerCutOnlyString","(1)").asString();      

  if(debug > 2) {
    std::cout << "\t\t\t binvar: " << binningVar << std::endl;
    std::cout << "\t\t\t catVar: " << categoryVar << std::endl;
    std::cout << "\t\t\t cutstring: " << jetCutString << std::endl;
    std::cout << "\t\t\t eventCutString: " << eventCutString << std::endl;
  }  


}

// create a new global jet probability for outtputing to a JSON
globalJetProbabilities::globalJetProbabilities(const std::string& label_, 
					       const std::string& stack,
					       const bool & isMC_ , 
					       const bool & isSig_ , 
					       const double& evWeight_, 
					       const double& xsec_, 
					       TTree*& tree,
					       jetSelector& jetSel, 
					       const std::string& divideString_,
					       const int & debug_) : 
  label(label_),
  nBins(jetSel.getHistBinning().size() - 1),
  isMC(isMC_), 
  isSig(isSig_), 
  evWeight(evWeight_),
  xsec(xsec_),
  //  jetSel(jetSel_),
  divideString(divideString_),
  debug(debug_)
{

  std::cout  << "[globalJetProbabilities] USING DIVIDE STRING: " << divideString << std::endl;
  
  // parse out the jet selection
  jetCutString	       = jetSel.getJetCutString();
  eventCutString       = jetSel.getEventCutString();
  baselineJetCutString = jetSel.getJetBaselineCutString();
  triggerCutOnlyString = jetSel.getTriggerCutOnlyString();
  // build the fake rate histograms from the binning contained in the jetSelection
  histBinVals	       = jetSel.getHistBinning();
  catBinVals	       = jetSel.getCatBinning();  

  // for now just make the single histogram for the fake Rate
  //  TH1D hist("genericHist", "genericHist", nBins, histBinVals); 
  std::string		    taggedJetHistName = "tagged_" + label;
  std::string		    allJetHistName    = "allJets_" + label;

  // get the names of the variables
  binningVar	   = jetSel.getBinningVarName();
  categoryVar	   = jetSel.getCategoryVarName();

  // form the draw string that will fill the histograms
  std::string	drawSelectedString = binningVar+">>" + taggedJetHistName+"(1000,-100,100)";
  std::string	drawAllString      = binningVar+">>" + allJetHistName + "(1000,-100,100)"; 

  if(debug > -1) { 
    std::cout << "-----------------------" << std::endl;
    std::cout << "binning variable string: " << binningVar << std::endl;
    std::cout << "jet selection string: " << jetCutString << std::endl;
    std::cout << "baseline jet  selection string: " << baselineJetCutString << std::endl;
    std::cout << "event selection string: " << eventCutString << std::endl;
  }

  // fill the histograms with the appropriate draw command
  // -- the tagged jets  
  tree->Draw(drawSelectedString.c_str(), ("(" + eventCutString + ") && (" + jetCutString + ")" + 
					  "&& (" + baselineJetCutString + ")").c_str(), "goff");
  // -- the not tagged jets
  tree->Draw(drawAllString.c_str(), ("(" + eventCutString + ") && (" + baselineJetCutString + ")").c_str(), "goff");

  // get the histograms from the pipe 
  taggedJetHist = (TH1D)*(TH1D*)gDirectory->Get(taggedJetHistName.c_str());
  allJetHist	= (TH1D)*(TH1D*)gDirectory->Get(allJetHistName.c_str());

  // calculate the flat fake rate
  nJetsTotal	   = allJetHist.Integral();
  nJetsTaggedTotal = taggedJetHist.Integral();
  fakeRate	   = float(nJetsTaggedTotal) / float(nJetsTotal);

  // check they are filled
  taggedJetHist.Print();
  allJetHist.Print();

  // rebin the histograms so they can be divided
  taggedJetHist = (TH1D)*(TH1D*)taggedJetHist.Rebin(nBins, "", &(histBinVals[0]));
  allJetHist    = (TH1D)*(TH1D*)allJetHist.Rebin(nBins, "", &(histBinVals[0]));

  // build the corresponding efficieny graph
  ratioGraph.Divide(&taggedJetHist, &allJetHist, divideString.c_str()); 

  std::string graphName = "efficieny_" + label;
  ratioGraph.SetTitle(graphName.c_str());
  ratioGraph.SetName(graphName.c_str());

  if(debug > -1) std::cout << "Building a histogram based on the tgraph " << std::endl;         
  // build histographs based on the ratioGraph for the fake rate
  // copy the all jets histogram  and clear the values
  // make copies off the histogram for the efficiency with the same binning
  ratioHistEff	    = (TH1D)*(TH1D*)allJetHist.Clone(effHistName.c_str());
  ratioHistEff.Reset();
  ratioHistEffErrUp = (TH1D)*(TH1D*)allJetHist.Clone((effHistNameUp).c_str());
  ratioHistEffErrUp.Reset();
  ratioHistEffErrDn = (TH1D)*(TH1D*)allJetHist.Clone((effHistNameDn).c_str());
  ratioHistEffErrDn.Reset();

  // parse the central values and errors of the ratio Graph
  double *  xEffVals	    = ratioGraph.GetX();
  double *  yEffVals	    = ratioGraph.GetY();
  double *  yEffErrorUpVals = ratioGraph.GetEYhigh();
  double *  yEffErrorDnVals = ratioGraph.GetEYlow();
  int	    nPoints	    = ratioGraph.GetN();

  if(debug > -1) std::cout << "Filling the parsed values form the tgraph into histograms " << std::endl;         
  // fill the histogram with the values from the graph
  for(int ii = 0; ii < nPoints; ++ii) {
    // parse the efficiency 
    float   eff = yEffVals[ii];
    float   var = xEffVals[ii];

    // find and fill the bin in the histogram translation
    int	bin   = ratioHistEff.FindBin(var);    
    int	binUp = ratioHistEffErrUp.FindBin(var);    
    int	binDn = ratioHistEffErrDn.FindBin(var);    

    ratioHistEff.SetBinContent(bin, eff);    
    ratioHistEffErrDn.SetBinContent(binDn, yEffErrorDnVals[ii]);    
    ratioHistEffErrUp.SetBinContent(binUp, yEffErrorUpVals[ii]);    
  }   
  
  // check the histogram was printed
  ratioHistEff.Print();
}


Json::Value globalJetProbabilities::getProbabilitiesJSON() {
  // create the json 
  Json::Value event;     

  // loop over all categories for the fake rate
  for(int cat = 0; cat < int(catBinVals.size()); ++cat) {
    // build the arrays for the bins and values
    Json::Value binning(Json::arrayValue);
    Json::Value tagged(Json::arrayValue);
    Json::Value all(Json::arrayValue);
    Json::Value eff(Json::arrayValue);
    Json::Value effErrUp(Json::arrayValue);    
    Json::Value effErrDn(Json::arrayValue);    

    // loop over bins in the histogram for that category 
    for(int bin = 0; bin < int(histBinVals.size()); ++bin) {      
      // histbinvals is the array of the binning and begins with index=0 so use bin
      int   histBin     = bin+1;
      float binVal      = histBinVals[bin];
      // get the values from the histograms 
      // histograms first bin is numbered 1 (0 is the underflow bin)
      float tagVal	= taggedJetHist.GetBinContent(histBin);
      float allJetVal	= allJetHist.GetBinContent(histBin);
      float effVal	= ratioHistEff.GetBinContent(histBin);
      float effValErrUp = ratioHistEffErrUp.GetBinContent(histBin);
      float effValErrDn = ratioHistEffErrDn.GetBinContent(histBin);
      
      // add the into the array 
      binning.append(Json::Value(binVal));
      tagged.append(Json::Value(tagVal));      
      all.append(Json::Value(allJetVal));      
      eff.append(Json::Value(effVal));
      effErrUp.append(Json::Value(effValErrUp));
      effErrDn.append(Json::Value(effValErrDn));
    }
    // add the last bin value for the binning array (which has nbins + 1 entries)
    binning.append(Json::Value(histBinVals[histBinVals.size()-1]));

    std::string catName = "cat" + std::to_string(cat); 

    if(debug > 0) std::cout << "[globalJetProbabities] Using sample xsec for JSON = " << xsec << " pb " << std::endl;
    
    // set info about how this was generated
    event["eventWeight"]		   = evWeight;
    event["xsec"]			   = xsec;
    event[catName]["fakeRate"]	           = fakeRate;
    event[catName]["binningVar"]	   = binningVar;
    event[catName]["binning"]		   = binning;
    event[catName]["jetTagString"]	   = jetCutString;
    event[catName]["eventTagString"]	   = eventCutString;
    event[catName]["triggerCutOnlyString"] = triggerCutOnlyString;
    event[catName]["baselineJetCutString"] = baselineJetCutString;
    // encode the values of the histograms
    event[catName]["tagged"]		   = tagged;
    event[catName]["all"]		   = all;
    event[catName]["eff"]		   = eff;
    event[catName]["effErrUp"]		   = effErrUp;
    event[catName]["effErrDn"]		   = effErrDn;
  }
  
  return event;
} 

// does the lookup in the ratiograph for a specific jet configuration
double globalJetProbabilities::getJetFakeProbability(float binVariable, float catVar) {
  // find the bin and return the central value
  //  ratioHistEff.Print();
  if(debug > 20) std::cout << "[globalJetProb] binVariable: " << binVariable << " catVar " << catVar << std::endl;
  int bin = ratioHistEff.FindBin(binVariable);
  if(debug > 20) std::cout << "[globalJetProb] bin: " << bin << std::endl;
  double value = ratioHistEff.GetBinContent(bin);
  if(debug > 20) std::cout << "[globalJetProb] bin Content: " << value << std::endl;

  if(debug > 5 && binVariable == 1) std::cout << "JET HAS 1 TRACK -- binvar = " << 
				      binVariable << " fake probability = " << value << std::endl;

  return value;
}


// does the lookup in the ratiograph for a specific jet configuration
const std::pair<double, double> globalJetProbabilities::getJetFakeProbabilityError(const float binVariable, const float catVar) {
  // find the bin and return the central value
  //  ratioHistEff.Print();
  if(debug > 5) std::cout << "[globalJetProb] binVariable: " << binVariable << " catVar " << catVar << std::endl;
  int bin = ratioHistEff.FindBin(binVariable);
  if(debug > 5) std::cout << "[globalJetProb] bin: " << bin << std::endl;
  double errUp = ratioHistEffErrUp.GetBinContent(bin);
  double errDn = ratioHistEffErrDn.GetBinContent(bin);
  if(debug > 5) std::cout << "[globalJetProb] bin errors up: " << errUp << " error down: "  <<  errDn << std::endl;

  const std::pair<double, double> errors(errUp, errDn);
  return errors;
}

// does the lookup in the ratiograph for a specific jet configuration
// std::vector<float> globalJetProbabilities::getJetFakeProbabilityVector(TTree* tree, int evNum ) {
//   std::vector<float> probabilityVector;

//   // find the bin and return the central value
//   int bin = ratioHistEff->GetBin(binVariable);
//   // for now dont worry about the category variable
//   return ratioHistEff->GetBinContent(bin);
// }

void globalJetProbabilities::printHistStatus() {  
  // local variables
  std::cout << "----------global jet probabilities to apply object status----------" << std::endl;
  std::cout << "tagged histogram and all jet histogram" << std::endl;
  taggedJetHist.Print();
  allJetHist.Print();
  std::cout << "ratio graph for efficiency" << std::endl;
  ratioGraph.Print();
  std::cout << "efficiency histogram" << std::endl;
  ratioHistEff.Print();
  std::cout << "efficiency error histograms" << std::endl;
  ratioHistEffErrUp.Print();
  ratioHistEffErrDn.Print();
}

void globalJetProbabilities::removeSignalRegion(TTree*& tree,jetSelector & jetSel, bool isContam, float norm) {
  float applyNorm    = 1;

  // change the norm for the signal injection subtraction
  if(isContam) {
    float   signalEvents = tree->GetEntries();    
    applyNorm		 = norm / signalEvents;
  }

  // build names for the new histgrams 
  std::string	appendHistName = isContam ? "isContam" : "noContam";
  std::string	taggedName     = "taggedSignalRegion_" + appendHistName;
  std::string	allName	       = "allSignalRegion_" + appendHistName;

  // initailize the histograms going to be used for the signal region subtraction
  TH1D	signalRegionTagged = (TH1D)*(TH1D*)taggedJetHist.Clone(taggedName.c_str());
  TH1D	signalRegionAll	   = (TH1D)*(TH1D*)allJetHist.Clone(allName.c_str());

  // reset the clones
  signalRegionTagged.Reset();
  signalRegionAll.Reset();

  long int nEvents = tree->GetEntries();
  // loop over all events to find the signal region events to remove
  for(long int event = 0; event < nEvents; ++event) {
    if(debug > 6) std::cout << "Checking event selection... "  <<  std::endl;
    if(event % 5000 == 0) std::cout << "Processing Signal Region  Removal from Probabilities # --- "  << event << std::endl;
    tree->GetEntry(event);

    int nJets = tree->GetLeaf("nCaloJets")->GetValue(0);
    int evNum = tree->GetLeaf("evNum")->GetValue(0);

    // make sure the event passes the event selection with the correct event index
    // or we dont care about the event index because we aren't doing cross validaiton (isContam or divisions =1)
    bool passProbIndex = (evNum % jetSel.nDivisions == jetSel.probIndex) || jetSel.nDivisions == 1 || isContam;
    // if we dont pass the index check and we are, dont use for this probability sample
    if(!passProbIndex) continue; 

    // check the kinematic selection is satisfied
    bool eventPassSelection = jetSel.doesEventPassSelection(tree, event);
    if(!eventPassSelection) continue;
    
    // check each jet if it is tagged
    std::vector<bool> taggedVector = jetSel.getJetTaggedVector(tree, event, false, false);

    // get the total number of taggs in the events
    int	 totalTags = 0;
    for(int ii = 0; ii < nJets; ++ii) totalTags += taggedVector[ii] ? 1 : 0;

    // skip all events not in the signal region
    // we only want to remove the contribution that is in the signal region from the probabilities
    if(totalTags <= 1) continue;

    // this is a vector of the values of the variable used to parameterize the fake rate
    std::vector<float> binningVarVector = jetSel.getJetBinningVarVector(tree, event);
    for(int ii = 0; ii < nJets; ++ii) {
      // zero tracks we can skip
      if(binningVarVector[ii] < 1) continue; 

      // if the jet is taagged fill the tagged hist
      if(taggedVector[ii]) signalRegionTagged.Fill(binningVarVector[ii], applyNorm);      

      // always fil the all hist
      signalRegionAll.Fill(binningVarVector[ii], applyNorm);
    } // endloop over jets in the event               
  } // end loop over events

  // scale the hists to subtract
  signalRegionAll.Scale(-1);
  signalRegionAll.Print();
  signalRegionTagged.Scale(-1);
  signalRegionTagged.Print();

  // check the size of the removal
  std::cout << "[globalProbabilities] N Jets being removed from allJetsHist: " << signalRegionAll.Integral() << std::endl;
  std::cout << "[globalProbabilities] N jets being removed from  tagJetHist: " << signalRegionTagged.Integral() << std::endl;
  std::cout << "allJets before removal: " << allJetHist.Integral() << std::endl;
  std::cout << "tagJets before removal: " << taggedJetHist.Integral() << std::endl;
    

  // add the negative hists to remove 
  taggedJetHist.Add(&signalRegionTagged);
  allJetHist.Add(&signalRegionAll);

  // check for bins that have gone negative
  std::cout << "[globalProb] allJets Hist after signal region removal" << allJetHist.Integral() << std::endl;
  std::cout << "[globalProb] tagJets after signal region removal" << taggedJetHist.Integral() << std::endl;
  
  //sanity check on the bincontents
  for(int ii = 1; ii <= taggedJetHist.GetNbinsX(); ++ii ) {
    float   taggedVal = taggedJetHist.GetBinContent(ii);
    float   allVal    = allJetHist.GetBinContent(ii);


    // remove rounding errors near zero
    if(taggedVal < 0.001 && allVal < 0.001) {
      taggedJetHist.SetBinContent(ii,0);
      allJetHist.SetBinContent(ii,0);
    }      


    if(taggedVal > allVal) {
      std::cout << "tagged bin: " << taggedVal << " allVal " << allVal << std::endl;
      std::cout << "tagged val larger? bin: " << ii << "difference: " << taggedVal - allVal << std::endl;
    }
  }
  

  // reset the ratio, re-perform the division
  ratioGraph.Set(0);

  ratioGraph.Divide(&taggedJetHist, &allJetHist, divideString.c_str()); //clopper pearson

  ratioGraph.SetTitle(("signalRegionRemoved_"+label).c_str());

  //////// NOW RE-FILL THE RATIO HISTOGRAMS ////////////

  // make copies off the histogram for the efficiency with the same binning
  ratioHistEff = (TH1D)*(TH1D*)allJetHist.Clone(effHistName.c_str());
  ratioHistEff.Reset();
  ratioHistEffErrUp = (TH1D)*(TH1D*)allJetHist.Clone((effHistName+"Up").c_str());
  ratioHistEffErrUp.Reset();
  ratioHistEffErrDn = (TH1D)*(TH1D*)allJetHist.Clone((effHistName+"Dn").c_str());
  ratioHistEffErrDn.Reset();

  // parse the central values and errors of the ratio Graph
  double *  xEffVals	    = ratioGraph.GetX();
  double *  yEffVals	    = ratioGraph.GetY();
  double *  yEffErrorUpVals = ratioGraph.GetEYhigh();
  double *  yEffErrorDnVals = ratioGraph.GetEYlow();
  int	    nPoints	    = ratioGraph.GetN();

  if(debug > -1) std::cout << "Filling the parsed values " << std::endl;         
  // fill the histogram with the values from the graph
  for(int ii = 0; ii < nPoints; ++ii) {
    // parse the efficiency 
    float   eff = yEffVals[ii];
    float   var = xEffVals[ii];

    // find and fill the bin in the histogram translation
    int	bin   = ratioHistEff.FindBin(var);    
    int	binUp = ratioHistEffErrUp.FindBin(var);    
    int	binDn = ratioHistEffErrDn.FindBin(var);    

    ratioHistEff.SetBinContent(bin, eff);    
    ratioHistEffErrDn.SetBinContent(binDn, yEffErrorDnVals[ii]);    
    ratioHistEffErrUp.SetBinContent(binUp, yEffErrorUpVals[ii]);    
  } // end loop over points in the ratiograph
}

// add contamination from a signal sample
void globalJetProbabilities::addSignalContamination(TTree*& tree, jetSelector & jetSel, float norm) {

  if(debug > -1) std::cout << "[GlobalProbaibliits] Adding Signal Contamination.....xSecXLumi = " 
			   << norm << std::endl;

  // names for the tagged and untagged jets
  std::string	taggedJetHistName = "sigCont_tagged_" + label;
  std::string	allJetHistName    = "sigCont_allJets_" + label;

  // get the names of the variables
  binningVar			  = jetSel.getBinningVarName();
  categoryVar			  = jetSel.getCategoryVarName();

  // form the draw string that will fill the histograms
  std::string	drawSelectedString = binningVar+">>" + taggedJetHistName+"(1000,-100,100)";
  std::string	drawAllString      = binningVar+">>" + allJetHistName + "(1000,-100,100)"; 

  // draw a pipe into the histograms
  tree->Draw(drawSelectedString.c_str(), ("(" + eventCutString + ") && (" + jetCutString + ")" + 
					  "&& (" + baselineJetCutString + ")").c_str(), "goff");
  tree->Draw(drawAllString.c_str(), ("(" + eventCutString + ") && (" + 
				     baselineJetCutString + ")").c_str(), "goff");

  // go and get the histograms build in the draw command 
  TH1D sigCont_taggedJetHist    = (TH1D)*(TH1D*)gDirectory->Get(taggedJetHistName.c_str());
  TH1D sigCont_allJetHist	= (TH1D)*(TH1D*)gDirectory->Get(allJetHistName.c_str()); 

  // rebin the contamination hists
  sigCont_taggedJetHist = (TH1D)*(TH1D*)sigCont_taggedJetHist.Rebin(nBins, "", &(histBinVals[0]));
  sigCont_allJetHist	= (TH1D)*(TH1D*)sigCont_allJetHist.Rebin(nBins, "", &(histBinVals[0])); 

  // sigCont_taggedJetHist = (TH1D)*(TH1D*)sigCont_taggedJetHist.Rebin(nBins, "", &(histBinVals[0]));
  // sigCont_allJetHist	= (TH1D)*(TH1D*)sigCont_allJetHist.Rebin(nBins, "", &(histBinVals[0])); 


  // store the statistical errors
  sigCont_taggedJetHist.Sumw2();
  sigCont_allJetHist.Sumw2();

  // check contamination hists are filled
  if(debug > -1) std::cout << "histograms after rebinning for signal contamination" << std::endl;
  sigCont_taggedJetHist.Print();
  sigCont_allJetHist.Print();

  // check the historams the contamination is being added to
  if(debug > -1) std::cout << "histograms which the contaminaiton will be added to" << std::endl;
  taggedJetHist.Print();
  allJetHist.Print();

  // print out the bin values for all jet histograms for the contamination
  if(debug > 0) {
    std::cout << "[signalContamination] bin values for tagged hist in contamination (PRE NORMALIZATION)\n " << std::endl;
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << sigCont_taggedJetHist.GetBinContent(bin);
    }
    std::cout << std::endl << "[signal Contamination bin values for all jet hist in contamination (PRE NORMALIZATION)\n ";
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << sigCont_allJetHist.GetBinContent(bin);
    }
    std::cout << "\n";
  }


  // normalize the histograms 
  double nEvents      = tree->GetEntries();
  double scale_factor = norm / nEvents; 
  if(debug > 0) {
    std::cout << "[signalcontamination] applying scale factor to contmaination:" << scale_factor << std::endl;
  }
  sigCont_taggedJetHist.Scale(scale_factor);
  sigCont_allJetHist.Scale(scale_factor);


  // print out the bin values for all jet histograms for the contamination
  if(debug > 0) {
    std::cout << "[signalContamination] bin values for tagged hist in contamination\n " << std::endl;
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << sigCont_taggedJetHist.GetBinContent(bin);
    }
    std::cout << std::endl << "[signal Contamination bin values for all jet hist in contamination\n ";
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << sigCont_allJetHist.GetBinContent(bin);
    }
    std::cout << "\n";
  }


  // print out the bin values for all jet histograms 
  if(debug > 0) {
    std::cout << "[signalContamination] bin values for tagged hist sample (before contamination is added) \n " << std::endl;
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << taggedJetHist.GetBinContent(bin);
    }
    std::cout << std::endl << "[signal Contamination bin values for all jet hist (before contamination is added) \n ";
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << allJetHist.GetBinContent(bin);
    }
    std::cout << "\n";
  }

  // add the contamination histograms to the current histograms
  taggedJetHist.Add(&sigCont_taggedJetHist);
  allJetHist.Add(&sigCont_allJetHist);

  // print out the bin values for all jet histograms for the contamination + sample (cross check of addition working)
  if(debug > 0) {
    std::cout << "[signalContamination] bin values for tagged hist contamination + sample \n" << std::endl;
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << taggedJetHist.GetBinContent(bin);
    }
    std::cout << std::endl << "[signal Contamination bin values for all jet hist in contamination + sample \n ";
    for(int bin = 1; bin <= nBins; ++bin) {
      std::cout << " " << allJetHist.GetBinContent(bin);
    }
    std::cout << "\n";
  }

  // re-calculate the efficiency graph and histogram
  ratioGraph.Set(0);
  ratioGraph.Divide(&taggedJetHist, &allJetHist, divideString.c_str()); 
  std::string graphName = "sigCont_efficieny_" + label;
  ratioGraph.SetTitle(graphName.c_str());
  ratioGraph.SetName(graphName.c_str());

  if(debug > -1) std::cout << "[Signal Contamation]Building a histogram based on the tgraph " << std::endl;         
  // build histographs based on the ratioGraph for the fake rate
  // copy the all jets histogram  and clear the values
  std::string effHistName = "sigCont_effHist_" + label;

  // make copies off the histogram for the efficiency with the same binning
  ratioHistEff	    = (TH1D)*(TH1D*)allJetHist.Clone(effHistName.c_str());
  ratioHistEff.Reset();
  ratioHistEffErrUp = (TH1D)*(TH1D*)allJetHist.Clone((effHistName+"Up").c_str());
  ratioHistEffErrUp.Reset();
  ratioHistEffErrDn = (TH1D)*(TH1D*)allJetHist.Clone((effHistName+"Dn").c_str());
  ratioHistEffErrDn.Reset();

  // parse the central values and errors of the ratio Graph
  double *  xEffVals	    = ratioGraph.GetX();
  double *  yEffVals	    = ratioGraph.GetY();
  double *  yEffErrorUpVals = ratioGraph.GetEYhigh();
  double *  yEffErrorDnVals = ratioGraph.GetEYlow();
  int	    nPoints	    = ratioGraph.GetN();

  if(debug > -1) std::cout << "[Signal Contamation] Filling the parsed values " << std::endl;         
  // fill the histogram with the values from the graph
  for(int ii = 0; ii < nPoints; ++ii) {
    // parse the efficiency 
    float   eff = yEffVals[ii];
    float   var = xEffVals[ii];

    // find and fill the bin in the histogram translation
    int	bin   = ratioHistEff.FindBin(var);    
    
    ratioHistEff.SetBinContent(bin, eff);    
    ratioHistEffErrDn.SetBinContent(bin, yEffErrorDnVals[ii]);    
    ratioHistEffErrUp.SetBinContent(bin, yEffErrorUpVals[ii]);    
  }   
  
  // check the histogram was filled
  ratioHistEff.Print();
}


// accessors
TGraphAsymmErrors globalJetProbabilities::getRatioGraph() { return ratioGraph; }
TH1D globalJetProbabilities::getTaggedHist() { return taggedJetHist; }
TH1D globalJetProbabilities::getAllHist() { return allJetHist; }
TH1D globalJetProbabilities::getCentralEffHist() { return ratioHistEff;}
TH1D globalJetProbabilities::getCentralEffHistErrUp() { return ratioHistEffErrUp;}
TH1D globalJetProbabilities::getCentralEffHistErrDn() { return ratioHistEffErrDn;}
// build from preset global Probabilities
// globalJetProbabilities::globalJetProbabilities(const Json::Value & assignedProb) {  
// }
