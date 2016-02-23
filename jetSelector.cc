#include "jetSelector.h"

// construct the selector 
jetSelector::jetSelector(const Json::Value & selectorJSON, 
			 const bool & runChop_, 
			 const int & nDivisions_,
			 const int & probIndex_, 
			 const int & valiIndex_,
			 const int & debug_) :
  runChop(runChop_),
  nDivisions(nDivisions_),
  probIndex(probIndex_),
  valiIndex(valiIndex_),
  debug(debug_)
{  

  jetCutString	       = selectorJSON.get("jetSelectionCutString","1").asString();
  eventCutString       = selectorJSON.get("eventSelectionCutString","1").asString();
  baselineJetCutString = selectorJSON.get("baselineJetSelectionCutString","(1)").asString();
  triggerCutOnlyString = selectorJSON.get("triggerCutOnlyString","(1)").asString();


  if(runChop) {
    std::string probIndexPassString = "(evNum %" + std::to_string(nDivisions) + " == " + std::to_string(probIndex) + ")";
    eventCutString		    = "(" + eventCutString + ") &&" + probIndexPassString;
    std::cout << "Appending cut for probability index for Cross Validition:" << probIndexPassString << std::endl; 
    std::cout << "New Cut" << eventCutString << std::endl; 
  }

  // 
  // FAKE RATE SETUP
  //
  
  // get the fake rate binning piece of the json
  const Json::Value fakeRateBinning = selectorJSON["fakeRateBinning"];
  // parse  the names of the fake rate variables
  binningVar	    = fakeRateBinning.get("binningVar","ERROR").asString();
  categoryVar	    = fakeRateBinning.get("categoryVar","none").asString();
  // pull out the binning arrays from the json
  const Json::Value histBins	    = fakeRateBinning["bins"];
  const Json::Value catBins	    = fakeRateBinning["catBins"];


  // fill vectors with the values from the JSON
  // print the values for the binned variable for the fake rate 
  std::cout << "[jetSelector] bin values for the fake rate variable: ";
  for(int ii = 0; ii < histBins.size(); ++ii)  {
    double binValue = histBins[ii].asDouble();
    std::cout << binValue << " ";
    histBinVals.push_back(binValue);
  }
  std::cout << std::endl;

  // print out the bin values for the category variables
  std::cout << "[jetSelector] bin values for the category variable: ";
  for(int ii = 0; ii < catBins.size(); ++ii)  { 
    double binValue = catBins[ii].asDouble();
    std::cout << binValue << " ";
    catBinVals.push_back(binValue);
  }
  std::cout << std::endl;

  //
  // VARIABLES TO SAVE SETUP (FOR SHALLOW TREE COPY)
  //

  // parse the variables from the JSON
  eventVariablesToSave = selectorJSON["eventVariablesToSave"];
  jetVariablesToSave   = selectorJSON["jetVariablesToSave"];

  //
  // PARSE THE JET AND EVENT SELECTION
  //   
  doEven            = selectorJSON.get("doEven", false).asBool();
  jetSelection	    = selectorJSON["jetSelection"];
  eventSelection    = selectorJSON["eventSelection"];
  //  triggerThresholds = selectorJSON["triggerThresholds"];

}

bool  jetSelector::doesEventPassSelection(TTree * tree, long int event) { 

 if(debug > 5) std::cout << "[jetSelector]  getting tree event for event selection " << std::endl; 
 tree->GetEntry(event);
 bool   eventPass = true;
 int	evNum	  = tree->GetLeaf("evNum")->GetValue(0);
 bool	isEven	  = (evNum % 2 == 0);

 bool passValidationIndex  = (evNum % nDivisions) == valiIndex;

 if(!passValidationIndex && runChop) return false;

 if(debug > 2) std::cout << "[jetSelector]  is Event EVEN?  " << isEven << " evNum " << evNum << std::endl; 
 if(debug > 2) std::cout << "[jetSelector]  starting event variable loop with event  " << event << std::endl; 
 // check the trigger thresholds first
 // for(int ii = 0; ii < triggerThresholds.size(); ++ii) {
 //   std::string triggerName = triggerThresholds[ii].asString();
   
 // }
 
 // assume all event variables are a single float variable
 for(int ii = 0; ii < eventSelection.size(); ++ii) {
   // check if the variable is an OR (for triggers mostly) 
   bool	isOR   = eventSelection[ii].get("isTriggerOR", false).asBool();

   if(isOR) {
     bool didPass = false;
     // get the list of variables, mins, and maxs for each
     Json::Value    variables	 = eventSelection[ii]["variables"];
     Json::Value    mins	 = eventSelection[ii]["mins"];
     Json::Value    maxs	 = eventSelection[ii]["maxs"];
     Json::Value    htThresholds = eventSelection[ii]["htThresholds"];

     // loop over each variable in the OR
     for(int var = 0; var < variables.size(); ++var){
       std::string triggerName = variables[var].asString();

       tree->GetEntry(event); // call this AGAIN!
       float	val    = tree->GetLeaf(triggerName.c_str())->GetValue(0);
       double   ht     = tree->GetLeaf("eventCaloHT")->GetValue(0);

       // parse the thresholds
       float htThreshold = htThresholds[var].asFloat();
	
       // bool is500 = variables[var].asString() == "passDisplaced500_40";
       // bool is350 = variables[var].asString() == "passDisplaced350_40";       

       // parse the individual boundaries for this piece of the OR
       float min = mins[var].asFloat();
       float max = maxs[var].asFloat();      

       bool pass = val >= min && val <= max && ht > htThreshold;
       if(doEven) pass = pass && isEven;

       if(debug > 5) std::cout << "[jetSelector] Checking HT threshold for trigger OR name=" << 
		       variables[var].asString() << " HT = " << ht << " trigger pass ? " << val  <<
		       " threshold and trigger pass? " << pass << std::endl;										


       // only one of the variables needs to pass in an OR
       if(pass) didPass = true;
     }

     if(debug >5) std::cout << "[jetSelector] did the event pass the trigger or with thresholds?" << didPass << std::endl;

     // if none of the variables in the loop passed...the OR fails
     if (!didPass) return false;
   }
   // otherwise the variable is an AND with the rest
   else {
     if(debug > 5) std::cout << "[jetSelector]  retreiving leaf values... " << std::endl; 
     std::string	var	= eventSelection[ii].get("variable","ERROR").asString();
     if(debug > 5) std::cout << "[jetSelector]  leaf name:... " << var << std::endl; 

     TLeaf *	    varLeaf = tree->GetLeaf(var.c_str());
     float	    val	    = varLeaf->GetValue(0);   
 


     // min and max values for the variable
     const float    min	    = eventSelection[ii].get("min","ERROR").asFloat();
     const float    max	    = eventSelection[ii].get("max","ERROR").asFloat();
     
     if(debug > 1) std::cout << "[jetSelector] Checking EVENT variable: " << var << " min " << 
		   min << " max " << max << " val " << val << std::endl;  
     
     bool fail = val < min || val > max;
     if(doEven) fail = fail || !isEven;

     if(debug > 5) std::cout << "[jetSelector]  Event pass....? " << !fail << std::endl;

     if(fail) return false;
   } // close if/else for performing OR or AND of requirement
 } // end loop over each piece of the event selection

 // if the event never fails it passes
 return true;
}


// given a tree and event produce a vector of whether the jet was tagged or not
std::vector<bool> jetSelector::getJetTaggedVector(TTree * tree, int event) {

  std::vector<bool> isTaggedVec; 
  tree->GetEntry(event);

  if(debug > 1) std::cout << "[jetSelector ] Getting N Jets from nCaloJets" << std::endl;  
  // determine the number of jets in the event to iterate
  TLeaf *   nJetLeaf = tree->GetLeaf("nCaloJets");
  int	    nJets    = nJetLeaf->GetValue(0);

  if(debug > 1) std::cout << "[jetSelector ] Begin looping event jets" << std::endl;  
  for(int jet = 0; jet < nJets; ++jet) {
    if(debug > 1) std::cout << "[jetSelector ] New Jet " << jet << "...egin looping selection variables" << std::endl;  
    for(int ii = 0; ii < jetSelection.size(); ++ii) {

      bool  isRatio = jetSelection[ii].get("isRatio",false).asBool();

      float val	    = -99999;
      
      if(isRatio) {
	// parse the numerator and denominator for the ratio
	std::string num	    = jetSelection[ii].get("num","ERROR").asString();
	std::string den	    = jetSelection[ii].get("den","ERROR").asString();
	if(debug > 2) std::cout << "[jetSelector ] Variable is ratio: " << num << "/" << den << std::endl;  

	TLeaf *     numLeaf = tree->GetLeaf(num.c_str());
	TLeaf *     denLeaf = tree->GetLeaf(den.c_str());
	
	// cross check that the two arrays have the same number of jets
	int           nJetsNum   = numLeaf->GetNdata();
	int           nJetsDen   = denLeaf->GetNdata();
	if (nJetsNum != nJetsDen) {
	  std::cout << "[jetSelector] ERROR -- jet variable arrays do not match for ratio" << std::endl;
	  exit(1);
	}
	
	float	numVal = numLeaf->GetValue(jet);
	float	denVal = denLeaf->GetValue(jet);
	val	       = numVal / denVal;
      }
      else {
	std::string var	    = jetSelection[ii].get("variable","ERROR").asString();
	if(debug > 2) std::cout << "[jetSelector] Variable is not ratio: " << var << std::endl;  
	TLeaf *	    varLeaf = tree->GetLeaf(var.c_str());
	int	    nJets   = varLeaf->GetNdata();
	
	val = varLeaf->GetValue(jet);
      } // end is  not ratio
      if(debug > 2) std::cout << "[jetSelector ] check if variable falls within min and max  " << std::endl;        

      // min and max values for the variable
      const float   min		   = jetSelection[ii].get("min","ERROR").asFloat();
      const float   max		   = jetSelection[ii].get("max","ERROR").asFloat();

      if(debug > 2) std::cout << "[jetSelector ] min  " << min << " max " << max << " val " <<val << std::endl;        
      const bool    isLastVariable = (ii == (jetSelection.size() - 1));      
      const bool    passCut	   = (val >= min) && (val <= max);

      // if it fails a cut, 
      if(!passCut) {
	if(debug > 2) std::cout << "[jetSelector ] Jet Fail #  " << jet <<  std::endl;        
	isTaggedVec.push_back(false);
	break;
      }

      if(passCut && !isLastVariable) continue;

      // we are at the last cut and all have passed      
      if(isLastVariable && passCut) { 
	if(debug > 2) std::cout << "[jetSelector ] Last Variable -- Jet pass #  " << jet <<  std::endl;        
	isTaggedVec.push_back(true);	 	
	break;
      }
    } // end loop over jet variables
    if(debug > 2) std::cout << "[jetSelector ] Variable Loop Complete...  " << std::endl;        
  } // end loop over jets
  if(debug > 2) std::cout << "[jetSelector ] Jet Loop Complete...  " << std::endl;        

  return isTaggedVec;
} // end getJetTaggedVector

 
 
TTree* jetSelector::shallowCopyTree(TTree* oldTree) {
  // dont copy everything
  oldTree->SetBranchStatus("*", 0);  
  // set the individual variables to be kept
  for(int ii = 0; ii < eventVariablesToSave.size(); ++ii) {
    oldTree->SetBranchStatus(eventVariablesToSave[ii].asString().c_str(), 1);
  } 
  for(int ii = 0; ii < jetVariablesToSave.size(); ++ii) {
    oldTree->SetBranchStatus(jetVariablesToSave[ii].asString().c_str(), 1);
  } 
  // clone the tree and return
  TTree * newTree = oldTree->CloneTree(0);  
  return newTree;
}



// check if a jet is seleceted against the jetSelector
//bool jetSelector::isJetSelected(const jetCandidate& jet) {
//  return false;
//}
