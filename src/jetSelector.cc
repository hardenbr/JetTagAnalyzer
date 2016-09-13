#include "jetSelector.h"
#include <random>

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

  // initialize the random see
  srand(123456);

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
  for(int ii = 0; ii < int(histBins.size()); ++ii)  {
    double binValue = histBins[ii].asDouble();
    std::cout << binValue << " ";
    histBinVals.push_back(binValue);
  }
  std::cout << std::endl;

  // print out the bin values for the category variables
  std::cout << "[jetSelector] bin values for the category variable: ";
  for(int ii = 0; ii < int(catBins.size()); ++ii)  { 
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


// returns a vector of tuples of track multiplicities (N_0, N+, N- ) indexed by the jets in the event
// if isPrompt then access the prompt pdfs and use the offline prompt track multip.
// if isDisp   then access the disp pdfs and use the offline disp track multip.
// requires input of the lifetime in units of [mm] to select the correct PDFs
std::vector<std::tuple<int,int,int> > jetSelector::buildOnlineTrackingFromJSON(TTree * tree, 
									       long int event, 
									       int lifetime,
									       const Json::Value & pdfJSON, 
									       const bool & isPrompt) {

  std::vector<std::tuple<int,int,int> > resultVector;
  if(debug > 5) std::cout << "[jetSelector]  Generating online tracking using RAWAOD+ PDFs from JSON " << std::endl;   
  // get the event
  tree->GetEntry(event);
  if(debug > 5) std::cout << "[jetSelector]  Tree Entry retrieved  " << event << std::endl;    

  if(debug > 5) std::cout << "[jetSelector] Getting Leaves" << std::endl;   
  // get the relavant leaves 
  TLeaf *   varLeaf    = tree->GetLeaf("nCaloJets");
  // prompt offline track leaves
  TLeaf *   promLeaf   = tree->GetLeaf("jetNTracksPrompt");
  // displaced offline track leaves
  TLeaf *   dispLeaf   = tree->GetLeaf("jetNTracksDisp");

  // loop over the jets in the event  
  float	    nJets   = varLeaf->GetValue(0);    
  if(debug > 5) std::cout << "[jetSelector]  Begin Jet Loop for online tracking nJets " << nJets << std::endl;   
  for(int jet = 0; jet < nJets; ++jet) {
    
    // set thisLeaf to the correct leaf pointer
    TLeaf * thisLeaf   = isPrompt ? promLeaf : dispLeaf;
    int	    nTracks   = thisLeaf->GetValue(jet);

    // get the correct pdfs (either displaced or prompt)
    std::string type   = isPrompt ? "prom" : "disp";
    std::string typeUp = isPrompt ? "promUp" : "dispUp";
    std::string typeDn = isPrompt ? "promDn" : "dispDn";

    // throw a toy given the pdfs and the number of tracks per jet
    if(debug > 5) std::cout << "[jetSelector]  Accessing json arrays and attempt to throw toy " << std::endl;   
    std::string lifetime_str  = std::to_string(lifetime);
    std::string nTracks_str   = std::to_string(nTracks);

    std::pair<int,float>    trackToyPair    = throwPDFToy(pdfJSON[type][lifetime_str][nTracks_str], -1);
    int			    nOnlineTracks   = trackToyPair.first;
    float		    toyValue	    = trackToyPair.second;
    // use the same toy value for the up and dn so the results are correlated
    // that is, the Up and Dn are strictly higher or lower than the nominal value
    std::pair<int,float>    trackToyPairUp  = throwPDFToy(pdfJSON[typeUp][lifetime_str][nTracks_str], toyValue);
    int			    nOnlineTracksUp = trackToyPairUp.first;
    std::pair<int,float>    trackToyPairDn  = throwPDFToy(pdfJSON[typeDn][lifetime_str][nTracks_str], toyValue);
    int			    nOnlineTracksDn = trackToyPairDn.first;

    std::tuple<int,int,int> thisTuple(nOnlineTracks, nOnlineTracksUp, nOnlineTracksDn);    
    if(debug > 5) std::cout << "[jetSelector]  access complete.....pushing back 3 tuple of systematic varied " << std::endl;   
    // add it to the result
    resultVector.push_back(thisTuple);
  } // end loop over jets  
  
  return resultVector;
} // end method for online tracking PDF delivering

bool jetSelector::doesEventPassPDGID(TTree * genTree, long int event, int pid1, int pid2) {
  genTree->GetEntry(event);
  
  TLeaf *   pidLeaf = genTree->GetLeaf("genPartPID");
  TLeaf *   momLeaf = genTree->GetLeaf("genMomPID");
  int	    nPart   = pidLeaf->GetNdata();

  // found the first particle index
  bool	foundOne = false;
  int	index1	 = -1; // index of the first particle found
  int   mom1     = -9999; 
  bool	foundTwo = false;

  // look for the first particle if we find it save the index
  for(int ii = 0; ii < nPart; ++ii) {
    int pid = fabs(pidLeaf->GetValue(ii));
    if (pid == pid1) {
      foundOne = true;
      index1 = ii;
      mom1 = momLeaf->GetValue(ii);
      continue;
    }
  }

  // if we dont find one, we are done 
  if(!foundOne) return false;

  // find the second particle and make sure its not the same index
  // also make sure it comes from separate mothers 
  for(int ii = 0; ii < nPart; ++ii) {
    float   pid	  = fabs(pidLeaf->GetValue(ii));
    float   momID = momLeaf->GetValue(ii);
    if (fabs(pid) == pid2 && (ii != index1) && (momID != mom1)) {
      foundTwo = true;
      continue;
    }
  }

  // return finding both particles
  return (foundOne && foundTwo);   
  
}

// check only the trigger ors 
bool  jetSelector::doesEventPassTriggers(TTree * tree, long int event) { 
 if(debug > 5) std::cout << "[jetSelector]  getting tree event for event selection " << std::endl; 
 tree->GetEntry(event);
 for(int ii = 0; ii < int(eventSelection.size()); ++ii) {
   // check if the variable is an OR (for triggers mostly) 
   bool	isOR   = eventSelection[ii].get("isTriggerOR", false).asBool();

   if(isOR) {
     // get the list of variables, mins, and maxs for each
     Json::Value    variables	 = eventSelection[ii]["variables"];
     Json::Value    mins	 = eventSelection[ii]["mins"];
     Json::Value    maxs	 = eventSelection[ii]["maxs"];
     
     // loop over each variable in the OR
     for(int var = 0; var < int(variables.size()); ++var){
       std::string triggerName = variables[var].asString();
       
       tree->GetEntry(event);	// call this AGAIN!
       float	val  = tree->GetLeaf(triggerName.c_str())->GetValue(0);
       float	min  = mins[var].asFloat();
       float	max  = maxs[var].asFloat();      
       bool	pass = val >= min && val <= max;       
       // only one of the variables needs to pass in an OR

       if(pass) return true;
     }
   }
 }
 return false;     
}
  	     
bool  jetSelector::doesEventPassSelection(TTree * tree, long int event) { 

 if(debug > 5) std::cout << "[jetSelector]  getting tree event for event selection " << std::endl; 
 tree->GetEntry(event);
 
 int	evNum	  = tree->GetLeaf("evNum")->GetValue(0);
 bool	isEven	  = (evNum % 2 == 0);

 // bool passValidationIndex  = (evNum % nDivisions) == valiIndex;
 // if(!passValidationIndex && runChop && nDivisions > 2) return false;

 if(debug > 6) std::cout << "[jetSelector]  is Event EVEN?  " << isEven << " evNum " << evNum << std::endl; 
 if(debug > 6) std::cout << "[jetSelector]  starting event variable loop with event  " << event << std::endl; 
 // check the trigger thresholds first
 // for(int ii = 0; ii < triggerThresholds.size(); ++ii) {
 //   std::string triggerName = triggerThresholds[ii].asString();
   
 // }
  // assume all event variables are a single float variable
 for(int ii = 0; ii < int(eventSelection.size()); ++ii) {
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
     for(int var = 0; var < int(variables.size()); ++var){
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
     
     if(debug > 6) std::cout << "[jetSelector] Checking EVENT variable: " << var << " min " << 
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

std::vector<float> jetSelector::getJetBinningVarVector(TTree * tree, int event) {
  std::vector<float> binningVarVec; 
  tree->GetEntry(event);
  TLeaf *   nJetLeaf = tree->GetLeaf("nCaloJets");
  int	    nJets    = nJetLeaf->GetValue(0);

  for(int jet = 0; jet < nJets; ++jet) {
    if(debug > 6) std::cout << "[jetSelector ] New Jet " << jet << "...egin looping selection variables" << std::endl;  
    
    TLeaf * binningVarLeaf = tree->GetLeaf(binningVar.c_str());
    float binningVarVal = binningVarLeaf->GetValue(jet);
    binningVarVec.push_back(binningVarVal);    
  } // loop over the number of jets

  return binningVarVec;
}

// given a tree and event produce a vector of whether the jet was tagged or not
// if isSmear then parse the smear value from the jet selection and compute accordingly
std::vector<bool> jetSelector::getJetTaggedVector(TTree * tree, int event, const bool & isSmear, const bool & smearUp) {

  std::vector<bool> isTaggedVec; 
  tree->GetEntry(event);

  if(debug > 6) std::cout << "[jetSelector ] Getting N Jets from nCaloJets" << std::endl;  
  // determine the number of jets in the event to iterate
  TLeaf *   nJetLeaf = tree->GetLeaf("nCaloJets");
  int	    nJets    = nJetLeaf->GetValue(0);

  if(debug > 6) std::cout << "[jetSelector ] Begin looping event jets" << std::endl;  
  for(int jet = 0; jet < nJets; ++jet) {
    if(debug > 6) std::cout << "[jetSelector ] New Jet " << jet << "...egin looping selection variables" << std::endl;  
    for(int ii = 0; ii < int(jetSelection.size()); ++ii) {

      bool  isRatio = jetSelection[ii].get("isRatio",false).asBool();
      float val	    = -99999;
      
      if(isRatio) {
	// parse the numerator and denominator for the ratio
	std::string num	    = jetSelection[ii].get("num","ERROR").asString();
	std::string den	    = jetSelection[ii].get("den","ERROR").asString();
	if(debug > 6) std::cout << "[jetSelector ] Variable is ratio: " << num << "/" << den << std::endl;  

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
	if(debug > 6) std::cout << "[jetSelector] Variable is not ratio: " << var << std::endl;  
	TLeaf *	    varLeaf = tree->GetLeaf(var.c_str());	
	val = varLeaf->GetValue(jet);
      } // end is  not ratio

      // if we are smearing the values
      if(isSmear) {
	float	smear_factor = jetSelection[ii].get("smear_factor", 0).asFloat();	  
	float	sign	     = smearUp ? 1 : -1;
	//std::cout << "[jetSelector ] smearing the value:   " << val << " by " << smear_factor << std::endl;        
	val		     = val * (1 + (sign * smear_factor));
        //std::cout << "[jetSelector ] new value:   " << val << std::endl;        
      }

      if(debug > 6) std::cout << "[jetSelector ] check if variable falls within min and max  " << std::endl;        

      // min and max values for the variable
      const float   min		   = jetSelection[ii].get("min","ERROR").asFloat();
      const float   max		   = jetSelection[ii].get("max","ERROR").asFloat();

      if(debug > 6) std::cout << "[jetSelector ] min  " << min << " max " << max << " val " <<val << std::endl;        
      const bool    isLastVariable = (ii == (int(jetSelection.size()) - 1));      
      const bool    passCut	   = (val >= min) && (val <= max);

      // if it fails a cut, 
      if(!passCut) {
	if(debug > 6) std::cout << "[jetSelector ] Jet Fail #  " << jet <<  std::endl;        
	isTaggedVec.push_back(false);
	break;
      }

      if(passCut && !isLastVariable) continue;

      // we are at the last cut and all have passed      
      if(isLastVariable && passCut) { 
	if(debug > 6) std::cout << "[jetSelector ] Last Variable -- Jet pass #  " << jet <<  std::endl;        
	isTaggedVec.push_back(true);	 	
	break;
      }
    } // end loop over jet variables
    if(debug > 6) std::cout << "[jetSelector ] Variable Loop Complete...  " << std::endl;        
  } // end loop over jets
  if(debug > 6) std::cout << "[jetSelector ] Jet Loop Complete...  " << std::endl;        

  return isTaggedVec;
} // end getJetTaggedVector
 
TTree* jetSelector::shallowCopyTree(TTree* oldTree) {
  // dont copy everything
  oldTree->SetBranchStatus("*", 0);  
  // set the individual variables to be kept
  for(int ii = 0; ii < int(eventVariablesToSave.size()); ++ii) {
    oldTree->SetBranchStatus(eventVariablesToSave[ii].asString().c_str(), 1);
  } 
  for(int ii = 0; ii < int(jetVariablesToSave.size()); ++ii) {
    oldTree->SetBranchStatus(jetVariablesToSave[ii].asString().c_str(), 1);
  } 
  // clone the tree and return
  TTree * newTree = oldTree->CloneTree(0);  
  return newTree;
}

 std::pair<int,float> jetSelector::throwPDFToy(const Json::Value & slice, const float & toyVal) {

  // throw a random toy and renomalize
  int	toy	 = std::rand();
  float unit_toy = float(toy) / float(RAND_MAX);

  
  if(toyVal > 0) unit_toy = toyVal;
  
  // proceed through the bins 0 -> MAX and pick the 
  // number of online tracks where the toy falls
  float sum = 0;  
  int nBins = slice.size();
  for(int bin = 0; bin < nBins; ++bin) {
    float val   = slice[bin].asDouble();    
    sum += val;
    if(sum > unit_toy) { 
      std::pair<int,float> pair(bin, unit_toy);
      return pair;
    }
  }              

  // this shouldnt happen as long as pdfs are properly normalized
  std::pair<int,float> bad_pair(-1, -1);
  return bad_pair;
}
