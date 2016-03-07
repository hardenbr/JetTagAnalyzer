#include "jetProbabilityMasterComputer.h"
#include <cmath>
// send a pointer to the tree containing the branches to use for the fake rate computation
jetProbabilityMasterComputer::jetProbabilityMasterComputer(globalJetProbabilities * jetProb_, TTree * jetTree_, const int & debug_) {
  std::cout << "Begin Constructor probabilities " << std::endl;
  jetProb = jetProb_;
  std::cout << "\t tree  " << std::endl;
  jetTree = jetTree_;
  std::cout << "\t debug  " << std::endl;
  debug	  = debug_;
}

jetProbabilityMasterComputer::~jetProbabilityMasterComputer() { 
  if(debug > 2) std::cout << " deleting probability master computer " << std::endl;
}

// propagate the statistical error from the fake rate to the n jet tag measurement 
std::vector<std::pair<double, double>> jetProbabilityMasterComputer::getNJetErrorVector(long int eventNumber, int maxJetsTagged) {
  std::string	binVar	= jetProb->getBinningVarName(); 
  if(debug > 2) std:: cout << "[JetProbMaster] error calc...parsing category name" << std::endl;
  std::string	catVar	= jetProb->getCategoryVarName(); 
  if(debug > 2) std:: cout << "[JetProbMaster] error calc...Getting leaf information for binVar and catVar: " 
			   << binVar << " " << catVar <<  std::endl;

  // leaves of the variables
  TLeaf *	binLeaf	= jetTree->GetLeaf(binVar.c_str());
  TLeaf *	catLeaf	= jetTree->GetLeaf(catVar.c_str());
  int		nJets   = binLeaf->GetNdata();

  // setup for parsing from the branch
  std::vector<double> binValues;
  std::vector<double> catValues;

  if(debug > 2) std:: cout << "[JetProbMaster] error calc...filling arrays forbinVar and catVar" << std::endl;
  if(debug > 2) std:: cout << "[JetProbMaster] error calc....nJets in Event = " << nJets << std::endl;
  // fill an array with the values in the ttree
  for(int ii = 0; ii < nJets; ++ii) {
    binValues.push_back(binLeaf->GetValue(ii));
    catValues.push_back(catLeaf->GetValue(ii));
  }
  // use the arrays to compute the configurational probabilities  
  std::vector<std::pair<double,double>> probabilityErrorPairVector;

  // check the tag probability in each scenario
  if(debug > 2) std:: cout << "[JetProbMaster] Looping nTags for event" << std::endl;
  for(int ii = 0; (ii <= nJets) && (ii < maxJetsTagged); ++ii) {
    std::pair<double,double> errors = getNJetProbabilityError(&binValues[0], &catValues[0], ii, nJets);
    probabilityErrorPairVector.push_back(errors);
  }

  // make checks that output is sensible 
  if(debug > 1) { 
    std::cout << " Probability Error Vector for event (jet idx, err_up, err_down): " << std::endl;
    for(int ii = 0; ii <= nJets; ++ii) {
      double up = probabilityErrorPairVector[ii].first;
      double dn = probabilityErrorPairVector[ii].second;
      std::cout << " (" << ii << "," << up << " , " << dn << ") ";
    }
    std::cout << std::endl;      
  }

  return probabilityErrorPairVector;  
}

std::vector<double> jetProbabilityMasterComputer::getJetProbabilityVector(long int eventNumber) {
  // maket he vector we will return with the jet probaiblities
  std::vector<double> jetProbabilityVector; 

  // get the event
  jetTree->GetEntry(eventNumber);

  // get th leaves from the event containing the values
  std::string	binVar	= jetProb->getBinningVarName(); 
  std::string	catVar	= jetProb->getCategoryVarName(); 
  TLeaf *	binLeaf	= jetTree->GetLeaf(binVar.c_str());
  TLeaf *	catLeaf	= jetTree->GetLeaf(catVar.c_str());

  int		nJets   = binLeaf->GetNdata();

  // setup for parsing from the branch
  std::vector<double> binValues;
  std::vector<double> catValues;

  // fill an array with the values in the ttree
  for(int ii = 0; ii < nJets; ++ii) {
    binValues.push_back(binLeaf->GetValue(ii));
    catValues.push_back(catLeaf->GetValue(ii));
  }

  // get the individual jet probabilities
  for(int ii = 0; ii < nJets; ++ii) {
    double prob = jetProb->getJetFakeProbability(binValues[ii], catValues[ii]);
    jetProbabilityVector.push_back(prob);
  }

  return jetProbabilityVector; 
}

// calculate the probability of N tags as a vector from 0 to N jets in the event
std::vector<double> jetProbabilityMasterComputer::getNJetTaggedVector(long int eventNumber, int maxJetsTagged) {
  if(debug > 2) std:: cout << "[JetProbMaster] getting event" << std::endl;
  jetTree->GetEntry(eventNumber);

  if(debug > 2) std:: cout << "[JetProbMaster] parsing binning information" << std::endl;
  // get the information in the tree being used for the fake rate parameterization
  if(jetProb == NULL) std::cout << "globalProb NOT SET" << std::endl;
  if(debug > 2) std:: cout << "[JetProbMaster] parsing binning variable name" << std::endl;
  std::string	binVar	= jetProb->getBinningVarName(); 
  if(debug > 2) std:: cout << "[JetProbMaster] parsing category name" << std::endl;
  std::string	catVar	= jetProb->getCategoryVarName(); 
  if(debug > 2) std:: cout << "[JetProbMaster] Getting leaf information for binVar and catVar: " << binVar << " " << catVar <<  std::endl;
  TLeaf *	binLeaf	= jetTree->GetLeaf(binVar.c_str());
  TLeaf *	catLeaf	= jetTree->GetLeaf(catVar.c_str());
  int		nJets   = binLeaf->GetNdata();

  // check for  outlier njet scenarios that will be computationally poor 
  if (nJets > 100) {
    std::cout << "ERROR EVENT HAS MORE THAN 100 JETS NJETS: " << nJets << std::endl;
    exit(1);
  }
  if(nJets > 50) {
    std:: cout << "[WARNING] EVENT HAS MORE THAN 50 JETS NJETS: " << nJets << std::endl;
  }

  // setup for parsing from the branch
  std::vector<double> binValues;
  std::vector<double> catValues;

  if(debug > 2) std:: cout << "[JetProbMaster] filling arrays forbinVar and catVar" << std::endl;
  if(debug > 2) std:: cout << "[JetProbMaster] nJets in Event = " << nJets << std::endl;
  // fill an array with the values in the ttree
  for(int ii = 0; ii < nJets; ++ii) {
    binValues.push_back(binLeaf->GetValue(ii));
    catValues.push_back(catLeaf->GetValue(ii));
    if(debug > 2) std:: cout << "[JetProbMaster] binVar " << binValues[ii] << " " << catValues[ii] << std::endl;
  }
  
  // use the arrays to compute the configurational probabilities  
  std::vector<double> nTaggedProbVector;
  // check the tag probability in each scenario
  if(debug > 2) std:: cout << "[JetProbMaster] Looping nTags for event" << std::endl;
  for(int ii = 0; (ii <= nJets) && (ii <= maxJetsTagged); ++ii) {
    double prob = getNJetProbability(&binValues[0], &catValues[0], ii, nJets);
    if(debug > 4) std:: cout << "[JetProbMaster] Adding Probability:" << prob << " to vector " << std::endl;
    nTaggedProbVector.push_back(prob);
  }

  // fill in the rest with blanks
  for(int ii = maxJetsTagged + 1; ii <= nJets; ++ii) {
    nTaggedProbVector.push_back(0);
  }

  // make checks that output is sensible 
  if(debug > 1) { 
    std::cout << " Probability Vector for event: " << std::endl;
    for(int ii = 0; ii <= nJets; ++ii) {     
      std::cout << " (" << ii << "," << nTaggedProbVector[ii] << ") ";
    }
    std::cout << std::endl;      
  }

  return nTaggedProbVector;
}

// helper method for the njet tagged vector. calculates a specific permutation
std::pair<double, double> jetProbabilityMasterComputer::getNJetProbabilityError(double * const binValues, 
									      double * const catValues, 
									      int nJetsTagged, int nJets) {
  if(debug > 4) std::cout << "\n\n------\n[jetProbabilityMasterComputer] Computing Probability Error for nTaggedJets " << 
		  nJetsTagged << " and N jets " << nJets << std::endl;

  double    total_error_up = 0;	// sum of all dp_djet terms with error fluctation up
  double    total_error_dn = 0;	// sum of all dp_djet terms with error fluctation down

  int nConfig = pow2[nJets] - 1;

  // calculate the term corresponding to each jet indexed by i in dP/dp_i
  for(int jet = 0; jet < nJets; ++jet ) {
    double dp_djet_term_up = 0;  // term corresponding to varying the probability for a given jet
    double dp_djet_term_dn = 0;  // term corresponding to varying the probability for a given jet

    // get the errors for the jet in question
    std::pair<double, double>	jetProbErrors = jetProb->getJetFakeProbabilityError(binValues[jet], catValues[jet]);
    double			jetErrUp      = jetProbErrors.first;
    double			jetErrDn      = jetProbErrors.second;

    // calculate each sub-term in the probability (1 for each configuration)
    for(long int ii = 0; ii <= nConfig; ++ii) {
      int ntags = getBinaryDigitSum(ii, nJets);
      if (ntags != nJetsTagged) continue; // only include terms with the correct configuration of tags

      double subTerm = 1;
      // check each bit in the 
      for(int subjet = 0; subjet < nJets; ++subjet) {
	bool isTagged = (ii & pow2[subjet]) > 0;
	if(debug > 5) std::cout << "[jetProbabilityMasterComputer] isTagged =" << isTagged << std::endl;
	// check the sign of the term which is determined by whether the jet is tagged or not 
	if(subjet == jet) {
	  double factor =  isTagged ? 1 : -1;
	  if(debug > 5) std::cout << "[jetProbabilityMasterComputer] subterm factor =" << factor << std::endl;
	  subTerm *= factor;	  
	}
	else { // otherwise multiply by the correct factor
	  double    p_jet  = jetProb->getJetFakeProbability(binValues[subjet], catValues[subjet]);
	  double    factor = isTagged ? p_jet : (1 - p_jet);
	  if(debug > 5) std::cout << "[jetProbabilityMasterComputer] subterm factor =" << factor << std::endl;
	  subTerm *= factor;
	}	
	if(debug > 5) std::cout << "[jetProbabilityMasterComputer] dP_dj subterm index=" << subjet << "progressive subTerm Val=" << subTerm << std::endl;	
	
      } // end loop over subjets in the given configuration
      // add the subterm to the total term in quadrature
      dp_djet_term_up += subTerm;
      dp_djet_term_dn += subTerm;

      if(debug > 5) std::cout << "[jetProbabilityMasterComputer] dp_djet_term_up val=" << dp_djet_term_up << std::endl;
      if(debug > 5) std::cout << "[jetProbabilityMasterComputer] dp_djet_term_dn val=" << dp_djet_term_dn << std::endl;
    } // end loop over all configurations
    
    total_error_up += dp_djet_term_up * dp_djet_term_up * jetErrUp * jetErrUp;
    total_error_dn += dp_djet_term_dn * dp_djet_term_dn * jetErrDn * jetErrDn;    

    if(debug > 5) std::cout << "[jetProbabilityMasterComputer] dp_djet_term=" << jet   <<  " jetErrUp=" << jetErrUp << std::endl;
    if(debug > 5) std::cout << "[jetProbabilityMasterComputer] dp_djet_term=" << jet   <<  " jetErrDn=" << jetErrDn << std::endl;


  } // end loop over the p_i in dP/dp_i  
  
  // apply a sqrt to the total
  total_error_up = std::sqrt(total_error_up);
  total_error_dn = std::sqrt(total_error_dn);

  if(debug > 5) std::cout << "\n\n[jetProbabilityMasterComputer] total_error_up=" << total_error_up << std::endl;
  if(debug > 5) std::cout << "[jetProbabilityMasterComputer] total_error_dn=" << total_error_dn << std::endl;


  std::pair<double, double> error_pair(total_error_up, total_error_dn);

  return error_pair;
} // end get n jet probability errors



// helper method for the njet tagged vector. calculates a specific permutation
double jetProbabilityMasterComputer::getNJetProbability(double * const binValues, 
						       double * const catValues, 
						       int nJetsTagged, int nJets) {
  
  if(debug > 2) std:: cout << "\n\n[JetProbMaster] computing ntags = " << nJetsTagged << " .... with nJets =  " << nJets << std::endl;
  // As a binary variable corresponds to all tagging combinations
  int nConfig = pow2[nJets] - 1;
  
  double probabilitySum = 0;
  if(debug > 2) std:: cout << "\n\n[JetProbMaster] begin loop over configurations" << std::endl;
  // loop over all posible configurations
  for(long int ii = 0; ii <= nConfig; ++ii) {
    // check if the binary representation has the correct number of tags
    if(debug > 9) std:: cout << "[JetProbMaster] computing configuration for " << ii << std::endl;
    if(debug > 9) {
      std:: cout << "[JetProbMaster]  binary representation ";
      for (int jj = nJets - 1; jj >= 0; --jj) {
	int bitVal = (ii & pow2[jj]) > 0;			
	if(bitVal) std::cout << "1";
	else std::cout << "0";
      } 
      std::cout << "\n";
    }
    int ntags = getBinaryDigitSum(ii, nJets);
    if (ntags != nJetsTagged) continue; // irrelavent contribution

    if(debug > 9) { 
      std::cout << "[JetProbMaster] configuration has correct number of tags: " << 
	nJetsTagged << " from calculated value " << ntags << " of possible " << nJets << "in event" << std::endl;
    }

    // this is the term for a specific tag configuration to be added to the total
    double localConfigProb = 1;
    // loop over each binary digit and check its value
    for(int jj = 0; jj < nJets; ++jj) {
      // parse the value of the bit in binary
      bool	bitTagged	       = (ii & pow2[jj]) > 0;
      if(debug > 9) std:: cout << "[JetProbMaster] getting jet fake rate probability... for  " << binValues[jj] << " " << catValues[jj] << std::endl;
      double	tagProbability = jetProb->getJetFakeProbability(binValues[jj], catValues[jj]);

      if (tagProbability > 1 || tagProbability < 0) {
	std:: cout << "[JetProbMaster] Invalid Jet Probability " << tagProbability << "setting to 0" << std::endl;
	tagProbability = 0;
      }

      if(debug > 9) std:: cout << "[JetProbMaster] jet probability # " << jj << " p = " << tagProbability << " isTagged? " << bitTagged << std::endl;;

      // multiply the probability byw whether the jet was tagged or not
      double term       = bitTagged ? tagProbability : (1 - tagProbability);

      localConfigProb *= term;

      if(debug > 9) {
	std:: cout << "[JetProbMaster] local term value = " << term << 
		      " local probability after application " << localConfigProb << std::endl;
      }
    } // end loop over digits of the configuration binary literal

    probabilitySum += localConfigProb;
    if(debug > 9) std:: cout << "[JetProbMaster] probability sum after adding term = " << probabilitySum << std::endl;
  } // end loop over all configurations

  if(debug > 9) std:: cout << "\t\t[JetProbMaster] final njet sum p = " << probabilitySum << " for nJetsTagged = " << nJetsTagged;
  return probabilitySum;

} // end method


// calculates the sum of the binary literal
// njets only provided to determine the highest order of magnitude
long int jetProbabilityMasterComputer::getBinaryDigitSum(long int num, int nJets){
  // we have pre-computed most scenarios
  if (num <= 1024) return hammingWeight[num];

  int sum = 0;
  // loop over the position of each binary digit
  for(int ii = 0; ii < nJets; ++ii) {
    long int value = (num & pow2[ii]) > 0;
    sum += value;
  }

  return sum;  
}


