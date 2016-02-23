
// Joshua Hardenbrook 18-11-2015
// 1. encapsulates the event level selection, as well as the jet variable selection
// 2. parses the variables / cuts from a JSON
#include <json/json.h>
#include <json/json-forwards.h>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLeaf.h"

class jetSelector {
 public:
  explicit jetSelector(const Json::Value & selectorJSON, 
		       const bool & runChop_, 
		       const int & nDivisions_, 
		       const int & probIndex_, 
		       const int & valiIndex_, 
		       const int & debug_);
 
  // prints a summary of all the current requirements on events and jets
  //void printSelectionSummary();
  
  // pull out the branches from the ntuples we want to keep
  TTree* shallowCopyTree(TTree* oldTree);
  // get a vector containing the tagging result
  std::vector<bool> getJetTaggedVector(TTree * tree, int event);
  bool doesEventPassSelection(TTree * tree, long int event);

  // accessors
  // -- for the fake rate histograms binning
  std::vector<double> getHistBinning() const { return histBinVals; }
  std::vector<double> getCatBinning() const { return catBinVals; }
  std::string getJetCutString() const { return jetCutString;}
  std::string getEventCutString() const { return eventCutString;}
  std::string getJetBaselineCutString() const { return baselineJetCutString;}
  std::string getTriggerCutOnlyString() const { return triggerCutOnlyString;}
  // -- for the name of the variable
  std::string getBinningVarName() const { return binningVar;}
  std::string getCategoryVarName() const { return categoryVar;}

  // variables for the tree shallow copy
  Json::Value eventVariablesToSave;
  Json::Value jetVariablesToSave; 

  // selection parameters for the tagging
  Json::Value jetSelection;  
  Json::Value eventSelection;  

  // containers for the values of the fake rate hist bins
  std::vector<double> histBinVals;
  std::vector<double> catBinVals;
  // cuts generated by the JSON for use with the draw command
  std::string jetCutString;
  std::string eventCutString;
  std::string baselineJetCutString;
  std::string triggerCutOnlyString;
  // variables for the fake rate encoded in the JSON
  std::string binningVar;
  std::string categoryVar;

  const bool runChop;
  const int nDivisions, probIndex, valiIndex;
  const int debug;
  bool doEven;
};
