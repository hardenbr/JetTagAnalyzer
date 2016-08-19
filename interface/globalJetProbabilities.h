// Joshua Hardenbrook 18-11-2015
// 1. encapsulates the relative tag efficiency given some selection
// 2. build the categories (or histograms) to perform the division for fake rate calculation
// 3. holds the binned pdf for signal/background allowing for the calculation of correlated likelihood variables 
// 3. ouputs an xml/JSON containing 'fixed' fake-rates to be used in the probability calculation
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "jetSelector.h"

class globalJetProbabilities {
 public:
  // constructor for a single probability object
  explicit globalJetProbabilities(const std::string &label, 
				  const std::string &stack,
				  const bool & isMC_ , 
				  const bool & isSig_ , 
				  const double& evWeight_, 
				  const double& xsec_, 
				  TTree*& tree,
				  jetSelector& jetSel,
				  const std::string& divideString_,
				  const int & debug_
				  );
  
  // constructor from a stored json 

  explicit globalJetProbabilities(const std::string& label, 
				  const std::string& stack,
				  const bool & isMC_ , 
				  const bool & isSig_ , 
				  const double& evWeight_, 
				  const double& xsec_, 
				  Json::Value probabilities,
				  const std::string& divideString_,
				  const int & debug_
				  );

  // add jets to the probabilities
  // void addJet(const jetCandidate&);  

  // when all the jets have been added, perform the division for the efficiency
  //  void generateGlobalProbabilities();

  // create the pdfs for the likelihood
  //void generateLikelihoodPDF();

  // called by jetTagAnalyzer to retrieve probabilities to be included in the outputted JSON
  //void getGlobalProbabilities();
  //float getNJetTagProbabilityVector();

  // calculates the local pdf likelihood
  //float getJetLikelihood(const jetCandidate&);

  double getJetFakeProbability(float binVariable, float catVar);
  const std::pair<double, double> getJetFakeProbabilityError(const float binVariable, const float catVar);
  void printHistStatus();
  void addSignalContamination(TTree*& tree, jetSelector & jetSel, float norm);
  void removeSignalRegion(TTree*& tree, jetSelector & jetSel, bool isContam, float norm);


  Json::Value getProbabilitiesJSON();

  // acessor
  std::string getBinningVarName() { return binningVar; }
  std::string getCategoryVarName() { return categoryVar; }
  std::vector<double> getBinningVector() { return histBinVals; }

  // accessors
  TGraphAsymmErrors getRatioGraph();
  TH1D		    getTaggedHist();
  TH1D		    getAllHist();
  TH1D              getCentralEffHist();
  TH1D              getCentralEffHistErrUp();
  TH1D              getCentralEffHistErrDn();

  // local variables
  TH1D		    taggedJetHist;
  TH1D		    allJetHist;
  TGraphAsymmErrors ratioGraph;
  // histograms of the assymetric errors
  TH1D              ratioHistEff;
  TH1D              ratioHistEffErrUp;
  TH1D              ratioHistEffErrDn;

  std::vector<double>   histBinVals;
  std::vector<double>   catBinVals;

  // sample label
  std::string   label;
  // bins and category variables
  std::string   binningVar;
  std::string   categoryVar;
  // selection strings
  std::string	jetCutString;
  std::string	eventCutString;
  std::string	baselineJetCutString;
  std::string	triggerCutOnlyString;

  // names for the histograms
  std::string	taggedJetHistName;
  std::string	allJetHistName;
  std::string	effHistName;
  std::string	effHistNameUp;
  std::string	effHistNameDn;

  // string for the division calculation
  //const std::string divide_string = "b(1,1) mode cl=.683";
  //const std::string divide_string = "fc  cl=.683";
  


  // configuration params
  int nBins;
  const bool isMC;
  const bool isSig;
  double evWeight;
  double xsec;
  const std::string divideString;
  const int debug;

};
