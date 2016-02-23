// Joshua Hardenbrook 18-11-2015
// description: main analysis class for the tag analysis and n-tag jet estimations. 

// Helper Classes Used
// -------------------
// 1. Computes the globalJetProbabilities object (fake rate / sig efficiency) 
// 2. Utilizes a jetSelector to apply jet (tag definition) and event level selection
// 3. Uses the probabilityMasterComputer to perform fast combinatorics of jet probabilities (N_jets choose N_tagged) 

// Actions Performed
//------------------
// 1. Parses sample configuration of ntuples containing trees with the relevant jet information
// 2. Samples are processed individually
// 3. Computes the globalJetProbabilities if not provided as input (sum over samples)
// 4. Uses the globalProbabilities to assign probabilities to jets in the sample
// 5. Uses the globalProbabilities with the jetProbabilitiyMasterComputer to build
// per event expectations for the number of tags.
// 6. Organize output

// TODO 
// ------------------
// 1. JES calculations
// 2. PDF systematics
// 3. Should the limit calculation be included here?

// ANALYZER INPUT
// --------------
// 1. List of samples to process (mc/data, signal, stack declaration)
// 2. Definition of variables to be used for tagging
// 3. Definition of variables to be used for fake-rate parameterization (if necessary)
// 4. globalProbabilities Object (if not provided)
// 5. Which stack is going to be used for the probabilities object? Or always do the full some of the samples?

// ANALYZER OUTPUT
// ---------------
// 1. Outputs 1 tree per sample containing the relevent jet tagging variables, jet prob. vars, 
//    in the parameterization of the fake rate, as well as the jet probabilities (p(fake), p(event has 1 tags), p(event has 2 tag) 
// 2. Output 1 tree with sum over stack with the appropriate event weight included in the tree
// 3. Output histograms (1 per sample) (1 including sum total) containing the counts for the specified jet tag definition
// 4. Outputs a fixed fake rate object (in xml/JSON) that can serve as input to jet probabilitiy calculations 

// OPTIONS
// -------
// 1. Only output jet tree events with N tags
// 2. Only output jets that are (tagged/not-tagged)

class jetTagAnalyzer {
 public:
  public jetTagAnalyzer();

  // SCHEDULING / SETUP
  // parsing of configurations 
  void parseSampleJSON();
  void parseJetSelection();
  // optional parsing of fixed probabilities object 
  void parseGlobalProbabilities();
  
  // SETUP
  // -----
  // make output files, trees, branches, branch addresses
  void setupOutput();
  // make histgrams for probabilities
  void setupProbabilities();

  // TO BE CALLED BY PROCESSING
  // ------------------------------
  // dump the necessary variables from the ntuples
  void fillJetVariables(); 
  // dump event specific variables
  void fillEventVariables(); 
  // call in the processTotal method to build globalProbabilities object
  void computeGlobalProbabilities(); 
  // use the global probabilities to determine the 
  void applyGlobalProbabilities();
  
  // PROCESSING
  // ----------
  void processSample();
  void processStack();
  void processTotal(); // merge the necessary samples?

  // OUTPUT
  // ------
  // output 1 to 1 the files inputted
  void outputSampleFiles();
  // output the tree as a combination of the 
  void outputStackFiles();
  // write the file containing the fixed derivation of the probabilities
  // include the selection used to derive the 
  void outputProbabilities(); 
  // write just the histograms (ntags, baseline variables?)
  void outputHistograms();

  static const int MAX_JETS = 50;

  int nJets;
  double * nJetTagProb[MAX_JETS];
  TH1D * pileupWeightHist;
}
