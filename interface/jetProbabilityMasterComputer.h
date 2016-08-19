// Joshua Hardenbrook 18-11-2015
// -- taking as input the globalJetProbabilities, this class computes the combinitorial possibilities of Ntags existing in a single event
// -- compute likelihood discriminants from the globalJetProbabilities 
#include "globalJetProbabilities.h"
#include "TTree.h"
#include "TLeaf.h"

class jetProbabilityMasterComputer {
 public:
  // send a pointer to the tree containing the branches to use for the fake rate computation
  explicit jetProbabilityMasterComputer(globalJetProbabilities * jetProb,  TTree *  jetTree, const int & debug);
  ~jetProbabilityMasterComputer();

  // calculate the probability of N tags as a vector from 0 to N jets in the event
  // limit the calculator to maxJetsTagged
  const std::vector<double>                   getNJetTaggedVector(const long int & eventNumber,const int & maxJetsTagged);
  const std::vector<double>                   getJetProbabilityVector(const long int & eventNumber);
  const std::vector<std::pair<double, double>> getNJetErrorVector(const long int & eventNumber,const int & maxJetsTagged);

  // helper method for the njet tagged vector. calculates a specific permutation with nJets in the event
  // and nJets being tagged. pass the variable as well as the category for the jets in the event
  double getNJetProbability(const std::vector<double> binValues, const std::vector<double> catValues, const int & nJetsTagged, const int & nJets);  
  const std::pair<double, double> getNJetProbabilityError(const std::vector<double> binValues, const std::vector<double> catValues, 
						  const int & nJetsTagged, const int & nJets);

  const int getNTaggedInConfig(const int & nJets, const std::vector<double> & binVals, const int & config, const double & binVal);
  const int binValMultiplicity(const std::vector<double> & binVals, const double & binVal) ;
  
  const std::vector<std::pair<double,double>> getNTrackTermsForNTagError(const long int & eventNumber, const int & nTagged);
  // calculates the sum of the digits of the binary literal (corresponding to the number of tags)
  const int getBinaryDigitSum(const long int & num, const int & njets);

  // precompute powers of 2 
  std::vector<long int> pow2 = { 1 , 2 , 4 , 8 , 16 , 32 , 64 , 128 , 256 , 512 , 1024 , 2048 , 4096 , 8192 , 16384 , 32768 , 65536 , 131072 , 262144 , 524288 , 1048576};			 // , 2097152 , 4194304 , 8388608 , 16777216 , 33554432 , 67108864 , 134217728 , 268435456 , 536870912 , 1073741824 , 2147483648 , 4294967296 , 8589934592 , 17179869184 , 34359738368 , 68719476736 , 137438953472 , 274877906944 , 549755813888 , 1099511627776 , 2199023255552 , 4398046511104 , 8796093022208 , 17592186044416 , 35184372088832 , 70368744177664 , 140737488355328 , 281474976710656 , 562949953421312 , 1125899906842624 , 2251799813685248 , 4503599627370496 , 9007199254740992 , 18014398509481984 , 36028797018963968 , 72057594037927936 , 144115188075855872 , 288230376151711744 , 576460752303423488 , 1152921504606846976 , 2305843009213693952 , 4611686018427387904 , 9223372036854775808 , 18446744073709551616 , 36893488147419103232 , 73786976294838206464 , 147573952589676412928 , 295147905179352825856 , 590295810358705651712 , 1180591620717411303424 , 2361183241434822606848 , 4722366482869645213696 , 9444732965739290427392 , 18889465931478580854784 , 37778931862957161709568 , 75557863725914323419136 , 151115727451828646838272 , 302231454903657293676544 , 604462909807314587353088 , 1208925819614629174706176 , 2417851639229258349412352 , 4835703278458516698824704 , 9671406556917033397649408 , 19342813113834066795298816 , 38685626227668133590597632 , 77371252455336267181195264 , 154742504910672534362390528 , 309485009821345068724781056 , 618970019642690137449562112 , 1237940039285380274899124224 , 2475880078570760549798248448 , 4951760157141521099596496896 , 9903520314283042199192993792 , 19807040628566084398385987584 , 39614081257132168796771975168 , 79228162514264337593543950336 , 158456325028528675187087900672 , 316912650057057350374175801344 , 633825300114114700748351602688 };

  //  int * hammingWeight = new int[1025];
  const std::vector<int> hammingWeight = {0, 1 , 1 , 2 , 1 , 2 , 2 , 3 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 6 , 7 , 7 , 8 , 7 , 8 , 8 , 9 , 7 , 8 , 8 , 9 , 8 , 9 , 9 , 10 }; 

  globalJetProbabilities*   jetProb;
  TTree *		    jetTree;
  int     		    debug;
};
