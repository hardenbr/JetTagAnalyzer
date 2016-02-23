all: jetTagAnalyzer

jetTagAnalyzer: jetTagAnalyzer.cc jetSelector.o globalJetProbabilities.o jetProbabilityMasterComputer.o
	g++ jetTagAnalyzer.cc -o jetTagAnalyzer jsoncpp.o jetSelector.o globalJetProbabilities.o jetProbabilityMasterComputer.o -I ../interface/ `root-config --cflags --glibs` --std=c++11	

jetProbabilityMasterComputer.o: globalJetProbabilities.cc jetSelector.o jsoncpp.o jetProbabilityMasterComputer.cc
	g++ jetProbabilityMasterComputer.cc globalJetProbabilities.o jsoncpp.o -c -I ../interface/ `root-config --cflags --glibs` --std=c++11

globalJetProbabilities.o : globalJetProbabilities.cc jetSelector.o jsoncpp.o
	g++ globalJetProbabilities.cc jetSelector.o jsoncpp.o -c -I ../interface/ `root-config --cflags --glibs` --std=c++11

jetSelector.o: jsoncpp.o jetSelector.cc
	g++ jetSelector.cc jsoncpp.o -c -I ../interface/ --std=c++11 `root-config --cflags --glibs` 

jsoncpp.o: jsoncpp.cpp
	g++ jsoncpp.cpp -c -I ../interface/  -std=c++11

clean:
	rm *.o jetTagAnalyzer