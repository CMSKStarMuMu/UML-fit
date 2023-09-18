#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore -lMinuit
ROOTCINT  := $(shell which rootcint)

#directories
SOURCEDIR   := ./src
INCLUDEDIR  := ./interface

#exe_files
EXECUTABLE0 := simfit_recoMC_singleComponent
EXECUTABLE1 := simfit4d_recoMC_singleComponent
EXECUTABLE2 := simfit_recoMC_fullAngular
EXECUTABLE3 := simfit_genMC
EXECUTABLE4 := simfit_genMC_multiFit
EXECUTABLE5 := plotMultiResults
EXECUTABLE6 := simfit_toy_fullAngular
EXECUTABLE7 := simfit_recoMC_fullAngularMass
EXECUTABLE8 := simfit_recoMC_fullMass
EXECUTABLE9 := simfit_recoMC_fullAngularMass_toybkg
EXECUTABLE10 := simfit_data_fullAngularMass
EXECUTABLE11 := simfit_data_fullAngularMass_Swave
EXECUTABLE12 := plot_simfit_data_fullAngularMass_Swave
EXECUTABLE13 := plot_simfit_recoMC_fullAngularMass_toybkg
EXECUTABLE14 := plot_simfit_recoMC_fullAngularMass
EXECUTABLE15 := simfit_recoMC_fullAngularMass_toybkg_toyeff
EXECUTABLE16 := simfit_data_fullAngularMass_Swave_toyeff

EXTRACLASS := RooDataHist.cxx
CLASS0     := PdfRT
CLASS1     := PdfWT
CLASS2     := DecayRate
CLASS3     := PdfSigAng
CLASS4     := RooDoubleCBFast
CLASS5     := BoundCheck
CLASS6     := Penalty
CLASS7     := BoundDist
CLASS8     := PdfSigAngMass
CLASS9     := PdfSigMass
CLASS10    := ShapeSigAng
CLASS11    := Fitter
CLASS12    := RooBernsteinSideband

CLASSDICT  := AngDict
CLASSDICT2 := RooDoubleCBDict

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(SOURCEDIR)/$(CLASS0).cc $(SOURCEDIR)/$(CLASS1).cc $(SOURCEDIR)/$(CLASS2).cc $(SOURCEDIR)/$(CLASS3).cc \
        $(SOURCEDIR)/$(CLASS5).cc $(SOURCEDIR)/$(CLASS6).cc $(SOURCEDIR)/$(CLASS7).cc $(SOURCEDIR)/$(CLASS8).cc \
        $(SOURCEDIR)/$(CLASS9).cc $(SOURCEDIR)/$(CLASS10).cc $(SOURCEDIR)/$(CLASS11).cc $(SOURCEDIR)/$(CLASS12).cxx $(CLASSDICT).cc $(SOURCEDIR)/$(EXTRACLASS)

$(CLASSDICT): $(INCLUDEDIR)/$(CLASS0).h $(INCLUDEDIR)/$(CLASS1).h $(INCLUDEDIR)/$(CLASS2).h $(INCLUDEDIR)/$(CLASS3).h \
              $(INCLUDEDIR)/$(CLASS5).h $(INCLUDEDIR)/$(CLASS6).h $(INCLUDEDIR)/$(CLASS7).h $(INCLUDEDIR)/$(CLASS8).h \
              $(INCLUDEDIR)/$(CLASS9).h $(INCLUDEDIR)/$(CLASS10).h $(INCLUDEDIR)/$(CLASS11).h $(INCLUDEDIR)/$(CLASS12).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^

$(CLASSDICT2): $(INCLUDEDIR)/$(CLASS4).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^ -I./vdt	

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE1): $(EXECUTABLE1).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE2): $(EXECUTABLE2).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE3): $(EXECUTABLE3).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE4): $(EXECUTABLE4).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE5): $(EXECUTABLE5).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE6): $(EXECUTABLE6).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE7): $(EXECUTABLE7).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE8): $(EXECUTABLE8).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE9): $(EXECUTABLE9).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE10): $(EXECUTABLE10).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE11): $(EXECUTABLE11).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE12): $(EXECUTABLE12).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE13): $(EXECUTABLE13).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE14): $(EXECUTABLE14).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE15): $(EXECUTABLE15).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

$(EXECUTABLE16): $(EXECUTABLE16).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR) 

#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE0) $(EXECUTABLE1) $(EXECUTABLE7)
