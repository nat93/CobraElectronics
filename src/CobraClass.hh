//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  7 10:23:57 2017 by ROOT version 6.10/02
// from TTree SigTree/Signal Data
//////////////////////////////////////////////////////////

#ifndef CobraClass_hh
#define CobraClass_hh

// C++
#include "iostream"
#include "string"
#include "fstream"
#include "vector"
#include "stdlib.h"
#include "time.h"
#include "iomanip"
#include "assert.h"

// ROOT
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSystem.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph2DErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TRandom3.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TSpline.h"

// MY
#include "./include/Constants.hh"

// Header file for the classes stored in the TTree if any.

class CobraClass {
public :
   TChain         *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           amplitude[Constants::_nPoints];
   Int_t           sample_rate;
   Short_t         operation_mode;
   Int_t           trigger_depth;
   Short_t         trigger_slope;
   Short_t         trigger_source;
   Short_t         trigger_level;
   Int_t           sample_depth;
   Double_t        captured_gain_v;
   Short_t         captured_data;
   Int_t           current_mem_ptr;
   Int_t           starting_adress;
   Int_t           trigger_adress;
   Int_t           ending_adress;
   Int_t           trigger_time_date;
   Short_t         trigger_coupling;
   Double_t        trigger_gain_v;
   Short_t         prob;
   Short_t         inverted_data;
   Short_t         board_type;
   Short_t         resolution_12_bits;
   Short_t         multiple_record;
   Short_t         trigger_probe;
   Short_t         sample_bits;
   Int_t           extended_trigger_time;
   Short_t         imped_a;
   Short_t         imped_b;
   Float_t         external_tbs;
   Float_t         external_clock_rate;
   Short_t         sample_offset;
   Short_t         sample_resolution;
   Short_t         dc_offset;
   Float_t         amplitude_f[Constants::_nPoints];
   Int_t           size_of_data;
   Int_t           samp_per_msec;
   Double_t        divide;

   // List of branches
   TBranch        *b_amplitude;   //!
   TBranch        *b_sample_rate;   //!
   TBranch        *b_operation_mode;   //!
   TBranch        *b_trigger_depth;   //!
   TBranch        *b_trigger_slope;   //!
   TBranch        *b_trigger_source;   //!
   TBranch        *b_trigger_level;   //!
   TBranch        *b_sample_depth;   //!
   TBranch        *b_captured_gain_v;   //!
   TBranch        *b_captured_data;   //!
   TBranch        *b_current_mem_ptr;   //!
   TBranch        *b_starting_adress;   //!
   TBranch        *b_trigger_adress;   //!
   TBranch        *b_ending_adress;   //!
   TBranch        *b_trigger_time_date;   //!
   TBranch        *b_trigger_coupling;   //!
   TBranch        *b_trigger_gain_v;   //!
   TBranch        *b_prob;   //!
   TBranch        *b_inverted_data;   //!
   TBranch        *b_board_type;   //!
   TBranch        *b_resolution_12_bits;   //!
   TBranch        *b_multiple_record;   //!
   TBranch        *b_trigger_probe;   //!
   TBranch        *b_sample_bits;   //!
   TBranch        *b_extended_trigger_time;   //!
   TBranch        *b_imped_a;   //!
   TBranch        *b_imped_b;   //!
   TBranch        *b_external_tbs;   //!
   TBranch        *b_external_clock_rate;   //!
   TBranch        *b_sample_offset;   //!
   TBranch        *b_sample_resolution;   //!
   TBranch        *b_dc_offset;   //!
   TBranch        *b_amplitude_f;   //!
   TBranch        *b_size_of_data_r;   //!
   TBranch        *b_samp_per_msec_r;   //!
   TBranch        *b_divide;   //!

   CobraClass(TString _fileName);
   virtual ~CobraClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop(TTree *_tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //Common
   virtual void     PlotEvent(Long64_t _entryID, Long64_t _min_i, Long64_t _max_i);
   virtual void     TransformAmplitude(Long64_t _entryID, Double_t *OutAmpl, Long64_t _min_i, Long64_t _max_i);
   virtual int      GetNumPeaks(Double_t *InAmpl, Long64_t _min_i, Long64_t _max_i);
   virtual int      GetTrueNumPeaks(Double_t *InAmpl, Long64_t _min_i, Long64_t _max_i);
   virtual int      GetMaxAmplPerPeak(Double_t *InAmpl, Long64_t _min_i, Long64_t _max_i, Double_t *OutAmpl, Double_t *OutTime);
   virtual double   GetTimeCF(Double_t *InAmpl, Double_t MaxAmpl, Double_t TimeMax);
   virtual double   GetUnixTime(Int_t _date_in, TString &_date_out);
};

#endif
