#define CobraClass_cxx
#include "CobraClass.hh"

using namespace std;

CobraClass::CobraClass(TString _fileName) : fChain(0)
{
    cout<<endl<<"--> Input file name: "<<_fileName<<endl;
    fChain = new TChain("SigTree");
    fChain->Add(_fileName.Data());
    Init();
}

CobraClass::~CobraClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CobraClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CobraClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CobraClass::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!fChain) return;

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("amplitude", amplitude, &b_amplitude);
   fChain->SetBranchAddress("sample_rate", &sample_rate, &b_sample_rate);
   fChain->SetBranchAddress("operation_mode", &operation_mode, &b_operation_mode);
   fChain->SetBranchAddress("trigger_depth", &trigger_depth, &b_trigger_depth);
   fChain->SetBranchAddress("trigger_slope", &trigger_slope, &b_trigger_slope);
   fChain->SetBranchAddress("trigger_source", &trigger_source, &b_trigger_source);
   fChain->SetBranchAddress("trigger_level", &trigger_level, &b_trigger_level);
   fChain->SetBranchAddress("sample_depth", &sample_depth, &b_sample_depth);
   fChain->SetBranchAddress("captured_gain_v", &captured_gain_v, &b_captured_gain_v);
   fChain->SetBranchAddress("captured_data", &captured_data, &b_captured_data);
   fChain->SetBranchAddress("current_mem_ptr", &current_mem_ptr, &b_current_mem_ptr);
   fChain->SetBranchAddress("starting_adress", &starting_adress, &b_starting_adress);
   fChain->SetBranchAddress("trigger_adress", &trigger_adress, &b_trigger_adress);
   fChain->SetBranchAddress("ending_adress", &ending_adress, &b_ending_adress);
   fChain->SetBranchAddress("trigger_time_date", &trigger_time_date, &b_trigger_time_date);
   fChain->SetBranchAddress("trigger_coupling", &trigger_coupling, &b_trigger_coupling);
   fChain->SetBranchAddress("trigger_gain_v", &trigger_gain_v, &b_trigger_gain_v);
   fChain->SetBranchAddress("prob", &prob, &b_prob);
   fChain->SetBranchAddress("inverted_data", &inverted_data, &b_inverted_data);
   fChain->SetBranchAddress("board_type", &board_type, &b_board_type);
   fChain->SetBranchAddress("resolution_12_bits", &resolution_12_bits, &b_resolution_12_bits);
   fChain->SetBranchAddress("multiple_record", &multiple_record, &b_multiple_record);
   fChain->SetBranchAddress("trigger_probe", &trigger_probe, &b_trigger_probe);
   fChain->SetBranchAddress("sample_bits", &sample_bits, &b_sample_bits);
   fChain->SetBranchAddress("extended_trigger_time", &extended_trigger_time, &b_extended_trigger_time);
   fChain->SetBranchAddress("imped_a", &imped_a, &b_imped_a);
   fChain->SetBranchAddress("imped_b", &imped_b, &b_imped_b);
   fChain->SetBranchAddress("external_tbs", &external_tbs, &b_external_tbs);
   fChain->SetBranchAddress("external_clock_rate", &external_clock_rate, &b_external_clock_rate);
   fChain->SetBranchAddress("sample_offset", &sample_offset, &b_sample_offset);
   fChain->SetBranchAddress("sample_resolution", &sample_resolution, &b_sample_resolution);
   fChain->SetBranchAddress("dc_offset", &dc_offset, &b_dc_offset);
   fChain->SetBranchAddress("amplitude_f", amplitude_f, &b_amplitude_f);
   fChain->SetBranchAddress("size_of_data", &size_of_data, &b_size_of_data_r);
   fChain->SetBranchAddress("samp_per_msec", &samp_per_msec, &b_samp_per_msec_r);
   fChain->SetBranchAddress("divide", &divide, &b_divide);
   Notify();
}

Bool_t CobraClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CobraClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CobraClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void CobraClass::Loop(TTree *_tree)
{
//   In a ROOT session, you can do:
//      root> .L CobraClass.C
//      root> CobraClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nEntries = fChain->GetEntries();
   cout<<"--> nEntries = "<<nEntries<<endl;

   Double_t* Amplitude = new Double_t[Constants::_nPoints];
   Double_t* MaxAmpl = new Double_t[Constants::_nPeakMax];
   Double_t* TimeMaxAmpl = new Double_t[Constants::_nPeakMax];

   TString  _date;

   fChain->GetEntry(0);

   Double_t _timeCf;
   Double_t _maxAmpl;
   Double_t _unixTime;

   _unixTime = GetUnixTime(trigger_time_date,_date);
   cout<<_date<<endl;
   printf ("--> Unixtime: %20.f [sec]\n", _unixTime);

   _tree->Branch("UnixTime",     &_unixTime,     "_unixTime/D");
   _tree->Branch("TimeCF",       &_timeCf,       "_timeCf/D");
   _tree->Branch("MaxAmpl",      &_maxAmpl,      "_maxAmpl/D");

   for (Long64_t _entryID = 0; _entryID < nEntries; _entryID++)
   {
       cout<<"--> Event: "<<_entryID+1<<"/"<<nEntries<<endl;
       //-------------//
       //-- To plot --//
       //-------------//
       //Long64_t _min_i = 0;
       //Long64_t _max_i = Constants::_nPoints;
       //PlotEvent(_entryID,_min_i,_max_i);
       //-------------//

       TransformAmplitude(_entryID,Amplitude,0,Constants::_nPoints);
       cout<<"--> Number of all peaks with ampl.  > "<<Constants::_level<<" [mV]: "<<GetNumPeaks(Amplitude,0,Constants::_nPoints)<<endl;
       cout<<"--> Number of real peaks with ampl. > "<<Constants::_level<<" [mV]: "<<GetTrueNumPeaks(Amplitude,0,Constants::_nPoints)<<endl;
       Int_t num_peaks = GetMaxAmplPerPeak(Amplitude, 0, Constants::_nPoints, MaxAmpl, TimeMaxAmpl);

       if(num_peaks > 0)
       {
           for(Int_t k = 0; k < num_peaks; k++)
           {
               _unixTime    = GetUnixTime(trigger_time_date,_date);
               _maxAmpl     = MaxAmpl[k];
               _timeCf      = _entryID*Constants::_nPoints*Constants::_dtime + GetTimeCF(Amplitude,MaxAmpl[k],TimeMaxAmpl[k]);

               _tree->Fill();
           }
       }
   }
}

void CobraClass::PlotEvent(Long64_t _entryID, Long64_t _min_i, Long64_t _max_i)
{
    cout<<"--> ----------------- <--"<<endl;
    cout<<"--> Plotting WaveForm <--"<<endl;
    cout<<"--> ----------------- <--"<<endl;
    fChain->GetEntry(_entryID);
    TString _fileName = "plot_";
    _fileName += _entryID;
    _fileName += ".root";
    cout<<"--> Output file with plots: "<<_fileName<<endl;
    cout<<"--> Sample rate: "<<sample_rate/1.0e9<<" GS/s"<<endl;
    cout <<"--> nPoints = "<<Constants::_nPoints<<endl;
    Double_t* _time = new Double_t[Constants::_nPoints];
    Double_t* _ampl = new Double_t[Constants::_nPoints];

    TH1D* h1 = new TH1D("h1","Amplitude",2000,0,2000);
    TH1D* h2 = new TH1D("h2","Peaks Amplitude",2000,0,2000);
    TH1D* h3 = new TH1D("h3","Peaks TimeCF [sec]",1000000,_entryID*0.001,0.001*(_entryID+1));

    if(_min_i < 0)          {_min_i = 0;}
    if(_max_i > Constants::_nPoints)   {_max_i = Constants::_nPoints;}

    for(Long64_t i = _min_i; i < _max_i; i++)
    {
        _time[i] = _entryID*Constants::_nPoints*Constants::_dtime + i*Constants::_dtime;
        _ampl[i] = (-1.0)*amplitude_f[i];
        if(_ampl[i] > 9000.) _ampl[i] = 0.0;

        h1->Fill(_ampl[i]);

        printf("\r--> Entry %llu: %3.1f %%",_entryID,100*(Double_t)i/(_max_i-_min_i));
        fflush(stdout);
    }

    cout<<endl;

    Double_t* max_ampl_per_peak = new Double_t[Constants::_nPeakMax];
    Double_t* time_max_ampl_per_peak = new Double_t[Constants::_nPeakMax];
    Int_t num_peaks = GetMaxAmplPerPeak(_ampl, _min_i, _max_i, max_ampl_per_peak, time_max_ampl_per_peak);

    TFile* _file = new TFile(_fileName.Data(),"RECREATE");
    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    c1->cd();
    TGraph* gr = new TGraph((_max_i-_min_i),_time,_ampl);
    gr->SetName("gr");
    gr->Draw("APL");

    if(num_peaks > 0)
    {
        TLine *line1[num_peaks];
        for(Int_t k = 0; k < num_peaks; k++)
        {
            line1[k] = new TLine(time_max_ampl_per_peak[k],0,time_max_ampl_per_peak[k],2000);
            line1[k]->SetLineWidth(2);
            line1[k]->SetLineStyle(2);
            line1[k]->SetLineColor(kRed);
            line1[k]->Draw("same");

            h2->Fill(max_ampl_per_peak[k]);
            h3->Fill(_entryID*Constants::_nPoints*Constants::_dtime + GetTimeCF(_ampl,max_ampl_per_peak[k],time_max_ampl_per_peak[k]));
        }
    }

    c1->Write();

    h1->Write();
    h2->Write();
    h3->Write();

    _file->Close();

    delete gr;
    delete c1;

    delete [] _time;
    delete [] _ampl;
    delete [] max_ampl_per_peak;
    delete [] time_max_ampl_per_peak;
    delete h1;
    delete h2;
    delete h3;

    cout<<"--> ----------------- <--"<<endl;
}

void CobraClass::TransformAmplitude(Long64_t _entryID, Double_t *OutAmpl, Long64_t _min_i, Long64_t _max_i)
{
    fChain->GetEntry(_entryID);

    for(Long64_t i = _min_i; i < _max_i; i++)
    {
        OutAmpl[i] = (-1.0)*amplitude_f[i];
        if(OutAmpl[i] > 9000.) OutAmpl[i] = 0.0;
    }
}

int CobraClass::GetNumPeaks(Double_t *InAmpl, Long64_t _min_i, Long64_t _max_i)
{
    Int_t num_peaks = 0;

    for (Long64_t i = _min_i; i < _max_i; i++)
    {
        if(InAmpl[i] > Constants::_level)
        {
            num_peaks++;
            while(InAmpl[i] > Constants::_level)
            {
                i++;
            }
            i--;
        }
    }
    //cout<<"--> GetNumPeaks() <--"<<endl;
    return num_peaks;
}

int CobraClass::GetTrueNumPeaks(Double_t *InAmpl, Long64_t _min_i, Long64_t _max_i)
{
    Int_t num_peaks = 0;
    Double_t max = 0.0;
    Double_t time_max = -Constants::_separTime;

    for (Long64_t i = _min_i; i < _max_i; i++)
    {
        while(i*Constants::_dtime < (time_max + Constants::_separTime) && i < _max_i)
        {
            i++;
        }

        if(InAmpl[i] > Constants::_level)
        {
            max = InAmpl[i];
            time_max = i*Constants::_dtime;

            while(InAmpl[i] > Constants::_level && i < _max_i)
            {
                i++;
                if(max < InAmpl[i])
                {
                    max = InAmpl[i];
                    time_max = i*Constants::_dtime;
                }
            }
            i--;
            num_peaks++;
        }
    }
    //cout<<"--> GetTrueNumPeaks() <--"<<endl;
    return num_peaks;
}

int CobraClass::GetMaxAmplPerPeak(Double_t *InAmpl, Long64_t _min_i, Long64_t _max_i, Double_t *OutAmpl, Double_t *OutTime)
{
    Int_t num_peaks = 0;
    Double_t max = 0.0;
    Double_t t = 0;
    for (Int_t i = _min_i; i < _max_i; i++)
    {
        if (InAmpl[i] > Constants::_level)
        {
            while(InAmpl[i] > Constants::_level && i < _max_i)
            {
                if (InAmpl[i] > max)
                {
                    max = InAmpl[i];
                    t = i*Constants::_dtime;
                }
                i++;
            }
            i--;

            OutAmpl[num_peaks] = max;
            OutTime[num_peaks] = t;
            max = 0.0;

            num_peaks++;
        }
    }
    //cout<<"--> GetMaxAmplPerPeak() <--"<<endl;
    return num_peaks;
}

double CobraClass::GetTimeCF(Double_t *InAmpl, Double_t MaxAmpl, Double_t TimeMax)
{
    Long64_t i = TimeMax/Constants::_dtime;
    Double_t OutTime = -999.999;
    for (Long64_t t = i; t > 0; t--)
    {
        if(InAmpl[t] < Constants::_CF*MaxAmpl)
        {
            Double_t x1 = t*Constants::_dtime + Constants::_dtime;
            Double_t x2 = t*Constants::_dtime;
            Double_t y1 = InAmpl[t+1];
            Double_t y2 = InAmpl[t];
            Double_t k = (y1-y2)/(x1-x2);
            Double_t b = (y2*x1-y1*x2)/(x1-x2);
            Double_t y3 = Constants::_CF*MaxAmpl;
            OutTime = (y3-b)/k;
            break;
        }
    }
    //cout<<"--> GetTimeCF() <--"<<endl;
    return OutTime;
}

double CobraClass::GetUnixTime(Int_t _date_in, TString &_date_out)
{
    int year = ((_date_in >> 25) & 0x7f) + 1980;
    int month = (_date_in >> 21) & 0x0f;
    int day = (_date_in >> 16) & 0x1f;
    int hour = (_date_in >> 11) & 0x1f;
    int minute = (_date_in >> 5) & 0x3f;
    int second = (_date_in << 1) & 0x3e;

    _date_out = "--> Date: ";
    _date_out += day;
    _date_out += "/";
    _date_out += month;
    _date_out += "/";
    _date_out += year;
    _date_out += " | Time: ";
    _date_out += hour;
    _date_out += ":";
    _date_out += minute;
    _date_out += ":";
    _date_out += second;

    time_t rawtime;
    struct tm * timeinfo;

    /* get current timeinfo: */
    time ( &rawtime ); //or: rawtime = time(0);
    /* convert to struct: */
    timeinfo = localtime ( &rawtime );

    /* now modify the timeinfo to the given date: */
    timeinfo->tm_year   = year - 1900;
    timeinfo->tm_mon    = month - 1;    //months since January - [0,11]
    timeinfo->tm_mday   = day;          //day of the month - [1,31]
    timeinfo->tm_hour   = hour;         //hours since midnight - [0,23]
    timeinfo->tm_min    = minute;          //minutes after the hour - [0,59]
    timeinfo->tm_sec    = second;          //seconds after the minute - [0,59]

    return (double)timegm ( timeinfo );
}
