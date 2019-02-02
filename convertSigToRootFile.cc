//root headers
#include <TH1D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

//C, C++ headers
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <float.h>
#include <stdint.h>
#include <limits.h>

#pragma warning(disable: 4996)

using namespace std;

//const int size_of_data_max = 1073741568; //max number of samples
//const int size_of_data_max = 2147483392;
const long size_of_data_max = 4008192000;

const int samp_per_msec = 1000000;       //number of samples per milisecond
const int n_sampleRateVal = 48;          //size of sample_rate_table[] array
const int n_board_type = 22;         //size of board_type_int[] and board_type_char[] arrays 
void fill_sample_rate_table(long int srtab[n_sampleRateVal]);
void fill_board_type_int(int btint[n_board_type]);
void fill_board_type_char(TString btchar[n_board_type]);
//function that converts .sig to .root
void convertSigToRootFile(TString inputDataFileList, TString inputDataFileFolder, TString outputRootFileFolderFilder);

int main(int argc, char* argv[]){
  TString inputDataFileList; // file with list of .sig files
  TString inputDataFileFolder; // folder with .sig files
  TString outputRootFileFolder; // folder with .root files
  if(argc == 5 && atoi(argv[1])==0){
    inputDataFileList = argv[2];
    inputDataFileFolder = argv[3];
    outputRootFileFolder = argv[4];
    cout<<"Input data file list : "<<inputDataFileList<<endl
        <<"Input data path      : "<<inputDataFileFolder<<endl
        <<"Output root path     : "<<outputRootFileFolder<<endl;
    convertSigToRootFile(inputDataFileList, inputDataFileFolder, outputRootFileFolder);
  } // if all parameters are inputted
  else{
    cout<<" ERROR --->  in input arguments "<<endl
        <<"  runID [1] = 0  "<<endl
        <<"        [2] - input data file list"<<endl
        <<"        [3] - input data path"<<endl
        <<"        [4] - output data path"<<endl;
  }  // if don't
  return 0;
}

void convertSigToRootFile(TString inputDataFileList, TString inputDataFileFolder, TString outputRootFileFolder){

  ifstream inputdataList; // stream that will contain the information from inputDataFileList
  inputdataList.open(inputDataFileList.Data()); // opening stream  
  assert(inputdataList.is_open()); // assertion if not opened
  
  string nameOfDataFile; // will contain the names of .sig files
  
  //header variables
  char file_version[14];     // file version
  char ch_name[9];           // channel name
  char comment[256];         // comments
  short control_z;           // artificial end of file
  short sample_rate_index;   // index to the sample rate table
  short operation_mode;      // 1 = single channel, 2 = dual channel
  int trigger_depth;         // number of samples after the trigger point
  short trigger_slope;       // 1 = positive slope, 2 = negative slope
  short trigger_source;      // 1 = chan A, 2 = chan B, 3 = external, 4 = automatic, 5 = keyboard
  short trigger_level;       // stored as an int, actually a byte with the same format as the data
  int sample_depth;          // number of samples stored in the signal section of the file
  short captured_gain;       // index to the input range table
  short captured_data;       // stored as an int, actually a byte with the same format as the data
  int current_mem_ptr;       // where display started when signal was stored
  int starting_adress;       // the first point in the data
  int trigger_adress;        // the point in the data where trigger occurred
  int ending_adress;         // the last point of the captured data
  int trigger_time_date;     // the time and date when the trigger event occurred
  short trigger_coupling;    // 1 = DC, 2 = AC; for the external trigger input
  short trigger_gain;        // index to the input range table
  short prob;                // index to the probe table
  short inverted_data;       // 0 = normal data, 1 = inverted data (CS220), 2 = inverted and flipped data (CS220)
  short board_type;          // the CompuScope board type on which the saved data was captured
  short resolution_12_bits;  // 0 = 8 bit file format, 1 = 12/16 bit file format
  short multiple_record;     // the mode that the saved data was captured in: 0 = normal mode, 
                             // 1 = Hardware multiple record, 2 = Software multiple record
  short trigger_probe;       // index to the probe table
  short sample_offset;       // Used to offset the data for display andconversion to real voltages. 
                             // Normally 127 for 8-bit CompuScopes and -1 for 12-bit CompuScopes. 
                             // File versions before GS V. 2.95 will have sample offset=0 for 12-bit 
                             // CompuScopes and 128 for 8-bit CompuScopes.
  short sample_resolution;   // Used to scale the data for display and conversion to real voltages. 
                             // Normally 128 for 8-bit CompuScopes and 2048 for 12-bit CompuScopes.
  short sample_bits;         // Number of bits in the sampled data. Normally 8 for 8-bit CompuScopes and 12 for 12-bit CompuScopes.
  int extended_trigger_time; // The time when trigger event occurred
  short imped_a;             // Impedance for Channel A. 0 = 1MegaOhm, 0x10 = 50 Ohm.
  short imped_b;             // Impedance for Channel B. 0 = 1MegaOhm, 0x10 = 50 Ohm.
  float external_tbs;        // Time between samples in nanoseconds when using external clock.
  float external_clock_rate; // Minimum sample rate when using external clock.
  short dc_offset;
  char unknwn[256];          // stores trash
  //end of header variables

  // sample rate table
  long int sample_rate_table[n_sampleRateVal];
  fill_sample_rate_table(sample_rate_table);
  //for(int i = 0; i < n_sampleRateVal; i++){
  //cout<<"sample_rate_table["<<i<<"] = "<<sample_rate_table[i]<<endl;
  //}
  //cout<<LLONG_MIN<<endl;
  //cout<<LLONG_MAX<<endl;

  // board type values
  int board_type_int[n_board_type];
  fill_board_type_int(board_type_int);
  //for(int i = 0; i < n_board_type; i++){
  //cout<<std::dec<<"board_type_int["<<i<<"] = "<<std::hex<<board_type_int[i]<<endl;
  //}

  // board types
  TString board_type_char[n_board_type];
  fill_board_type_char(board_type_char);
  //for(int i = 0; i < n_board_type; i++){
  //cout<<"board_type_char["<<i<<"] = "<<board_type_char[i]<<endl;
  //}

  double input_range[7] = {10000, 5000, 2000, 1000, 500, 200, 100}; // input range table
  string probe_table[8] = {"x1", "x10", "x20", "x50", "x100", "x200", "x500", "x1000"};
  
  TString inputDataFile;

  while(inputdataList >> nameOfDataFile){
		
    FILE *pFile;
    
    inputDataFile = inputDataFileFolder + nameOfDataFile;
    pFile = fopen(inputDataFile.Data(), "rb");
    cout << endl;
    cout << "--------------------------------------------------------" << endl << endl;
    cout << "Reading the file     : " << inputDataFile.Data() << endl;
    
    fseek(pFile, 0, SEEK_END);
    
    //Create ROOT file
    TFile* hfile = new TFile(outputRootFileFolder + nameOfDataFile + ".root", "RECREATE", "SigData", 1);
    if(hfile->IsZombie()){
      cout << "PROBLEM with the initialization of the output ROOT ntuple "
	   << outputRootFileFolder << ": check that the path is correct!!!"
	   << endl;
      assert(0);
    }
    
    //Create ROOT tree
    TTree *tree = new TTree("SigTree", "Signal Data");
    hfile->SetCompressionLevel(2);
    tree->SetAutoSave(1000000);
    
    TTree::SetBranchStyle(0);
    
    int size_of_data = ftell(pFile) - 512;
    cout<<"size_of_data: "<<size_of_data<<endl;
    cout<<"size_of_data_max: "<<size_of_data_max<<endl;
    //size_of_data=100000;
    assert(size_of_data <= size_of_data_max);
    
    rewind(pFile);
    
    Int_t amplitude[samp_per_msec];
    Float_t amplitude_f[samp_per_msec];
    TString leaf_name = "amplitude[";
    leaf_name += samp_per_msec;
    leaf_name += "]/I";
        
    tree->Branch("amplitude", amplitude, leaf_name.Data());
    
    unsigned char num_data_h;
    unsigned long long currentPos = 0;
    cout << endl << "Reading of header:" << endl;
    
    ofstream header_info;
    header_info.open(outputRootFileFolder + nameOfDataFile + "_header_info.txt");
    
    //Read header
    currentPos = fread(&file_version,1 , 14, pFile);
    for(int i = 0; i < 14; i++){
      if(file_version[i] == 0x00){break;}
      else{header_info << file_version[i];}
    }
    currentPos += fread(&unknwn, 1, 2, pFile); header_info << endl;
    currentPos += fread(&ch_name, 1, 9, pFile);
    for(int i = 0; i < 9; i++){
      if(ch_name[i] == 0x00){break;}
      else{header_info << ch_name[i];}
    }
    char chn[] = {ch_name[3], ch_name[4]};
    int ch_number = atoi(chn);
    
    currentPos += fread(&unknwn, 1, 2, pFile); header_info << endl;
    currentPos += fread(&comment, 1, 256, pFile); 
    for(int i = 0; i < 256; i++){
      if(comment[i] == 0x00){break;}
      else{header_info << comment[i];}
    }
    currentPos += fread(&unknwn, 1, 2, pFile); header_info << endl; 
    currentPos += fread(&control_z, 1, 2, pFile); header_info << "control_z = " << control_z << " - artificial and of file" << endl; 
    //tree->Branch("control_z", &control_z, "control_z/S");
    currentPos += fread(&sample_rate_index, 1, 2, pFile); sample_rate_table[47] = 1e-09/external_tbs;
    Int_t sample_rate;
    if(sample_rate_index == 47){header_info << "---For time between samples see data from external clock---" << endl;}
    else
      {
	header_info << sample_rate_table[sample_rate_index] << " - sample rate" << endl;
	sample_rate = sample_rate_table[sample_rate_index];
	tree->Branch("sample_rate", &sample_rate, "sample_rate/I");
      }
    
    currentPos += fread(&operation_mode, 1, 2, pFile);
    switch (operation_mode)
      {
      case 1:
	header_info << "Operation mode: single channel" << endl;
	break;
      case 2:
	header_info << "Operation mode: dual channel" << endl;
	break;
      }
    tree->Branch("operation_mode", &operation_mode, "operation_mode/S");
    
    currentPos += fread(&trigger_depth, 1, 4, pFile);  header_info << "Number of samples after the trigger point: " << trigger_depth << endl;
    tree->Branch("trigger_depth", &trigger_depth, "trigger_depth/I");
    currentPos += fread(&trigger_slope, 1, 2, pFile);
    switch (trigger_slope)
      {
      case 1:
	header_info << "Trigger slope: positive slope" << endl;
	break;
      case 2:
	header_info << "Trigger slope: negative slope" << endl;
	break;
      }
    tree->Branch("trigger_slope", &trigger_slope, "trigger_slope/S");
    currentPos += fread(&trigger_source, 1, 2, pFile);
    header_info << "Trigger source: ";
    switch(trigger_source)
      {
      case 1:
	header_info << "channel A";
	break;
      case 2:
	header_info << "channel B";
	break;
      case 3:
	header_info << "external";
	break;
      case 4:
	header_info << "automatic";
	break;
      case 5:
	header_info << "keyboard";
	break;
      }
    header_info << endl;
    
    tree->Branch("trigger_source", &trigger_source, "trigger_source/S");
    
    currentPos += fread(&trigger_level, 1, 2, pFile); 
    header_info << "Trigger level:" << trigger_level*2000./255 - 1000. << " mV" << endl;
    tree->Branch("trigger_level", &trigger_level, "trigger_level/S");
    
    currentPos += fread(&sample_depth, 1, 4, pFile);
    header_info << "Number of samples stored in the signal section of the file: " << sample_depth << endl;
    tree->Branch("sample_depth", &sample_depth, "sample_depth/I");
    
    currentPos += fread(&captured_gain, 1, 2, pFile);
    header_info << "Input range (captured gain): " << "+/-" << input_range[captured_gain] << "mv" << endl;
    Float_t range = 2*input_range[captured_gain];
    double captured_gain_v = 2*input_range[captured_gain];
    tree->Branch("captured_gain_v", &captured_gain_v, "captured_gain_v/D");
    
    currentPos += fread(&captured_data, 1, 2, pFile); 
    header_info << "Captured data: " << -captured_data*2000./255 + 1000. << " mV" << endl;
    tree->Branch("captured_data", &captured_data, "captured_data/S");
    
    currentPos += fread(&current_mem_ptr, 1, 4, pFile); 
    header_info << "Where display started when signal was stored: " << current_mem_ptr << endl;
    tree->Branch("current_mem_ptr", &current_mem_ptr, "current_mem_ptr/I");
    
    currentPos += fread(&starting_adress, 1, 4, pFile); 
    header_info << "The first point in the data (starting_adress): " << starting_adress << endl;
    tree->Branch("starting_adress", &starting_adress, "starting_adress/I");
    
    currentPos += fread(&trigger_adress, 1, 4, pFile); 
    header_info << "The point in the data where trigger occurred (trigger_adress): " << trigger_adress << endl;
    tree->Branch("trigger_adress", &trigger_adress, "trigger_adress/I");
    
    currentPos += fread(&ending_adress, 1, 4, pFile); 
    header_info << "The last point of the captured data (ending_adress): " << ending_adress << endl;
    tree->Branch("ending_adress", &ending_adress, "ending_adress/I");
    
    currentPos += fread(&trigger_time_date, 1, 4, pFile); 
    tree->Branch("trigger_time_date", &trigger_time_date, "trigger_time_date/I");
    int year = ((trigger_time_date >> 25) & 0x7f) + 1980;
    int month = (trigger_time_date >> 21) & 0x0f;
    int day = (trigger_time_date >> 16) & 0x1f;
    int hour = (trigger_time_date >> 11) & 0x1f;
    int minute = (trigger_time_date >> 5) & 0x3f;
    int second = (trigger_time_date << 1) & 0x3e;
    header_info << "Date (DD.MM.YYYY): " << setw(2) << setfill('0') << day<< "."
		<< setw(2) << setfill('0') << month << "." 
		<< year << endl;
    header_info << "Time: "
		<< setw(2) << setfill('0') << hour << ":" 
		<< setw(2) << setfill('0') << minute << ":" 
		<< setw(2) << setfill('0') << second << endl;
    
    currentPos += fread(&trigger_coupling, 1, 2, pFile); 
    header_info << "For the external trigger input: ";
    switch(trigger_coupling)
      {
      case 1:
	header_info << "DC" << endl;
	break;
      case 2:
	header_info << "AC" << endl;
	break;
      }
    tree->Branch("trigger_coupling", &trigger_coupling, "trigger_coupling/S");
    
    currentPos += fread(&trigger_gain, 1, 2, pFile);
    header_info << "Input range (trigger gain): "<< "+/-" << input_range[trigger_gain] << "mv" << endl;
    double trigger_gain_v = 2*input_range[trigger_gain];
    tree->Branch("trigger_gain_v", &trigger_gain_v, "trigger_gain_v/D");
    
    currentPos += fread(&prob, 1, 2, pFile); 
    header_info << "Probe multiplier: " << probe_table[prob] << endl;;
    
    tree->Branch("prob", &prob, "prob/S");
    
    currentPos += fread(&inverted_data, 1, 2, pFile);
    switch(inverted_data)
      {
      case 0:
	header_info << "Normal data";
	break;
      case 1:
	header_info << "Inverted data";
	break;
      case 2:
	header_info << "Inverted and flipped data";
	break;
      }
    header_info << endl;
    tree->Branch("inverted_data", &inverted_data, "inverted_data/S");
    
    currentPos += fread(&board_type, 1, 2, pFile); 
    header_info << "Board type: ";
    for(int i = 0; i < 22; i++)
      {
	if(board_type_int[i] == board_type)
	  {
	    header_info << board_type_char[i];
	  }
      }
    header_info << endl;
    tree->Branch("board_type", &board_type, "board_type/S");
    
    currentPos += fread(&resolution_12_bits, 1, 2, pFile); 
    switch(resolution_12_bits)
      {
      case 0:
	header_info << "8 bit file format";
	break;
      case 1:
	header_info << "12/16 bit file format";
	break;
      }
    header_info << endl;
    tree->Branch("resolution_12_bits", &resolution_12_bits, "resolution_12_bits/S");
    
    currentPos += fread(&multiple_record, 1, 2, pFile); 
    header_info << "Saved data was captured in ";
    switch(multiple_record)
      {
      case 0:
	header_info << "normal mode";
	break;
      case 1:
	header_info << "Hardware multiple record mode";
	break;
      case 2:
	header_info << "Software multiple record mode";
	break;
      }
    header_info << endl;
    tree->Branch("multiple_record", &multiple_record, "multiple_record/S");
    
    currentPos += fread(&trigger_probe, 1, 2, pFile); 
    header_info << "Trigger probe: " << probe_table[prob] << endl;
    tree->Branch("trigger_probe", &trigger_probe, "trigger_probe/S");
    
    currentPos += fread(&unknwn, 1, 4, pFile);
    
    currentPos += fread(&sample_bits, 1, 2, pFile); 
    switch(sample_bits)
      {
      case 8:
	header_info << "8 bit";
	break;
      case 12:
	header_info << "12 bit";
	break;
      }
    header_info << " in the sampled data" << endl;
    tree->Branch("sample_bits", &sample_bits, "sample_bits/S");
    
    currentPos += fread(&extended_trigger_time, 1, 4, pFile); 
    header_info << "Trigger event occured at " << extended_trigger_time << " // format - hhhhhmmmmmmssssssddddddd" << endl;
    tree->Branch("extended_trigger_time", &extended_trigger_time, "extended_trigger_time/I");
    
    currentPos += fread(&imped_a, 1, 2, pFile); 
    header_info << "Impedance for Channel A: ";
    switch(imped_a)
      {
      case 0:
	header_info << "1 MegaOhm";
	break;
      case 0x10:
	header_info << "50 Ohm";
	break;
      }
    header_info << endl;
    tree->Branch("imped_a", &imped_a, "imped_a/S");
    
    currentPos += fread(&imped_b, 1, 2, pFile);
    header_info << "Impedance for Channel B: ";
    switch(imped_b)
      {
      case 0:
	header_info << "1 MegaOhm";
	break;
      case 0x10:
	header_info << "50 Ohm";
	break;
      }
    header_info << endl;
    tree->Branch("imped_b", &imped_b, "imped_b/S");
    
    currentPos += fread(&external_tbs, 1, 4, pFile); 
    header_info << "Time between samples: " << external_tbs << " nsec" << endl;
    tree->Branch("external_tbs", &external_tbs, "external_tbs/F");
    
    currentPos += fread(&external_clock_rate, 1, 4, pFile); 
    header_info << "Minimum sample rate when using external clock: " 
		<< external_clock_rate/1.e+09 << " sec" << endl;
    tree->Branch("external_clock_rate", &external_clock_rate, "external_clock_rate/F");
    
    currentPos += fread(&unknwn, 1, 18, pFile);
    
    currentPos += fread(&sample_offset, 1, 2, pFile); //cout << "sample_offset = " << sample_offset << endl;
    tree->Branch("sample_offset", &sample_offset, "sample_offset/S");
    
    currentPos += fread(&unknwn, 1, 2, pFile);
    
    currentPos += fread(&sample_resolution, 1, 2, pFile); //cout << "sample_resolution = " << sample_resolution << endl;
    tree->Branch("sample_resolution", &sample_resolution, "sample_resolution/S");
    
    currentPos += fread(&unknwn, 1, 6, pFile);
    
    currentPos += fread(&dc_offset, 1, 2, pFile);
    dc_offset = dc_offset*range/(1000.*2.);
    header_info << "DC offset: " << dc_offset << " mV." << endl;
    tree->Branch("dc_offset", &dc_offset, "dc_offset/S");
    
    currentPos += fread(&unknwn, 1, 113, pFile);
    
    cout << "Done!" << endl;
    
    //end of header
    cout << endl << "Reading of data (progress in %):";
    double percentage;
    
    TString leaf_name1 = "amplitude_f[";
    leaf_name1 += samp_per_msec;
    leaf_name1 += "]/F";
    tree->Branch("amplitude_f", amplitude_f, leaf_name1.Data());
    Int_t size_of_data_r = size_of_data;
    tree->Branch("size_of_data", &size_of_data_r,"size_of_data_r/I");
    Int_t samp_per_msec_r = samp_per_msec;
    tree->Branch("samp_per_msec", &samp_per_msec_r,"samp_per_msec_r/I");
    Double_t divide = size_of_data/samp_per_msec;
    tree->Branch("divide", &divide, "divide/D");
    
    
    int limit;
    
    if(size_of_data%samp_per_msec == 0)
      limit = size_of_data/samp_per_msec;
    else
      limit = size_of_data/samp_per_msec + 1;
    
    for(int j = 0; j < limit; j++)
      {
	for (int i = 0; i < samp_per_msec; i++)
	  {
	    if(i+j*samp_per_msec < size_of_data)
	      {
		currentPos += fread(&num_data_h, 1, 1, pFile);
		amplitude[i] = num_data_h;
		
		amplitude_f[i] = range/2. + dc_offset - (255 - amplitude[i])*range/256.;
		
	      }
	    else
	      {
		amplitude[i] = -999;
		amplitude_f[i] = -9999.;
	      }
	    
	    if ((i+j*samp_per_msec) % (int)(size_of_data/10.) == 0 && size_of_data >= samp_per_msec)
	      {
		percentage = 100 * ((double)(i+j*samp_per_msec) / (double)size_of_data);
		//cout << amplitude[0] << endl;
		cout << endl;
		
		if((i+j*samp_per_msec)==0){cout << " ";}
		cout << (int)(percentage + 0.5) << flush;
	      }
	    
	    if ((i+j*samp_per_msec) % (int)(size_of_data / 500.) == 0 && (i+j*samp_per_msec) % (int)(size_of_data/10.) != 0 && size_of_data >= samp_per_msec)
	      {
		cout << "." << flush;
	      }
	  }
	
	tree->Fill();
      }
    if(size_of_data < samp_per_msec){cout << "100%" << endl;}
    cout << endl << "Reading of file " << nameOfDataFile << " of size " 
	 << size_of_data << " is done. Tree has been filled. ";
    hfile = tree->GetCurrentFile();
    cout << "Write tree to file... ";
    
    hfile->Write();
    hfile->Close();
    cout << "Done." << endl;
    fclose(pFile);
    header_info.close();
  }
  
  inputdataList.close();
  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "----------------------All done!-------------------------" << endl;
  cout << "--------------------------------------------------------" << endl << endl;
}

void fill_sample_rate_table(long int srtab[n_sampleRateVal]){
   srtab[0] = 1; 
   srtab[1] = 2;
   srtab[2] = 5;
   srtab[3] = 10;
   srtab[4] = 20;
   srtab[5] = 50;
   srtab[6] = 100; 
   srtab[7] = 200;
   srtab[8] = 500;
   srtab[9] = 1000;
  srtab[10] = 2000;
  srtab[11] = 5000;
  srtab[12] = 10000; 
  srtab[13] = 20000;
  srtab[14] = 50000;
  srtab[15] = 100000;
  srtab[16] = 200000;
  srtab[17] = 500000;
  srtab[18] = 1e+06;
  srtab[19] = 2e+06;
  srtab[20] = 2.5e+06; 
  srtab[21] = 5e+06;
  srtab[22] = 1e+07;
  srtab[23] = 1.25e+07; 
  srtab[24] = 2e+07;
  srtab[25] = 2.5e+07; 
  srtab[26] = 3e+07;
  srtab[27] = 4e+07;
  srtab[28] = 5e+07;
  srtab[29] = 6e+07;
  srtab[30] = 6.5e+07; 
  srtab[31] = 8e+07;
  srtab[32] = 1e+08;
  srtab[33] = 1.2e+08; 
  srtab[34] = 1.25e+08;
  srtab[35] = 1.3e+08;
  srtab[36] = 1.5e+08;
  srtab[37] = 2e+08;
  srtab[38] = 2.5e+08; 
  srtab[39] = 3e+08;
  srtab[40] = 5e+08;
  srtab[41] = 1e+09;
  srtab[42] = 2e+09;
  srtab[43] = 4e+09;
  srtab[44] = 5e+09; 
  srtab[45] = 8e+09;
  srtab[46] = 1e+10; 
  srtab[47] = 0;
}

void fill_board_type_int(int btint[n_board_type]){
  //for(int i = 0;i<n_board_type_int;i++){
  //btint[i] = i;
  //}
  btint[0] = 0x0000;
  btint[1] = 0x0001;
  btint[2] = 0x0002;
  btint[3] = 0x0004;
  btint[4] = 0x0008;
  btint[5] = 0x000E;
  btint[6] = 0x0010;
  btint[7] = 0x001C;
  btint[8] = 0x0020;
  btint[9] = 0x0040;
  btint[10] = 0x0080;
  btint[11] = 0x0100;
  btint[12] = 0x0200;
  btint[13] = 0x0400;
  btint[14] = 0x0800;
  btint[15] = 0x1000;
  btint[16] = 0x1010;
  btint[17] = 0x2000;
  btint[18] = 0x2010;
  btint[19] = 0x4000;
  btint[20] = 0x8000;
  btint[21] = 0x8006;
}

void fill_board_type_char(TString btchar[n_board_type]){
   btchar[0] = "Unknown(pre \"file_version\" GS v2.25)";
   btchar[1] = "Compuscope 265";
   btchar[2] = "Compuscope 8500";
   btchar[3] = "Compuscope 8012-TYPE (CompuScope 8012-TYPE refers to CSx012 cards, version 2.0+.)";
   btchar[4] = "Compuscope 8012";
   btchar[5] = "Compuscope 12100, Compuscope 12130, Compuscope 1250";
   btchar[6] = "Compuscope PCI";
   btchar[7] = "Compuscope 8012/PCI";
   btchar[8] = "Compuscope 512";
   btchar[9] = "Compuscope 225";
  btchar[10] = "Compuscope 250 v1.8+";
  btchar[11] = "Compuscope LITE (pre-hardware version 1.5)";
  btchar[12] = "Compuscope 220";
  btchar[13] = "Compuscope 250";
  btchar[14] = "Compuscope LITE (hardware version 1.5 & up)";
  btchar[15] = "Compuscope 1012";
  btchar[16] = "Compuscope 1012/PCI";
  btchar[17] = "Compuscope 6012";
  btchar[18] = "Compuscope 6012/PCI";
  btchar[19] = "Compuscope 2125";
  btchar[20] = "Compuscope 1016";
  btchar[21] = "Compuscope 1610";
}
