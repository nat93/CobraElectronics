#include "./src/CobraClass.hh"

using namespace std;

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        cout<<endl;
        cout<<"--> ERROR: Wrong number of input parameters!"<<endl;
        cout<<"--> [0] ./scriptname"<<endl;
        cout<<"--> [1] File with input file list"<<endl;
        cout<<"--> [2] Folder with input files"<<endl;
        cout<<"--> [3] File name for output data"<<endl;
        cout<<endl;
        assert(0);
    }

    TString inputDataFileList = argv[1];

    ifstream inputdataList;
    inputdataList.open(inputDataFileList.Data());
    assert(inputdataList.is_open());

    TString inputFileName, nameOfDataFile, outputFileName;

    outputFileName  =   argv[2];
    outputFileName  +=  argv[3];
    outputFileName  +=  ".root";

    TFile* outputfile = new TFile(outputFileName.Data(),"RECREATE");
    if (outputfile->IsZombie())
    {
        cout<<"--> ERROR: The output file is zombie!"<<endl;
        assert(0);
    }
    TTree* tree = new TTree("Tree", "CpFM data");

    while(inputdataList >> nameOfDataFile)
    {
        inputFileName  = argv[2] + nameOfDataFile;
        inputFileName  += ".root";

        CobraClass* fCobra = new CobraClass(inputFileName.Data());
        fCobra->Loop(tree);
        fCobra->~CobraClass();
        delete fCobra;
    }

    tree->Write();

    cout<<endl<<"--> The tree is written to the file: "<<outputfile->GetName()<<endl;
    outputfile->Close();

    return 0;
}
