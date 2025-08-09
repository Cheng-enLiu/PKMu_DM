#include <math.h>
#include <fstream>
using namespace std;

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"

#include <iostream>
#include <cmath>
#include <vector>

Int_t calc(

           //TString dm1ROOTfilename="smeared_rec_poca_DM_1.0e-02GeV_XS_1.0e-05cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_5.0e-02GeV_XS_5.0e-05cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_1.0e-01GeV_XS_1.0e-04cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_5.0e-01GeV_XS_5.0e-04cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_1.0e_00GeV_XS_1.0e-03cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_5.0e_00GeV_XS_5.0e-03cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_1.0e_01GeV_XS_1.0e-02cm2.root"
           //TString dm1ROOTfilename="smeared_rec_poca_DM_5.0e_01GeV_XS_5.0e-02cm2.root"
           TString dm1ROOTfilename="smeared_rec_poca_DM_1.0e_02GeV_XS_1.0e-01cm2.root"
           )
{
    //---读取用于重建图像的数据------------------------------
    TFile *f_dm1 = TFile::Open (dm1ROOTfilename,"read");
    TTree *Trec_dm1 = (TTree *)f_dm1->Get("Trec");
    Int_t nentries_dm1 = Trec_dm1->GetEntries();
    Double_t dm1x1,dm1y1,dm1z1,dm1x2,dm1y2,dm1z2,dm1x3,dm1y3,dm1z3,dm1x4,dm1y4,dm1z4;
    Double_t dm1x,dm1y,dm1z,dm1ang,dm1d;
    Int_t dm1Mark,dm1Pid;
    Double_t dm1xscatter, dm1yscatter, dm1zscatter;
    Trec_dm1->SetBranchAddress("x",&dm1x); 	  
    Trec_dm1->SetBranchAddress("y",&dm1y); 
    Trec_dm1->SetBranchAddress("z",&dm1z);
    //Trec_dm1->SetBranchAddress("xscatter",&dm1xscatter); 	  
    //Trec_dm1->SetBranchAddress("yscatter",&dm1yscatter); 
    //Trec_dm1->SetBranchAddress("zscatter",&dm1zscatter);
    Trec_dm1->SetBranchAddress("x1",&dm1x1);
    Trec_dm1->SetBranchAddress("y1",&dm1y1);
    Trec_dm1->SetBranchAddress("z1",&dm1z1);
    Trec_dm1->SetBranchAddress("x2",&dm1x2);
    Trec_dm1->SetBranchAddress("y2",&dm1y2);
    Trec_dm1->SetBranchAddress("z2",&dm1z2);
    Trec_dm1->SetBranchAddress("x3",&dm1x3);
    Trec_dm1->SetBranchAddress("y3",&dm1y3); 
    Trec_dm1->SetBranchAddress("z3",&dm1z3);   
    Trec_dm1->SetBranchAddress("x4",&dm1x4);
    Trec_dm1->SetBranchAddress("y4",&dm1y4);
    Trec_dm1->SetBranchAddress("z4",&dm1z4);
    Trec_dm1->SetBranchAddress("ang",&dm1ang);
    Trec_dm1->SetBranchAddress("d",&dm1d);
    Trec_dm1->SetBranchAddress("Pid",&dm1Pid);
    Trec_dm1->SetBranchAddress("Mark",&dm1Mark);
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------
    //--------------------------------------------------




    int dm1_count=0;
    int dm1_count_all=0;
    int i;

    for(i=0;i<nentries_dm1;i++){
        Trec_dm1->GetEntry(i);
        if(abs(dm1Pid)!=13) continue;
        if(dm1Mark!=1) continue;
        if(dm1x1>145||dm1x1<-145||dm1y1>145||dm1y1<-145||dm1x2>145||dm1x2<-145||dm1y2>145||dm1y2<-145) continue;
        if(dm1z>200||dm1z<-200) continue;
        if(dm1x>145||dm1x<-145) continue;
        if(dm1y>145||dm1y<-145) continue;
        if(dm1ang>0.5||dm1ang<0.05) continue;
        dm1_count_all++;
        if(dm1x3>145||dm1x3<-145||dm1x4>145||dm1x4<-145||dm1y3>145||dm1y3<-145||dm1y4>145||dm1y4<-145) continue;
        dm1_count++;
    }

    //int dm1_counts = dm1_count * (10000000./1000000);
    //int dm1_counts_all = dm1_count_all * (10000000./1000000);

    int dm1_counts = dm1_count;
    int dm1_counts_all = dm1_count_all;

    cout << "输入数：" << dm1_counts_all << endl;
    cout << "有效数：" << dm1_counts << endl;
    cout << "效率：" << std::fixed << std::setprecision(2) << ((dm1_counts * 1.0)/(dm1_counts_all * 1.0)) * 100 << endl;
    cout << "误差：" << std::fixed << std::setprecision(2) << (sqrt(((dm1_counts * 1.0) * ((dm1_counts * 1.0) + (dm1_counts_all * 1.0))) / (pow((dm1_counts_all * 1.0), 3)))) * 100 << endl;

    return dm1_count;
}