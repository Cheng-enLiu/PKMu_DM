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

Int_t hist_copy_0730(TString expROOTfilename="250212_63d_air.root",
           TString simROOTfilename="smeared_rec_CryMuAna_0_749.root",
           TString dm1ROOTfilename="rec_poca_DM_1.0e-02GeV_XS_1.0e-05cm2.root",
           TString dm2ROOTfilename="rec_poca_DM_1.0e-01GeV_XS_1.0e-04cm2.root",
           TString dm3ROOTfilename="rec_poca_DM_1.0e_00GeV_XS_1.0e-03cm2.root",
           TString dm4ROOTfilename="rec_poca_DM_1.0e_01GeV_XS_1.0e-02cm2.root",
           TString dm5ROOTfilename="rec_poca_DM_1.0e_02GeV_XS_1.0e-01cm2.root",
           TString dm6ROOTfilename="rec_poca_DM_5.0e-02GeV_XS_5.0e-05cm2.root",
           TString dm7ROOTfilename="rec_poca_DM_5.0e-01GeV_XS_5.0e-04cm2.root",
           TString dm8ROOTfilename="rec_poca_DM_5.0e_00GeV_XS_5.0e-03cm2.root",
           TString dm9ROOTfilename="rec_poca_DM_5.0e_01GeV_XS_5.0e-02cm2.root")
{
    //---读取用于重建图像的数据------------------------------
    TFile *f_exp = TFile::Open (expROOTfilename,"read");
    TTree *Trec_fix = (TTree *)f_exp->Get("Trec");
    Int_t nentries_exp = Trec_fix->GetEntries(); 
    Double_t x,y,z,ang,d;
    Float_t x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    Int_t times[32];
    Trec_fix->SetBranchAddress("x",&x); 	  
    Trec_fix->SetBranchAddress("y",&y); 
    Trec_fix->SetBranchAddress("z",&z); 
    Trec_fix->SetBranchAddress("x1",&x1);
    Trec_fix->SetBranchAddress("y1",&y1);
    Trec_fix->SetBranchAddress("z1",&z1);
    Trec_fix->SetBranchAddress("x2",&x2);
    Trec_fix->SetBranchAddress("y2",&y2);
    Trec_fix->SetBranchAddress("z2",&z2);
    Trec_fix->SetBranchAddress("x3",&x3);
    Trec_fix->SetBranchAddress("y3",&y3); 
    Trec_fix->SetBranchAddress("z3",&z3);   
    Trec_fix->SetBranchAddress("x4",&x4);
    Trec_fix->SetBranchAddress("y4",&y4);
    Trec_fix->SetBranchAddress("z4",&z4);
    Trec_fix->SetBranchAddress("ang",&ang);
    Trec_fix->SetBranchAddress("d",&d);
    Trec_fix->SetBranchAddress("time", &times);
    //--------------------------------------------------
    Double_t x1s,y1s,z1s,x2s,y2s,z2s,x3s,y3s,z3s,x4s,y4s,z4s;
    Double_t xs,ys,zs,angs,ds;
    Int_t Pid;
    TFile *f_sim = TFile::Open (simROOTfilename,"read");
    TTree *Trec = (TTree *)f_sim->Get("Trec");
    Int_t nentries_sim = Trec->GetEntries(); 
    Trec->SetBranchAddress("x",&xs); 	  
    Trec->SetBranchAddress("y",&ys); 
    Trec->SetBranchAddress("z",&zs); 
    Trec->SetBranchAddress("x1",&x1s);
    Trec->SetBranchAddress("y1",&y1s);
    Trec->SetBranchAddress("z1",&z1s);
    Trec->SetBranchAddress("x2",&x2s);
    Trec->SetBranchAddress("y2",&y2s);
    Trec->SetBranchAddress("z2",&z2s);
    Trec->SetBranchAddress("x3",&x3s);
    Trec->SetBranchAddress("y3",&y3s); 
    Trec->SetBranchAddress("z3",&z3s);   
    Trec->SetBranchAddress("x4",&x4s);
    Trec->SetBranchAddress("y4",&y4s);
    Trec->SetBranchAddress("z4",&z4s);
    Trec->SetBranchAddress("ang",&angs);
    Trec->SetBranchAddress("d",&ds);
    Trec->SetBranchAddress("Pid",&Pid);
    //--------------------------------------------------
    TFile *f_dm1 = TFile::Open (dm1ROOTfilename,"read");
    TTree *Trec_dm1 = (TTree *)f_dm1->Get("Trec");
    Int_t nentries_dm1 = Trec_dm1->GetEntries();
    Double_t dm1x1,dm1y1,dm1z1,dm1x2,dm1y2,dm1z2,dm1x3,dm1y3,dm1z3,dm1x4,dm1y4,dm1z4;
    Double_t dm1x,dm1y,dm1z,dm1ang,dm1d;
    Int_t dm1Mark,dm1Pid;
    Trec_dm1->SetBranchAddress("x",&dm1x); 	  
    Trec_dm1->SetBranchAddress("y",&dm1y); 
    Trec_dm1->SetBranchAddress("z",&dm1z); 
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
    TFile *f_dm2 = TFile::Open (dm2ROOTfilename,"read");
    TTree *Trec_dm2 = (TTree *)f_dm2->Get("Trec");
    Int_t nentries_dm2 = Trec_dm2->GetEntries();
    Double_t dm2x1,dm2y1,dm2z1,dm2x2,dm2y2,dm2z2,dm2x3,dm2y3,dm2z3,dm2x4,dm2y4,dm2z4;
    Double_t dm2x,dm2y,dm2z,dm2ang,dm2d;
    Int_t dm2Mark,dm2Pid;
    Trec_dm2->SetBranchAddress("x",&dm2x); 	  
    Trec_dm2->SetBranchAddress("y",&dm2y); 
    Trec_dm2->SetBranchAddress("z",&dm2z); 
    Trec_dm2->SetBranchAddress("x1",&dm2x1);
    Trec_dm2->SetBranchAddress("y1",&dm2y1);
    Trec_dm2->SetBranchAddress("z1",&dm2z1);
    Trec_dm2->SetBranchAddress("x2",&dm2x2);
    Trec_dm2->SetBranchAddress("y2",&dm2y2);
    Trec_dm2->SetBranchAddress("z2",&dm2z2);
    Trec_dm2->SetBranchAddress("x3",&dm2x3);
    Trec_dm2->SetBranchAddress("y3",&dm2y3); 
    Trec_dm2->SetBranchAddress("z3",&dm2z3);   
    Trec_dm2->SetBranchAddress("x4",&dm2x4);
    Trec_dm2->SetBranchAddress("y4",&dm2y4);
    Trec_dm2->SetBranchAddress("z4",&dm2z4);
    Trec_dm2->SetBranchAddress("ang",&dm2ang);
    Trec_dm2->SetBranchAddress("d",&dm2d);
    Trec_dm2->SetBranchAddress("Pid",&dm2Pid);
    Trec_dm2->SetBranchAddress("Mark",&dm2Mark);
    //--------------------------------------------------
    TFile *f_dm3 = TFile::Open (dm3ROOTfilename,"read");
    TTree *Trec_dm3 = (TTree *)f_dm3->Get("Trec");
    Int_t nentries_dm3 = Trec_dm3->GetEntries();
    Double_t dm3x1,dm3y1,dm3z1,dm3x2,dm3y2,dm3z2,dm3x3,dm3y3,dm3z3,dm3x4,dm3y4,dm3z4;
    Double_t dm3x,dm3y,dm3z,dm3ang,dm3d;
    Int_t dm3Mark,dm3Pid;
    Trec_dm3->SetBranchAddress("x",&dm3x); 	  
    Trec_dm3->SetBranchAddress("y",&dm3y); 
    Trec_dm3->SetBranchAddress("z",&dm3z); 
    Trec_dm3->SetBranchAddress("x1",&dm3x1);
    Trec_dm3->SetBranchAddress("y1",&dm3y1);
    Trec_dm3->SetBranchAddress("z1",&dm3z1);
    Trec_dm3->SetBranchAddress("x2",&dm3x2);
    Trec_dm3->SetBranchAddress("y2",&dm3y2);
    Trec_dm3->SetBranchAddress("z2",&dm3z2);
    Trec_dm3->SetBranchAddress("x3",&dm3x3);
    Trec_dm3->SetBranchAddress("y3",&dm3y3); 
    Trec_dm3->SetBranchAddress("z3",&dm3z3);   
    Trec_dm3->SetBranchAddress("x4",&dm3x4);
    Trec_dm3->SetBranchAddress("y4",&dm3y4);
    Trec_dm3->SetBranchAddress("z4",&dm3z4);
    Trec_dm3->SetBranchAddress("ang",&dm3ang);
    Trec_dm3->SetBranchAddress("d",&dm3d);
    Trec_dm3->SetBranchAddress("Pid",&dm3Pid);
    Trec_dm3->SetBranchAddress("Mark",&dm3Mark);
    //--------------------------------------------------
    TFile *f_dm4 = TFile::Open (dm4ROOTfilename,"read");
    TTree *Trec_dm4 = (TTree *)f_dm4->Get("Trec");
    Int_t nentries_dm4 = Trec_dm4->GetEntries();
    Double_t dm4x1,dm4y1,dm4z1,dm4x2,dm4y2,dm4z2,dm4x3,dm4y3,dm4z3,dm4x4,dm4y4,dm4z4;
    Double_t dm4x,dm4y,dm4z,dm4ang,dm4d;
    Int_t dm4Mark,dm4Pid;
    Trec_dm4->SetBranchAddress("x",&dm4x); 	  
    Trec_dm4->SetBranchAddress("y",&dm4y); 
    Trec_dm4->SetBranchAddress("z",&dm4z); 
    Trec_dm4->SetBranchAddress("x1",&dm4x1);
    Trec_dm4->SetBranchAddress("y1",&dm4y1);
    Trec_dm4->SetBranchAddress("z1",&dm4z1);
    Trec_dm4->SetBranchAddress("x2",&dm4x2);
    Trec_dm4->SetBranchAddress("y2",&dm4y2);
    Trec_dm4->SetBranchAddress("z2",&dm4z2);
    Trec_dm4->SetBranchAddress("x3",&dm4x3);
    Trec_dm4->SetBranchAddress("y3",&dm4y3); 
    Trec_dm4->SetBranchAddress("z3",&dm4z3);   
    Trec_dm4->SetBranchAddress("x4",&dm4x4);
    Trec_dm4->SetBranchAddress("y4",&dm4y4);
    Trec_dm4->SetBranchAddress("z4",&dm4z4);
    Trec_dm4->SetBranchAddress("ang",&dm4ang);
    Trec_dm4->SetBranchAddress("d",&dm4d);
    Trec_dm4->SetBranchAddress("Pid",&dm4Pid);
    Trec_dm4->SetBranchAddress("Mark",&dm4Mark);
    //--------------------------------------------------
    TFile *f_dm5 = TFile::Open (dm5ROOTfilename,"read");
    TTree *Trec_dm5 = (TTree *)f_dm5->Get("Trec");
    Int_t nentries_dm5 = Trec_dm5->GetEntries();
    Double_t dm5x1,dm5y1,dm5z1,dm5x2,dm5y2,dm5z2,dm5x3,dm5y3,dm5z3,dm5x4,dm5y4,dm5z4;
    Double_t dm5x,dm5y,dm5z,dm5ang,dm5d;
    Int_t dm5Mark,dm5Pid;
    Trec_dm5->SetBranchAddress("x",&dm5x); 	  
    Trec_dm5->SetBranchAddress("y",&dm5y); 
    Trec_dm5->SetBranchAddress("z",&dm5z); 
    Trec_dm5->SetBranchAddress("x1",&dm5x1);
    Trec_dm5->SetBranchAddress("y1",&dm5y1);
    Trec_dm5->SetBranchAddress("z1",&dm5z1);
    Trec_dm5->SetBranchAddress("x2",&dm5x2);
    Trec_dm5->SetBranchAddress("y2",&dm5y2);
    Trec_dm5->SetBranchAddress("z2",&dm5z2);
    Trec_dm5->SetBranchAddress("x3",&dm5x3);
    Trec_dm5->SetBranchAddress("y3",&dm5y3); 
    Trec_dm5->SetBranchAddress("z3",&dm5z3);   
    Trec_dm5->SetBranchAddress("x4",&dm5x4);
    Trec_dm5->SetBranchAddress("y4",&dm5y4);
    Trec_dm5->SetBranchAddress("z4",&dm5z4);
    Trec_dm5->SetBranchAddress("ang",&dm5ang);
    Trec_dm5->SetBranchAddress("d",&dm5d);
    Trec_dm5->SetBranchAddress("Pid",&dm5Pid);
    Trec_dm5->SetBranchAddress("Mark",&dm5Mark);
    //--------------------------------------------------
    TFile *f_dm6 = TFile::Open (dm6ROOTfilename,"read");
    TTree *Trec_dm6 = (TTree *)f_dm6->Get("Trec");
    Int_t nentries_dm6 = Trec_dm6->GetEntries();
    Double_t dm6x1,dm6y1,dm6z1,dm6x2,dm6y2,dm6z2,dm6x3,dm6y3,dm6z3,dm6x4,dm6y4,dm6z4;
    Double_t dm6x,dm6y,dm6z,dm6ang,dm6d;
    Int_t dm6Mark,dm6Pid;
    Trec_dm6->SetBranchAddress("x",&dm6x); 	  
    Trec_dm6->SetBranchAddress("y",&dm6y); 
    Trec_dm6->SetBranchAddress("z",&dm6z); 
    Trec_dm6->SetBranchAddress("x1",&dm6x1);
    Trec_dm6->SetBranchAddress("y1",&dm6y1);
    Trec_dm6->SetBranchAddress("z1",&dm6z1);
    Trec_dm6->SetBranchAddress("x2",&dm6x2);
    Trec_dm6->SetBranchAddress("y2",&dm6y2);
    Trec_dm6->SetBranchAddress("z2",&dm6z2);
    Trec_dm6->SetBranchAddress("x3",&dm6x3);
    Trec_dm6->SetBranchAddress("y3",&dm6y3); 
    Trec_dm6->SetBranchAddress("z3",&dm6z3);   
    Trec_dm6->SetBranchAddress("x4",&dm6x4);
    Trec_dm6->SetBranchAddress("y4",&dm6y4);
    Trec_dm6->SetBranchAddress("z4",&dm6z4);
    Trec_dm6->SetBranchAddress("ang",&dm6ang);
    Trec_dm6->SetBranchAddress("d",&dm6d);
    Trec_dm6->SetBranchAddress("Pid",&dm6Pid);
    Trec_dm6->SetBranchAddress("Mark",&dm6Mark);
    //--------------------------------------------------
    TFile *f_dm7 = TFile::Open (dm7ROOTfilename,"read");
    TTree *Trec_dm7 = (TTree *)f_dm7->Get("Trec");
    Int_t nentries_dm7 = Trec_dm7->GetEntries();
    Double_t dm7x1,dm7y1,dm7z1,dm7x2,dm7y2,dm7z2,dm7x3,dm7y3,dm7z3,dm7x4,dm7y4,dm7z4;
    Double_t dm7x,dm7y,dm7z,dm7ang,dm7d;
    Int_t dm7Mark,dm7Pid;
    Trec_dm7->SetBranchAddress("x",&dm7x); 	  
    Trec_dm7->SetBranchAddress("y",&dm7y); 
    Trec_dm7->SetBranchAddress("z",&dm7z); 
    Trec_dm7->SetBranchAddress("x1",&dm7x1);
    Trec_dm7->SetBranchAddress("y1",&dm7y1);
    Trec_dm7->SetBranchAddress("z1",&dm7z1);
    Trec_dm7->SetBranchAddress("x2",&dm7x2);
    Trec_dm7->SetBranchAddress("y2",&dm7y2);
    Trec_dm7->SetBranchAddress("z2",&dm7z2);
    Trec_dm7->SetBranchAddress("x3",&dm7x3);
    Trec_dm7->SetBranchAddress("y3",&dm7y3); 
    Trec_dm7->SetBranchAddress("z3",&dm7z3);   
    Trec_dm7->SetBranchAddress("x4",&dm7x4);
    Trec_dm7->SetBranchAddress("y4",&dm7y4);
    Trec_dm7->SetBranchAddress("z4",&dm7z4);
    Trec_dm7->SetBranchAddress("ang",&dm7ang);
    Trec_dm7->SetBranchAddress("d",&dm7d);
    Trec_dm7->SetBranchAddress("Pid",&dm7Pid);
    Trec_dm7->SetBranchAddress("Mark",&dm7Mark);
    //--------------------------------------------------
    TFile *f_dm8 = TFile::Open (dm8ROOTfilename,"read");
    TTree *Trec_dm8 = (TTree *)f_dm8->Get("Trec");
    Int_t nentries_dm8 = Trec_dm8->GetEntries();
    Double_t dm8x1,dm8y1,dm8z1,dm8x2,dm8y2,dm8z2,dm8x3,dm8y3,dm8z3,dm8x4,dm8y4,dm8z4;
    Double_t dm8x,dm8y,dm8z,dm8ang,dm8d;
    Int_t dm8Mark,dm8Pid;
    Trec_dm8->SetBranchAddress("x",&dm8x); 	  
    Trec_dm8->SetBranchAddress("y",&dm8y); 
    Trec_dm8->SetBranchAddress("z",&dm8z); 
    Trec_dm8->SetBranchAddress("x1",&dm8x1);
    Trec_dm8->SetBranchAddress("y1",&dm8y1);
    Trec_dm8->SetBranchAddress("z1",&dm8z1);
    Trec_dm8->SetBranchAddress("x2",&dm8x2);
    Trec_dm8->SetBranchAddress("y2",&dm8y2);
    Trec_dm8->SetBranchAddress("z2",&dm8z2);
    Trec_dm8->SetBranchAddress("x3",&dm8x3);
    Trec_dm8->SetBranchAddress("y3",&dm8y3); 
    Trec_dm8->SetBranchAddress("z3",&dm8z3);   
    Trec_dm8->SetBranchAddress("x4",&dm8x4);
    Trec_dm8->SetBranchAddress("y4",&dm8y4);
    Trec_dm8->SetBranchAddress("z4",&dm8z4);
    Trec_dm8->SetBranchAddress("ang",&dm8ang);
    Trec_dm8->SetBranchAddress("d",&dm8d);
    Trec_dm8->SetBranchAddress("Pid",&dm8Pid);
    Trec_dm8->SetBranchAddress("Mark",&dm8Mark);
    //--------------------------------------------------
    TFile *f_dm9 = TFile::Open (dm9ROOTfilename,"read");
    TTree *Trec_dm9 = (TTree *)f_dm9->Get("Trec");
    Int_t nentries_dm9 = Trec_dm9->GetEntries();
    Double_t dm9x1,dm9y1,dm9z1,dm9x2,dm9y2,dm9z2,dm9x3,dm9y3,dm9z3,dm9x4,dm9y4,dm9z4;
    Double_t dm9x,dm9y,dm9z,dm9ang,dm9d;
    Int_t dm9Mark,dm9Pid;
    Trec_dm9->SetBranchAddress("x",&dm9x); 	  
    Trec_dm9->SetBranchAddress("y",&dm9y); 
    Trec_dm9->SetBranchAddress("z",&dm9z); 
    Trec_dm9->SetBranchAddress("x1",&dm9x1);
    Trec_dm9->SetBranchAddress("y1",&dm9y1);
    Trec_dm9->SetBranchAddress("z1",&dm9z1);
    Trec_dm9->SetBranchAddress("x2",&dm9x2);
    Trec_dm9->SetBranchAddress("y2",&dm9y2);
    Trec_dm9->SetBranchAddress("z2",&dm9z2);
    Trec_dm9->SetBranchAddress("x3",&dm9x3);
    Trec_dm9->SetBranchAddress("y3",&dm9y3); 
    Trec_dm9->SetBranchAddress("z3",&dm9z3);   
    Trec_dm9->SetBranchAddress("x4",&dm9x4);
    Trec_dm9->SetBranchAddress("y4",&dm9y4);
    Trec_dm9->SetBranchAddress("z4",&dm9z4);
    Trec_dm9->SetBranchAddress("ang",&dm9ang);
    Trec_dm9->SetBranchAddress("d",&dm9d);
    Trec_dm9->SetBranchAddress("Pid",&dm9Pid);
    Trec_dm9->SetBranchAddress("Mark",&dm9Mark);
    //--------------------------------------------------

    TCanvas *canvas = new TCanvas("canvas", "ang", 800, 600);

    TH1F *h_const_0p01 = new TH1F("h_const_0p01", "h_const_0p01", 50, 0.05, 0.5);
    TH1F *h_const_0p1 = new TH1F("h_const_0p1", "h_const_0p1", 50, 0.05, 0.5);
    TH1F *h_const_1 = new TH1F("h_const_1", "h_const_1", 50, 0.05, 0.5);
    TH1F *h_const_10 = new TH1F("h_const_10", "h_const_10", 50, 0.05, 0.5);
    TH1F *h_const_100 = new TH1F("h_const_100", "h_const_100", 50, 0.05, 0.5);
    TH1F *h_const_0p05 = new TH1F("h_const_0p05", "h_const_0p05", 50, 0.05, 0.5);
    TH1F *h_const_0p5 = new TH1F("h_const_0p5", "h_const_0p5", 50, 0.05, 0.5);
    TH1F *h_const_5 = new TH1F("h_const_5", "h_const_5", 50, 0.05, 0.5);
    TH1F *h_const_50 = new TH1F("h_const_50", "h_const_50", 50, 0.05, 0.5);
    
    TH1F *h_air = new TH1F("h_air", "distribution of scattering angle #theta", 50, 0.05, 0.5);
    TH1F *h_mu = new TH1F("h_mu", "distribution of scattering angle #theta", 50, 0.05, 0.5);
    TH1F *h_e = new TH1F("h_e", "distribution of scattering angle #theta", 50, 0.05, 0.5);
    TH1F *h_other = new TH1F("h_other", "distribution of scattering angle #theta", 50, 0.05, 0.5);

    TH1F *h_obs = new TH1F("h_obs", "distribution of scattering angle #theta", 50, 0.05, 0.5);

    int exp_count=0;
    int sim_count=0;
    int dm1_count=0;
    int dm2_count=0;
    int dm3_count=0;
    int dm4_count=0;
    int dm5_count=0;
    int dm6_count=0;
    int dm7_count=0;
    int dm8_count=0;
    int dm9_count=0;
    int exp_counts=0;
    int sim_counts=0;
    int dm1_counts=0;
    int dm2_counts=0;
    int dm3_counts=0;
    int dm4_counts=0;
    int dm5_counts=0;
    int dm6_counts=0;
    int dm7_counts=0;
    int dm8_counts=0;
    int dm9_counts=0;
    int mu_count=0;
    int e_count=0;
    int other_count=0;
    int mu_counts=0;
    int e_counts=0;
    int other_counts=0;
    int i;
    int dmflux=0;

    for(i=0;i<nentries_exp;i++){
        Trec_fix->GetEntry(i);
	    if(x1>145||x1<-145||x2>145||x2<-145||x3>145||x3<-145||x4>145||x4<-145||y1>145||y1<-145||y2>145||y2<-145||y3>145||y3<-145||y4>145||y4<-145) continue;
        if (times[1]<0 || times[10]<0 || times[11]<0 || times[12]<0 || times[13]<0 || times[14]<0 || times[15]<0 || times[16]<0 || times[17]<0 || times[18]<0 || times[19]<0 ||
            times[20]<0 || times[21]<0 || times[22]<0 || times[23]<0 || times[24]<0 || times[25]<0 || times[28]<0 || times[29]<0 || times[30]<0 || times[31]<0) 
            continue;
        if(z>110||z<-110) continue;
        if(x>110||x<-110) continue;
        if(y>110||y<-110) continue;
        if(ang<0.05||ang>0.5)continue;
        exp_count++;
        exp_counts++;
        h_obs->Fill(ang);
    }

    for(i=0;i<nentries_sim;i++){
        Trec->GetEntry(i);
        if(x1s>145||x1s<-145||x2s>145||x2s<-145||x3s>145||x3s<-145||x4s>145||x4s<-145||y1s>145||y1s<-145||y2s>145||y2s<-145||y3s>145||y3s<-145||y4s>145||y4s<-145) continue;
        if(zs>110||zs<-100) continue;
        if(xs>110||xs<-110) continue;
        if(ys>110||ys<-110) continue;
        if(angs<0.05||angs>0.5) continue;
        sim_count++;
        if(abs(Pid)==11) e_count++;
        if(abs(Pid)==13) mu_count++;
        if(abs(Pid)!=11&&abs(Pid)!=13) other_count++;
        sim_counts++;
        if(abs(Pid)==13){
            mu_counts++;
            h_mu->Fill(angs);
        }
        if(abs(Pid)==11){
            e_counts++;
            h_e->Fill(angs);
        }
        if(abs(Pid)!=13&&abs(Pid)!=11){
            other_counts++;
            h_other->Fill(angs);
        }
    }
    for (int bin = 1; bin <= 50; bin++){
        if (h_mu->GetBinContent(bin) == 0 ) h_mu->SetBinContent(bin,(h_mu->GetBinContent(bin - 1) + h_mu->GetBinContent(bin + 1)) * 1.0 / 2);
        if (h_e->GetBinContent(bin) == 0 ) h_e->SetBinContent(bin,(h_e->GetBinContent(bin - 1) + h_e->GetBinContent(bin + 1)) * 1.0 / 2);
        if (h_other->GetBinContent(bin) == 0 ) h_other->SetBinContent(bin,(h_other->GetBinContent(bin - 1) + h_other->GetBinContent(bin + 1)) * 1.0 / 2);
    }

    for(i=0;i<nentries_dm1;i++){
        Trec_dm1->GetEntry(i);
        if(dm1x1>145||dm1x1<-145||dm1x2>145||dm1x2<-145||dm1x3>145||dm1x3<-145||dm1x4>145||dm1x4<-145||dm1y1>145||dm1y1<-145||dm1y2>145||dm1y2<-145||dm1y3>145||dm1y3<-145||dm1y4>145||dm1y4<-145) continue;
        if(dm1z>110||dm1z<-110) continue;
        if(dm1x>110||dm1x<-110) continue;
        if(dm1y>110||dm1y<-110) continue;
        if(dm1ang<0.05||dm1ang>0.5) continue;
        if(abs(dm1Pid)==13) dm1_count++;
        if(dm1Mark==1) {
            if(abs(dm1Pid)==13){
                dm1_counts++;
                h_const_0p01->Fill(dm1ang);
            }
        }
    }

    for(i=0;i<nentries_dm2;i++){
        Trec_dm2->GetEntry(i);
        if(dm2x1>145||dm2x1<-145||dm2x2>145||dm2x2<-145||dm2x3>145||dm2x3<-145||dm2x4>145||dm2x4<-145||dm2y1>145||dm2y1<-145||dm2y2>145||dm2y2<-145||dm2y3>145||dm2y3<-145||dm2y4>145||dm2y4<-145) continue;
        if(dm2z>110||dm2z<-110) continue;
        if(dm2x>110||dm2x<-110) continue;
        if(dm2y>110||dm2y<-110) continue;
        if(dm2ang<0.05||dm2ang>0.5) continue;
        if(abs(dm2Pid)==13) dm2_count++;
        if(dm2Mark==1) {
            if(abs(dm2Pid)==13){
                dm2_counts++;
                h_const_0p1->Fill(dm2ang);
            }
        }
    }

    for(i=0;i<nentries_dm3;i++){
        Trec_dm3->GetEntry(i);
        if(dm3x1>145||dm3x1<-145||dm3x2>145||dm3x2<-145||dm3x3>145||dm3x3<-145||dm3x4>145||dm3x4<-145||dm3y1>145||dm3y1<-145||dm3y2>145||dm3y2<-145||dm3y3>145||dm3y3<-145||dm3y4>145||dm3y4<-145) continue;
        if(dm3z>110||dm3z<-110) continue;
        if(dm3x>110||dm3x<-110) continue;
        if(dm3y>110||dm3y<-110) continue;
        if(dm3ang<0.05||dm3ang>0.5) continue;
        if(abs(dm3Pid)==13) dm3_count++;
        if(abs(dm3Pid)==13) dmflux++;
        if(dm3Mark==1) {
            if(abs(dm3Pid)==13){
                dm3_counts++;
                h_const_1->Fill(dm3ang);
            }
        }
    }

    for(i=0;i<nentries_dm4;i++){
        Trec_dm4->GetEntry(i);
	    if(dm4x1>145||dm4x1<-145||dm4x2>145||dm4x2<-145||dm4x3>145||dm4x3<-145||dm4x4>145||dm4x4<-145||dm4y1>145||dm4y1<-145||dm4y2>145||dm4y2<-145||dm4y3>145||dm4y3<-145||dm4y4>145||dm4y4<-145) continue;
        if(dm4z>110||dm4z<-110) continue;
        if(dm4x>110||dm4x<-110) continue;
        if(dm4y>110||dm4y<-110) continue;
        if(dm4ang<0.05||dm4ang>0.5) continue;
        if(abs(dm4Pid)==13) dm4_count++;
        if(dm4Mark==1) {
            if(abs(dm4Pid)==13){
                dm4_counts++;
                h_const_10->Fill(dm4ang);
            }
        }
    }

    for(i=0;i<nentries_dm5;i++){
        Trec_dm5->GetEntry(i);
	    if(dm5x1>145||dm5x1<-145||dm5x2>145||dm5x2<-145||dm5x3>145||dm5x3<-145||dm5x4>145||dm5x4<-145||dm5y1>145||dm5y1<-145||dm5y2>145||dm5y2<-145||dm5y3>145||dm5y3<-145||dm5y4>145||dm5y4<-145) continue;
        if(dm5z>110||dm5z<-110) continue;
        if(dm5x>110||dm5x<-110) continue;
        if(dm5y>110||dm5y<-110) continue;
        if(dm5ang<0.05||dm5ang>0.5) continue;
        if(abs(dm5Pid)==13) dm5_count++;
        if(dm5Mark==1) {
            if(abs(dm5Pid)==13){
                dm5_counts++;
                h_const_100->Fill(dm5ang);
            }
        }
    }

    for(i=0;i<nentries_dm6;i++){
        Trec_dm6->GetEntry(i);
	    if(dm6x1>145||dm6x1<-145||dm6x2>145||dm6x2<-145||dm6x3>145||dm6x3<-145||dm6x4>145||dm6x4<-145||dm6y1>145||dm6y1<-145||dm6y2>145||dm6y2<-145||dm6y3>145||dm6y3<-145||dm6y4>145||dm6y4<-145) continue;
        if(dm6z>110||dm6z<-110) continue;
        if(dm6x>110||dm6x<-110) continue;
        if(dm6y>110||dm6y<-110) continue;
        if(dm6ang<0.05||dm6ang>0.5) continue;
        if(abs(dm6Pid)==13) dm6_count++;
        if(dm6Mark==1) {
            if(abs(dm6Pid)==13){
                dm6_counts++;
                h_const_0p05->Fill(dm6ang);
            }
        }
    }

    for(i=0;i<nentries_dm7;i++){
        Trec_dm7->GetEntry(i);
	    if(dm7x1>145||dm7x1<-145||dm7x2>145||dm7x2<-145||dm7x3>145||dm7x3<-145||dm7x4>145||dm7x4<-145||dm7y1>145||dm7y1<-145||dm7y2>145||dm7y2<-145||dm7y3>145||dm7y3<-145||dm7y4>145||dm7y4<-145) continue;
        if(dm7z>110||dm7z<-110) continue;
        if(dm7x>110||dm7x<-110) continue;
        if(dm7y>110||dm7y<-110) continue;
        if(dm7ang<0.05||dm7ang>0.5) continue;
        if(abs(dm7Pid)==13) dm7_count++;
        if(dm7Mark==1) {
            if(abs(dm7Pid)==13){
                dm7_counts++;
                h_const_0p5->Fill(dm7ang);
            }
        }
    }

    for(i=0;i<nentries_dm8;i++){
        Trec_dm8->GetEntry(i);
	    if(dm8x1>145||dm8x1<-145||dm8x2>145||dm8x2<-145||dm8x3>145||dm8x3<-145||dm8x4>145||dm8x4<-145||dm8y1>145||dm8y1<-145||dm8y2>145||dm8y2<-145||dm8y3>145||dm8y3<-145||dm8y4>145||dm8y4<-145) continue;
        if(dm8z>110||dm8z<-110) continue;
        if(dm8x>110||dm8x<-110) continue;
        if(dm8y>110||dm8y<-110) continue;
        if(dm8ang<0.05||dm8ang>0.5) continue;
        if(abs(dm8Pid)==13) dm8_count++;
        if(dm8Mark==1) {
            if(abs(dm8Pid)==13){
                dm8_counts++;
                h_const_5->Fill(dm8ang);
            }
        }
    }

    for(i=0;i<nentries_dm9;i++){
        Trec_dm9->GetEntry(i);
	    if(dm9x1>145||dm9x1<-145||dm9x2>145||dm9x2<-145||dm9x3>145||dm9x3<-145||dm9x4>145||dm9x4<-145||dm9y1>145||dm9y1<-145||dm9y2>145||dm9y2<-145||dm9y3>145||dm9y3<-145||dm9y4>145||dm9y4<-145) continue;
        if(dm9z>110||dm9z<-110) continue;
        if(dm9x>110||dm9x<-110) continue;
        if(dm9y>110||dm9y<-110) continue;
        if(dm9ang<0.05||dm9ang>0.5) continue;
        if(abs(dm9Pid)==13) dm9_count++;
        if(dm9Mark==1) {
            if(abs(dm9Pid)==13){
                dm9_counts++;
                h_const_50->Fill(dm9ang);
            }
        }
    }

    cout<<"exp_count="<<exp_count<<endl;
    cout<<"sim_count="<<sim_count<<endl;
    cout<<"mu_count="<<mu_count<<endl;
    cout<<"e_count="<<e_count<<endl;
    cout<<"other_count"<<other_count<<endl;
    cout<<"dm1_count="<<dm1_count<<endl;
    cout<<"dm2_count="<<dm2_count<<endl;
    cout<<"dm3_count="<<dm3_count<<endl;
    cout<<"dm4_count="<<dm4_count<<endl;
    cout<<"dm5_count="<<dm5_count<<endl;
    cout<<"dm6_count="<<dm6_count<<endl;
    cout<<"dm7_count="<<dm7_count<<endl;
    cout<<"dm8_count="<<dm8_count<<endl;
    cout<<"dm9_count="<<dm9_count<<endl;
    cout<<"exp_counts="<<exp_counts<<endl;
    cout<<"sim_counts="<<sim_counts<<endl;
    cout<<"mu_counts="<<mu_counts<<endl;
    cout<<"e_counts="<<e_counts<<endl;
    cout<<"other_counts"<<other_counts<<endl;
    cout<<"dm1_counts="<<dm1_counts<<endl;
    cout<<"dm2_counts="<<dm2_counts<<endl;
    cout<<"dm3_counts="<<dm3_counts<<endl;
    cout<<"dm4_counts="<<dm4_counts<<endl;
    cout<<"dm5_counts="<<dm5_counts<<endl;
    cout<<"dm6_counts="<<dm6_counts<<endl;
    cout<<"dm7_counts="<<dm7_counts<<endl;
    cout<<"dm8_counts="<<dm8_counts<<endl;
    cout<<"dm9_counts="<<dm9_counts<<endl;
    cout<<"dmflux="<<dmflux<<endl;

    TH1F *h_mu_alphaUp = new TH1F(*h_mu);
    h_mu_alphaUp->SetName("h_mu_alphaUp");
    TH1F *h_mu_alphaDown = new TH1F(*h_mu);
    h_mu_alphaDown->SetName("h_mu_alphaDown");
    TH1F *h_e_alphaUp = new TH1F(*h_e);
    h_e_alphaUp->SetName("h_e_alphaUp");
    TH1F *h_e_alphaDown = new TH1F(*h_e);
    h_e_alphaDown->SetName("h_e_alphaDown");

    TH1F *h_mu_sigmaUp = new TH1F(*h_mu);
    h_mu_sigmaUp->SetName("h_mu_sigmaUp");
    TH1F *h_mu_sigmaDown = new TH1F(*h_mu);
    h_mu_sigmaDown->SetName("h_mu_sigmaDown");
    TH1F *h_e_sigmaUp = new TH1F(*h_e);
    h_e_sigmaUp->SetName("h_e_sigmaUp");
    TH1F *h_e_sigmaDown = new TH1F(*h_e);
    h_e_sigmaDown->SetName("h_e_sigmaDown");
    TH1F *h_other_sigmaUp = new TH1F(*h_other);
    h_other_sigmaUp->SetName("h_other_sigmaUp");
    TH1F *h_other_sigmaDown = new TH1F(*h_other);
    h_other_sigmaDown->SetName("h_other_sigmaDown");

    TH1F *h_const_0p01_alphaUp = new TH1F(*h_const_0p01);
    h_const_0p01_alphaUp->SetName("h_const_0p01_alphaUp");
    TH1F *h_const_0p01_alphaDown = new TH1F(*h_const_0p01);
    h_const_0p01_alphaDown->SetName("h_const_0p01_alphaDown");
    TH1F *h_const_0p1_alphaUp = new TH1F(*h_const_0p1);
    h_const_0p1_alphaUp->SetName("h_const_0p1_alphaUp");
    TH1F *h_const_0p1_alphaDown = new TH1F(*h_const_0p1);
    h_const_0p1_alphaDown->SetName("h_const_0p1_alphaDown");
    TH1F *h_const_1_alphaUp = new TH1F(*h_const_1);
    h_const_1_alphaUp->SetName("h_const_1_alphaUp");
    TH1F *h_const_1_alphaDown = new TH1F(*h_const_1);
    h_const_1_alphaDown->SetName("h_const_1_alphaDown");
    TH1F *h_const_10_alphaUp = new TH1F(*h_const_10);
    h_const_10_alphaUp->SetName("h_const_10_alphaUp");
    TH1F *h_const_10_alphaDown = new TH1F(*h_const_10);
    h_const_10_alphaDown->SetName("h_const_10_alphaDown");
    TH1F *h_const_100_alphaUp = new TH1F(*h_const_100);
    h_const_100_alphaUp->SetName("h_const_100_alphaUp");
    TH1F *h_const_100_alphaDown = new TH1F(*h_const_100);
    h_const_100_alphaDown->SetName("h_const_100_alphaDown");
    TH1F *h_const_0p05_alphaUp = new TH1F(*h_const_0p05);
    h_const_0p05_alphaUp->SetName("h_const_0p05_alphaUp");
    TH1F *h_const_0p05_alphaDown = new TH1F(*h_const_0p05);
    h_const_0p05_alphaDown->SetName("h_const_0p05_alphaDown");
    TH1F *h_const_0p5_alphaUp = new TH1F(*h_const_0p5);
    h_const_0p5_alphaUp->SetName("h_const_0p5_alphaUp");
    TH1F *h_const_0p5_alphaDown = new TH1F(*h_const_0p5);
    h_const_0p5_alphaDown->SetName("h_const_0p5_alphaDown");
    TH1F *h_const_5_alphaUp = new TH1F(*h_const_5);
    h_const_5_alphaUp->SetName("h_const_5_alphaUp");
    TH1F *h_const_5_alphaDown = new TH1F(*h_const_5);
    h_const_5_alphaDown->SetName("h_const_5_alphaDown");
    TH1F *h_const_50_alphaUp = new TH1F(*h_const_50);
    h_const_50_alphaUp->SetName("h_const_50_alphaUp");
    TH1F *h_const_50_alphaDown = new TH1F(*h_const_50);
    h_const_50_alphaDown->SetName("h_const_50_alphaDown");

    TH1F *h_const_0p01_sigmaUp = new TH1F(*h_const_0p01);
    h_const_0p01_sigmaUp->SetName("h_const_0p01_sigmaUp");
    TH1F *h_const_0p01_sigmaDown = new TH1F(*h_const_0p01);
    h_const_0p01_sigmaDown->SetName("h_const_0p01_sigmaDown");
    TH1F *h_const_0p1_sigmaUp = new TH1F(*h_const_0p1);
    h_const_0p1_sigmaUp->SetName("h_const_0p1_sigmaUp");
    TH1F *h_const_0p1_sigmaDown = new TH1F(*h_const_0p1);
    h_const_0p1_sigmaDown->SetName("h_const_0p1_sigmaDown");
    TH1F *h_const_1_sigmaUp = new TH1F(*h_const_1);
    h_const_1_sigmaUp->SetName("h_const_1_sigmaUp");
    TH1F *h_const_1_sigmaDown = new TH1F(*h_const_1);
    h_const_1_sigmaDown->SetName("h_const_1_sigmaDown");
    TH1F *h_const_10_sigmaUp = new TH1F(*h_const_10);
    h_const_10_sigmaUp->SetName("h_const_10_sigmaUp");
    TH1F *h_const_10_sigmaDown = new TH1F(*h_const_10);
    h_const_10_sigmaDown->SetName("h_const_10_sigmaDown");
    TH1F *h_const_100_sigmaUp = new TH1F(*h_const_100);
    h_const_100_sigmaUp->SetName("h_const_100_sigmaUp");
    TH1F *h_const_100_sigmaDown = new TH1F(*h_const_100);
    h_const_100_sigmaDown->SetName("h_const_100_sigmaDown");
    TH1F *h_const_0p05_sigmaUp = new TH1F(*h_const_0p05);
    h_const_0p05_sigmaUp->SetName("h_const_0p05_sigmaUp");
    TH1F *h_const_0p05_sigmaDown = new TH1F(*h_const_0p05);
    h_const_0p05_sigmaDown->SetName("h_const_0p05_sigmaDown");
    TH1F *h_const_0p5_sigmaUp = new TH1F(*h_const_0p5);
    h_const_0p5_sigmaUp->SetName("h_const_0p5_sigmaUp");
    TH1F *h_const_0p5_sigmaDown = new TH1F(*h_const_0p5);
    h_const_0p5_sigmaDown->SetName("h_const_0p5_sigmaDown");
    TH1F *h_const_5_sigmaUp = new TH1F(*h_const_5);
    h_const_5_sigmaUp->SetName("h_const_5_sigmaUp");
    TH1F *h_const_5_sigmaDown = new TH1F(*h_const_5);
    h_const_5_sigmaDown->SetName("h_const_5_sigmaDown");
    TH1F *h_const_50_sigmaUp = new TH1F(*h_const_50);
    h_const_50_sigmaUp->SetName("h_const_50_sigmaUp");
    TH1F *h_const_50_sigmaDown = new TH1F(*h_const_50);
    h_const_50_sigmaDown->SetName("h_const_50_sigmaDown");

    h_mu->Scale(((1.0)*(exp_count))/(mu_count)*(0.351485));
    h_e->Scale(((1.0)*(exp_count))/(e_count)*0.525415);
    h_other->Scale(((1.0)*(exp_count))/(other_count)*0.1231);
    h_air->Add(h_mu);
    h_air->Add(h_e);
    h_air->Add(h_other);

    h_const_0p01->Scale(((1.0)*(exp_count))/(dm1_count)*(0.351485));
    h_const_0p1->Scale(((1.0)*(exp_count))/(dm2_count)*(0.351485));
    h_const_1->Scale(((1.0)*(exp_count))/(dm3_count)*(0.351485));
    h_const_10->Scale(((1.0)*(exp_count))/(dm4_count)*(0.351485));
    h_const_100->Scale(((1.0)*(exp_count))/(dm5_count)*(0.351485));
    h_const_0p05->Scale(((1.0)*(exp_count))/(dm6_count)*(0.351485));
    h_const_0p5->Scale(((1.0)*(exp_count))/(dm7_count)*(0.351485));
    h_const_5->Scale(((1.0)*(exp_count))/(dm8_count)*(0.351485));
    h_const_50->Scale(((1.0)*(exp_count))/(dm9_count)*(0.351485));

    h_mu_alphaUp->Scale(((1.0)*(exp_count))/(mu_count)*(0.351485)*(1.148));
    h_e_alphaUp->Scale(((1.0)*(exp_count))/(e_count)*(0.525415)*(1.046));

    h_mu_alphaDown->Scale(((1.0)*(exp_count))/(mu_count)*(0.351485)*(0.852));
    h_e_alphaDown->Scale(((1.0)*(exp_count))/(e_count)*(0.525415)*(0.954));

    h_const_0p01_alphaUp->Scale(((1.0)*(exp_count))/(dm1_count)*(0.351485)*(1.148));
    h_const_0p1_alphaUp->Scale(((1.0)*(exp_count))/(dm2_count)*(0.351485)*(1.148));
    h_const_1_alphaUp->Scale(((1.0)*(exp_count))/(dm3_count)*(0.351485)*(1.148));
    h_const_10_alphaUp->Scale(((1.0)*(exp_count))/(dm4_count)*(0.351485)*(1.148));
    h_const_100_alphaUp->Scale(((1.0)*(exp_count))/(dm5_count)*(0.351485)*(1.148));
    h_const_0p05_alphaUp->Scale(((1.0)*(exp_count))/(dm6_count)*(0.351485)*(1.148));
    h_const_0p5_alphaUp->Scale(((1.0)*(exp_count))/(dm7_count)*(0.351485)*(1.148));
    h_const_5_alphaUp->Scale(((1.0)*(exp_count))/(dm8_count)*(0.351485)*(1.148));
    h_const_50_alphaUp->Scale(((1.0)*(exp_count))/(dm9_count)*(0.351485)*(1.148));

    h_const_0p01_alphaDown->Scale(((1.0)*(exp_count))/(dm1_count)*(0.351485)*(0.852));
    h_const_0p1_alphaDown->Scale(((1.0)*(exp_count))/(dm2_count)*(0.351485)*(0.852));
    h_const_1_alphaDown->Scale(((1.0)*(exp_count))/(dm3_count)*(0.351485)*(0.852));
    h_const_10_alphaDown->Scale(((1.0)*(exp_count))/(dm4_count)*(0.351485)*(0.852));
    h_const_100_alphaDown->Scale(((1.0)*(exp_count))/(dm5_count)*(0.351485)*(0.852));
    h_const_0p05_alphaDown->Scale(((1.0)*(exp_count))/(dm6_count)*(0.351485)*(0.852));
    h_const_0p5_alphaDown->Scale(((1.0)*(exp_count))/(dm7_count)*(0.351485)*(0.852));
    h_const_5_alphaDown->Scale(((1.0)*(exp_count))/(dm8_count)*(0.351485)*(0.852));
    h_const_50_alphaDown->Scale(((1.0)*(exp_count))/(dm9_count)*(0.351485)*(0.852));

    h_mu_sigmaUp->Scale(((1.0)*(exp_count))/(mu_count)*(0.351485+0.0049));
    h_e_sigmaUp->Scale(((1.0)*(exp_count))/(e_count)*(0.525415+0.0049));
    h_other_sigmaUp->Scale(((1.0)*(exp_count))/(other_count)*(0.1231+0.0019));

    h_mu_sigmaDown->Scale(((1.0)*(exp_count))/(mu_count)*(0.351485-0.0049));
    h_e_sigmaDown->Scale(((1.0)*(exp_count))/(e_count)*(0.525415-0.0049));
    h_other_sigmaDown->Scale(((1.0)*(exp_count))/(other_count)*(0.1231-0.0019));

    h_const_0p01_sigmaUp->Scale(((1.0)*(exp_count))/(dm1_count)*(0.351485+0.0049));
    h_const_0p1_sigmaUp->Scale(((1.0)*(exp_count))/(dm2_count)*(0.351485+0.0049));
    h_const_1_sigmaUp->Scale(((1.0)*(exp_count))/(dm3_count)*(0.351485+0.0049));
    h_const_10_sigmaUp->Scale(((1.0)*(exp_count))/(dm4_count)*(0.351485+0.0049));
    h_const_100_sigmaUp->Scale(((1.0)*(exp_count))/(dm5_count)*(0.351485+0.0049));
    h_const_0p05_sigmaUp->Scale(((1.0)*(exp_count))/(dm6_count)*(0.351485+0.0049));
    h_const_0p5_sigmaUp->Scale(((1.0)*(exp_count))/(dm7_count)*(0.351485+0.0049));
    h_const_5_sigmaUp->Scale(((1.0)*(exp_count))/(dm8_count)*(0.351485+0.0049));
    h_const_50_sigmaUp->Scale(((1.0)*(exp_count))/(dm9_count)*(0.351485+0.0049));


    h_const_0p01_sigmaDown->Scale(((1.0)*(exp_count))/(dm1_count)*(0.351485-0.0049));
    h_const_0p1_sigmaDown->Scale(((1.0)*(exp_count))/(dm2_count)*(0.351485-0.0049));
    h_const_1_sigmaDown->Scale(((1.0)*(exp_count))/(dm3_count)*(0.351485-0.0049));
    h_const_10_sigmaDown->Scale(((1.0)*(exp_count))/(dm4_count)*(0.351485-0.0049));
    h_const_100_sigmaDown->Scale(((1.0)*(exp_count))/(dm5_count)*(0.351485-0.0049));
    h_const_0p05_sigmaDown->Scale(((1.0)*(exp_count))/(dm6_count)*(0.351485-0.0049));
    h_const_0p5_sigmaDown->Scale(((1.0)*(exp_count))/(dm7_count)*(0.351485-0.0049));
    h_const_5_sigmaDown->Scale(((1.0)*(exp_count))/(dm8_count)*(0.351485-0.0049));
    h_const_50_sigmaDown->Scale(((1.0)*(exp_count))/(dm9_count)*(0.351485-0.0049));

    cout<<"mu_fraction="<<((1.0)*(h_mu->Integral()))/((1.0)*(h_air->Integral()))*100.0<<endl;
    cout<<"e_fraction="<<((1.0)*(h_e->Integral()))/((1.0)*(h_air->Integral()))*100.0<<endl;
    cout<<"other_fraction="<<((1.0)*(h_other->Integral()))/((1.0)*(h_air->Integral()))*100.0<<endl;

    

    TFile *file = new TFile("Hist.root", "RECREATE");
    

    h_const_0p01->Write();
    h_const_0p1->Write();
    h_const_1->Write();
    h_const_10->Write();
    h_const_100->Write();
    h_const_0p05->Write();
    h_const_0p5->Write();
    h_const_5->Write();
    h_const_50->Write();

    h_mu->Write();
    h_e->Write();
    h_other->Write();

    h_const_0p01_alphaUp->Write();
    h_const_0p1_alphaUp->Write();
    h_const_1_alphaUp->Write();
    h_const_10_alphaUp->Write();
    h_const_100_alphaUp->Write();
    h_const_0p05_alphaUp->Write();
    h_const_0p5_alphaUp->Write();
    h_const_5_alphaUp->Write();
    h_const_50_alphaUp->Write();

    h_mu_alphaUp->Write();
    h_e_alphaUp->Write();

    h_const_0p01_alphaDown->Write();
    h_const_0p1_alphaDown->Write();
    h_const_1_alphaDown->Write();
    h_const_10_alphaDown->Write();
    h_const_100_alphaDown->Write();
    h_const_0p05_alphaDown->Write();
    h_const_0p5_alphaDown->Write();
    h_const_5_alphaDown->Write();
    h_const_50_alphaDown->Write();

    h_mu_alphaDown->Write();
    h_e_alphaDown->Write();

    h_const_0p01_sigmaUp->Write();
    h_const_0p1_sigmaUp->Write();
    h_const_1_sigmaUp->Write();
    h_const_10_sigmaUp->Write();
    h_const_100_sigmaUp->Write();
    h_const_0p05_sigmaUp->Write();
    h_const_0p5_sigmaUp->Write();
    h_const_5_sigmaUp->Write();
    h_const_50_sigmaUp->Write();

    h_mu_sigmaUp->Write();
    h_e_sigmaUp->Write();
    h_other_sigmaUp->Write();

    h_const_0p01_sigmaDown->Write();
    h_const_0p1_sigmaDown->Write();
    h_const_1_sigmaDown->Write();
    h_const_10_sigmaDown->Write();
    h_const_100_sigmaDown->Write();
    h_const_0p05_sigmaDown->Write();
    h_const_0p5_sigmaDown->Write();
    h_const_5_sigmaDown->Write();
    h_const_50_sigmaDown->Write();

    h_mu_sigmaDown->Write();
    h_e_sigmaDown->Write();
    h_other_sigmaDown->Write();
    h_obs->Write();

    file->Close();

    gPad->SetLogy();
	gStyle->SetOptStat(0);
	h_obs->SetMinimum(1e-1);
    h_obs->SetMaximum(1e6);
	h_obs->GetYaxis()->SetTitle("relative events");
	h_obs->GetXaxis()->SetTitle("#theta(rad)");
	h_obs->Draw("E");
	h_obs->SetLineColor(kBlack);
    h_obs->SetLineWidth(2);

	h_air->SetLineColor(kRed);
    h_air->SetLineWidth(2);
	h_air->Draw("HIST same");

    h_const_0p01->SetLineColor(kBlue);
    h_const_0p01->SetLineWidth(2);
    h_const_0p01->Draw("HIST same");

    h_const_0p1->SetLineColor(kCyan);
    h_const_0p1->SetLineWidth(2);
    h_const_0p1->Draw("HIST same");

    h_const_1->SetLineColor(kGreen);
    h_const_1->SetLineWidth(2);
    h_const_1->Draw("HIST same");

    h_const_10->SetLineColor(kMagenta);
    h_const_10->SetLineWidth(2);
    h_const_10->Draw("HIST same");

    h_const_100->SetLineColor(kOrange);
    h_const_100->SetLineWidth(2);
    h_const_100->Draw("HIST same");

    h_const_0p05->SetLineColor(kYellow);
    h_const_0p05->SetLineWidth(2);
    h_const_0p05->Draw("HIST same");

    h_const_0p5->SetLineColor(kPink);
    h_const_0p5->SetLineWidth(2);
    h_const_0p5->Draw("HIST same");

    h_const_5->SetLineColor(kSpring);
    h_const_5->SetLineWidth(2);
    h_const_5->Draw("HIST same");

    h_const_50->SetLineColor(kTeal);
    h_const_50->SetLineWidth(2);
    h_const_50->Draw("HIST same");

    TLegend *legend = new TLegend(0.4, 0.7, 0.9, 0.9);
	legend->SetX1NDC(0.6);
    legend->SetY1NDC(0.7);
    legend->SetX2NDC(0.9);
    legend->SetY2NDC(0.9);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_obs, "observed", "l");
    legend->AddEntry(h_air, "air", "l");
    legend->AddEntry(h_const_0p01, "1e-02GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_0p1, "1e-01GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_1, "1e+00GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_10, "1e+01GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_100, "1e+02GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_0p05, "5e-02GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_0p5, "5e-01GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_5, "5e+00GeVDM-mu_expected", "l");
    legend->AddEntry(h_const_50, "5e+01GeVDM-mu_expected", "l");
    legend->Draw();

    return exp_count;
}