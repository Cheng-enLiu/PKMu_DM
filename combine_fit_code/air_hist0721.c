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

#include <cstdlib>
#include <ctime>

Int_t air_hist0721(TString expROOTfilename="250212_63d_air.root",
    TString simROOTfilename="smeared_rec_CryMuAna_0_749.root",
    TString dm4ROOTfilename="smeared_rec_poca_DM_1GeV_0722.root")
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

    TCanvas *canvas = new TCanvas("canvas", "ang", 800, 600);

    TH1F *h_air = new TH1F("h_air", "distribution of scattering angle #it{#theta}", 50, 0.05, 0.5);
    TH1F *h_obs = new TH1F("h_obs", "distribution of scattering angle #it{#theta}", 50, 0.05, 0.5);
    TH1F *h_mu = new TH1F("h_mu", "distribution of scattering angle #it{#theta}", 50, 0.05, 0.5);
    TH1F *h_e = new TH1F("h_e", "distribution of scattering angle #it{#theta}", 50, 0.05, 0.5);
    TH1F *h_other = new TH1F("h_other", "distribution of scattering angle #it{#theta}", 50, 0.05, 0.5);

    TH1F *h_const_1 = new TH1F("h_const_1", "h_const_1", 50, 0.05, 0.5);
    TH1F *h_mu_1 = new TH1F("h_mu_1", "h_mu_1", 50, 0.05, 0.5);

    int exp_count=0;
    int sim_count=0;
    int mu_count=0;
    int e_count=0;
    int other_count=0;

    int dm4_count=0;
    int dm4_counts=0;

    int i;

    for(i=0;i<nentries_exp;i++){
        Trec_fix->GetEntry(i);
	    if(x1>145||x1<-145||x2>145||x2<-145||x3>145||x3<-145||x4>145||x4<-145||y1>145||y1<-145||y2>145||y2<-145||y3>145||y3<-145||y4>145||y4<-145) continue;
        if (times[1]<0 || times[10]<0 || times[11]<0 || times[12]<0 || times[13]<0 || times[14]<0 || times[15]<0 || times[16]<0 || times[17]<0 || times[18]<0 || times[19]<0 ||
            times[20]<0 || times[21]<0 || times[22]<0 || times[23]<0 || times[24]<0 || times[25]<0 || times[28]<0 || times[29]<0 || times[30]<0 || times[31]<0) 
            continue;
        if(z>110||z<-110) continue;
        if(x>110||x<-110) continue;
        if(y>110||y<-110) continue;
        if(ang>0.5||ang<0.05) continue;
        exp_count++;
        h_obs->Fill(ang);
    }

    for(i=0;i<nentries_sim;i++){
        Trec->GetEntry(i);
        if(x1s>145||x1s<-145||x2s>145||x2s<-145||x3s>145||x3s<-145||x4s>145||x4s<-145||y1s>145||y1s<-145||y2s>145||y2s<-145||y3s>145||y3s<-145||y4s>145||y4s<-145) continue;
        if(zs>110||zs<-110) continue;
        if(xs>110||xs<-110) continue;
        if(ys>110||ys<-110) continue;
        if(angs>0.5||angs<0.05) continue;
        sim_count++;
        //h_air->Fill(angs);
        if(abs(Pid)==11){
            e_count++;
            h_e->Fill(angs);
        }
        if(abs(Pid)==13){
            mu_count++;
            h_mu->Fill(angs);
        }
        if(abs(Pid)!=11&&abs(Pid)!=13){
            other_count++;
            h_other->Fill(angs);
        }
    }

    for(i=0;i<nentries_dm4;i++){
        Trec_dm4->GetEntry(i);
	    if(dm4x1>145||dm4x1<-145||dm4x2>145||dm4x2<-145||dm4x3>145||dm4x3<-145||dm4x4>145||dm4x4<-145||dm4y1>145||dm4y1<-145||dm4y2>145||dm4y2<-145||dm4y3>145||dm4y3<-145||dm4y4>145||dm4y4<-145) continue;
        if(dm4z>110||dm4z<-110) continue;
        if(dm4x>110||dm4x<-110) continue;
        if(dm4y>110||dm4y<-110) continue;
        if(dm4ang>0.5||dm4ang<0.05) continue;
        if(abs(dm4Pid)==13){
            h_mu_1->Fill(dm4ang);
            dm4_count++;
            if(dm4Mark==1){
                dm4_counts++;
                h_const_1->Fill(dm4ang);
            }
        }
    }




    cout<<"exp_count="<<exp_count<<endl;
    cout<<"sim_count="<<sim_count<<endl;
    cout<<"mu="<<mu_count<<endl;
    cout<<"e="<<e_count<<endl;
    cout<<"other="<<other_count<<endl;

    cout<<"dm4_count="<<dm4_count<<endl;
    cout<<"dm4_counts="<<dm4_counts<<endl;

    double sigma_mu = (100 * 1.0) * pow((((mu_count * 1.0) * ((mu_count * 1.0) + (sim_count * 1.0))) / (pow((sim_count * 1.0) , 3))) , 0.5);
    double sigma_e = (100 * 1.0) * pow((((e_count * 1.0) * ((e_count * 1.0) + (sim_count * 1.0))) / (pow((sim_count * 1.0) , 3))) , 0.5);
    double sigma_other = (100 * 1.0) * pow((((other_count * 1.0) * ((other_count * 1.0) + (sim_count * 1.0))) / (pow((sim_count * 1.0) , 3))) , 0.5);
    cout<<"sigma_mu="<<sigma_mu<<endl;
    cout<<"sigma_e="<<sigma_e<<endl;
    cout<<"sigma_other"<<sigma_other<<endl;

    //归一化

    //h_air->Scale((exp_count*1.0)/sim_count);
    /*h_obs->Scale((1.0/h_obs->Integral()));
    h_mu->Scale(((1.0)/sim_count));
    h_e->Scale(((1.0)/sim_count));
    h_other->Scale(((1.0)/sim_count));*/

    TH1F *h_air1 = new TH1F(*h_air);
    TH1F *h_mu1 = new TH1F(*h_mu);
    TH1F *h_e1 = new TH1F(*h_e);
    TH1F *h_other1 = new TH1F(*h_other);

    TH1F *h_air2 = new TH1F(*h_air);
    TH1F *h_mu2 = new TH1F(*h_mu);
    TH1F *h_e2 = new TH1F(*h_e);
    TH1F *h_other2 = new TH1F(*h_other);

    TH1F *h_air3 = new TH1F(*h_air);
    TH1F *h_mu3 = new TH1F(*h_mu);
    TH1F *h_e3 = new TH1F(*h_e);
    TH1F *h_other3 = new TH1F(*h_other);

    TH1F *h_air4 = new TH1F(*h_air);
    TH1F *h_mu4 = new TH1F(*h_mu);
    TH1F *h_e4 = new TH1F(*h_e);
    TH1F *h_other4 = new TH1F(*h_other);

    h_obs->Scale(((1.0)/h_obs->Integral()));
    h_mu->Scale(((1.0)/h_mu->Integral())*0.351485);
    h_e->Scale(((1.0)/h_e->Integral())*0.525415);
    h_other->Scale(((1.0)/h_other->Integral())*0.1231);

    /*h_obs->Scale(((1.0)/h_obs->Integral()));
    h_mu->Scale(((1.0)/h_mu->Integral()));
    h_e->Scale(((1.0)/h_e->Integral()));
    h_other->Scale(((1.0)/h_other->Integral()));*/

    h_air->Add(h_mu);
    h_air->Add(h_e);
    h_air->Add(h_other);

    //h_const_1->Scale(((1.0)*(exp_count))/(dm4_count)*0.351485);
    h_const_1->Scale((1.0)/(h_const_1->Integral()));

    h_mu1->Scale(((1.0)/h_mu1->Integral())*((0.351485)+(0.05225)));
    h_e1->Scale(((1.0)/h_e1->Integral())*((0.525415)+(0.02466)));
    h_other1->Scale(((1.0)/h_other1->Integral())*((0.1231)+(0.0019)));
    h_air1->Add(h_mu1);
    h_air1->Add(h_e1);
    h_air1->Add(h_other1);

    h_mu2->Scale(((1.0)/h_mu2->Integral())*((0.351485)-(0.05225)));
    h_e2->Scale(((1.0)/h_e2->Integral())*((0.525415)-(0.02466)));
    h_other2->Scale(((1.0)/h_other2->Integral())*((0.1231)-(0.0019)));
    h_air2->Add(h_mu2);
    h_air2->Add(h_e2);
    h_air2->Add(h_other2);
    
    h_mu3->Scale(((1.0)/h_mu3->Integral())*((0.351485)+2*(0.05225)));
    h_e3->Scale(((1.0)/h_e3->Integral())*((0.525415)+2*(0.12523)));
    h_other3->Scale(((1.0)/h_other3->Integral())*((0.1231)+2*(0.0019)));
    h_air3->Add(h_mu3);
    h_air3->Add(h_e3);
    h_air3->Add(h_other3);

    h_mu4->Scale(((1.0)/h_mu4->Integral())*((0.351485)-2*(0.05225)));
    h_e4->Scale(((1.0)/h_e4->Integral())*((0.525415)-2*(0.02466)));
    h_other4->Scale(((1.0)/h_other4->Integral())*((0.1231)-2*(0.0019)));
    h_air4->Add(h_mu4);
    h_air4->Add(h_e4);
    h_air4->Add(h_other4);
    
    /*TFile *file = new TFile("hist_of_air.root", "RECREATE");
    
    h_obs->Write();
    //h_air->Write();
    h_mu->Write();
    h_e->Write();
    h_other->Write();
    h_others->Write();

    file->Close();*/

    TH1F *h_air_copy = new TH1F("h_air_copy", "distribution of scattering angle #it{#theta}", 10, 0.05, 0.5);
    TH1F *h_obs_copy = new TH1F("h_obs_copy", "distribution of scattering angle #it{#theta}", 10, 0.05, 0.5);
    TH1F *h_mu_copy = new TH1F("h_mu_copy", "distribution of scattering angle #it{#theta}", 10, 0.05, 0.5);
    TH1F *h_e_copy = new TH1F("h_e_copy", "distribution of scattering angle #it{#theta}", 10, 0.05, 0.5);
    TH1F *h_other_copy = new TH1F("h_other_copy", "distribution of scattering angle #it{#theta}", 10, 0.05, 0.5);

    int exp_count_copy=0;
    int sim_count_copy=0;
    int mu_count_copy=0;
    int e_count_copy=0;
    int other_count_copy=0;

    for(i=0;i<nentries_exp;i++){
        Trec_fix->GetEntry(i);
	    if(x1>145||x1<-145||x2>145||x2<-145||x3>145||x3<-145||x4>145||x4<-145||y1>145||y1<-145||y2>145||y2<-145||y3>145||y3<-145||y4>145||y4<-145) continue;
        if (times[1]<0 || times[10]<0 || times[11]<0 || times[12]<0 || times[13]<0 || times[14]<0 || times[15]<0 || times[16]<0 || times[17]<0 || times[18]<0 || times[19]<0 ||
            times[20]<0 || times[21]<0 || times[22]<0 || times[23]<0 || times[24]<0 || times[25]<0 || times[28]<0 || times[29]<0 || times[30]<0 || times[31]<0) 
            continue;
        if(z>110||z<-110) continue;
        if(x>110||x<-110) continue;
        if(y>110||y<-110) continue;
        if(ang>0.5||ang<0.05) continue;
        exp_count_copy++;
        h_obs_copy->Fill(ang);
    }

    for(i=0;i<nentries_sim;i++){
        Trec->GetEntry(i);
        if(x1s>145||x1s<-145||x2s>145||x2s<-145||x3s>145||x3s<-145||x4s>145||x4s<-145||y1s>145||y1s<-145||y2s>145||y2s<-145||y3s>145||y3s<-145||y4s>145||y4s<-145) continue;
        if(zs>110||zs<-110) continue;
        if(xs>110||xs<-110) continue;
        if(ys>110||ys<-110) continue;
        if(angs>0.5||angs<0.05) continue;
        sim_count_copy++;
        if(abs(Pid)==11){
            e_count_copy++;
            h_e_copy->Fill(angs);
        }
        if(abs(Pid)==13){
            mu_count_copy++;
            h_mu_copy->Fill(angs);
        }
        if(abs(Pid)!=11&&abs(Pid)!=13){
            other_count_copy++;
            h_other_copy->Fill(angs);
        }
    }

    TH1F *h_air1_copy = new TH1F(*h_air_copy);
    TH1F *h_mu1_copy = new TH1F(*h_mu_copy);
    TH1F *h_e1_copy = new TH1F(*h_e_copy);
    TH1F *h_other1_copy = new TH1F(*h_other_copy);

    TH1F *h_air2_copy = new TH1F(*h_air_copy);
    TH1F *h_mu2_copy = new TH1F(*h_mu_copy);
    TH1F *h_e2_copy = new TH1F(*h_e_copy);
    TH1F *h_other2_copy = new TH1F(*h_other_copy);

    TH1F *h_air3_copy = new TH1F(*h_air_copy);
    TH1F *h_mu3_copy = new TH1F(*h_mu_copy);
    TH1F *h_e3_copy = new TH1F(*h_e_copy);
    TH1F *h_other3_copy = new TH1F(*h_other_copy);

    TH1F *h_air4_copy = new TH1F(*h_air_copy);
    TH1F *h_mu4_copy = new TH1F(*h_mu_copy);
    TH1F *h_e4_copy = new TH1F(*h_e_copy);
    TH1F *h_other4_copy = new TH1F(*h_other_copy);

    h_obs_copy->Scale((1.0/h_obs_copy->Integral()));
    h_mu_copy->Scale(((1.0)/h_mu_copy->Integral())*0.351485);
    h_e_copy->Scale(((1.0)/h_e_copy->Integral())*0.525415);
    h_other_copy->Scale(((1.0)/h_other_copy->Integral())*0.1231);

    /*h_obs_copy->Scale((1.0/h_obs_copy->Integral()));
    h_mu_copy->Scale(((1.0)/sim_count_copy));
    h_e_copy->Scale(((1.0)/sim_count_copy));
    h_other_copy->Scale(((1.0)/sim_count_copy));*/

    h_air_copy->Add(h_mu_copy);
    h_air_copy->Add(h_e_copy);
    h_air_copy->Add(h_other_copy);

    h_mu1_copy->Scale(((1.0)/h_mu1_copy->Integral())*((0.351485)+(0.05225)));
    h_e1_copy->Scale(((1.0)/h_e1_copy->Integral())*((0.525415)+(0.02466)));
    h_other1_copy->Scale(((1.0)/h_other1_copy->Integral())*((0.1231)+(0.0019)));
    h_air1_copy->Add(h_mu1_copy);
    h_air1_copy->Add(h_e1_copy);
    h_air1_copy->Add(h_other1_copy);

    h_mu2_copy->Scale(((1.0)/h_mu2_copy->Integral())*((0.351485)-(0.05225)));
    h_e2_copy->Scale(((1.0)/h_e2_copy->Integral())*((0.525415)-(0.02466)));
    h_other2_copy->Scale(((1.0)/h_other2_copy->Integral())*((0.1231)-(0.0019)));
    h_air2_copy->Add(h_mu2_copy);
    h_air2_copy->Add(h_e2_copy);
    h_air2_copy->Add(h_other2_copy);

    TPad* pad_top = new TPad("pad_top", "", 0, 0.35, 1, 1.0);
    TPad* pad_bot = new TPad("pad_bot", "", 0, 0.0, 1, 0.35);

    // 设置边距，使 pad 紧贴，不显示上图X轴
    pad_top->SetBottomMargin(0.0);
    pad_bot->SetTopMargin(0.0);
    pad_bot->SetBottomMargin(0.3);

    // 绘制 pad
    pad_top->Draw();
    pad_bot->Draw();

    pad_top->cd();
    gPad->SetLogy();
	gStyle->SetOptStat(0);
    h_obs->SetStats(0);
    h_obs->SetMinimum(3e-4);
    h_obs->SetMaximum(2e-1);
    h_obs->SetTitle("");
	h_obs->GetYaxis()->SetTitle("relative events");
    h_obs->GetXaxis()->SetLabelSize(0);     // 隐藏X轴标签
    h_obs->GetXaxis()->SetTitleSize(0);     // 隐藏X轴标题
    h_obs->GetYaxis()->SetLabelSize(0.05);     // Y轴标签
    h_obs->GetYaxis()->SetTitleSize(0.05);     // Y轴标题
    h_obs->SetMarkerStyle(20);
    h_obs->SetMarkerSize(1.2);
	h_obs->Draw("E");
	h_obs->SetLineColor(kBlack); //黑色
    h_obs->SetLineWidth(4);

    int N = h_air->GetNbinsX();
    TGraphAsymmErrors* gr_band1 = new TGraphAsymmErrors(N + 2);  // 多两个点
    TGraphAsymmErrors* gr_band2 = new TGraphAsymmErrors(N + 2);  // 多两个点

    double binWidth = h_air->GetBinWidth(1); // 假设等宽

    for (int i = 1; i <= N; ++i) {
        double x = h_air->GetBinCenter(i);
        double y = h_air->GetBinContent(i);
        double y_up1 = h_air1->GetBinContent(i);
        double y_down1 = h_air2->GetBinContent(i);
        double y_up2 = h_air3->GetBinContent(i);
        double y_down2 = h_air4->GetBinContent(i);

        double err_up1 = y_up1 - y;
        double err_down1 = y - y_down1;
        double err_up2 = y_up2 - y;
        double err_down2 = y - y_down2;

        gr_band1->SetPoint(i, x, y);
        gr_band1->SetPointError(i, 0, 0, err_down1, err_up1);
        gr_band2->SetPoint(i, x, y);
        gr_band2->SetPointError(i, 0, 0, err_down2, err_up2);
    }

    // 添加第一个点的延伸
    {
        double x = h_air->GetBinCenter(1) - binWidth;
        double y = h_air->GetBinContent(1);
        double y_up1 = h_air1->GetBinContent(1);
        double y_down1 = h_air2->GetBinContent(1);
        double y_up2 = h_air3->GetBinContent(1);
        double y_down2 = h_air4->GetBinContent(1);
        double err_up1 = y_up1 - y;
        double err_down1 = y - y_down1;
        double err_up2 = y_up2 - y;
        double err_down2 = y - y_down2;

        gr_band1->SetPoint(0, x, y);
        gr_band1->SetPointError(0, 0, 0, err_down1, err_up1);
        gr_band2->SetPoint(0, x, y);
        gr_band2->SetPointError(0, 0, 0, err_down2, err_up2);
    }

    // 添加最后一个点的延伸
    {
        double x = h_air->GetBinCenter(N) + binWidth;
        double y = h_air->GetBinContent(N);
        double y_up1 = h_air1->GetBinContent(N);
        double y_down1 = h_air2->GetBinContent(N);
        double y_up2 = h_air3->GetBinContent(N);
        double y_down2 = h_air4->GetBinContent(N);
        double err_up1 = y_up1 - y;
        double err_down1 = y - y_down1;
        double err_up2 = y_up2 - y;
        double err_down2 = y - y_down2;

        gr_band1->SetPoint(N + 1, x, y);
        gr_band1->SetPointError(N + 1, 0, 0, err_down1, err_up1);
        gr_band2->SetPoint(N + 1, x, y);
        gr_band2->SetPointError(N + 1, 0, 0, err_down2, err_up2);
    }

    gr_band1->SetFillColor(kYellow); // 填充颜色
    gr_band1->SetFillStyle(1001);    // 半透明填充
    gr_band2->SetFillColor(kGreen); // 填充颜色
    gr_band2->SetFillStyle(1001);    // 半透明填充

    //gr_band2->Draw("3 SAME");
    gr_band1->Draw("3 SAME");
    gPad->RedrawAxis();

	h_air->SetLineColor(kRed); //红色
    h_air->SetLineWidth(3);
	h_air->Draw("HIST same");
    h_obs->Draw("E same");


    h_mu->SetLineColor(kBlue); //蓝色
    h_mu->SetLineWidth(3);
	h_mu->Draw("HIST same"); 

    h_e->SetLineColor(kCyan); 
    h_e->SetLineWidth(3);
	h_e->Draw("HIST same");
    
    
    h_other->SetLineColor(kMagenta); 
    h_other->SetLineWidth(3);
	h_other->Draw("HIST same");

    h_const_1->SetLineColor(kOrange);
    h_const_1->SetLineWidth(3);
    h_const_1->Draw("HIST same");



    TLegend *legend = new TLegend(0.20, 0.65, 0.9, 0.9);
    legend->SetNColumns(2); 
	legend->SetX1NDC(0.6);
    legend->SetY1NDC(0.70);
    legend->SetX2NDC(0.9);
    legend->SetY2NDC(0.9);
    legend->SetTextSize(0.05);
    legend->SetBorderSize(0);   // 去除边框
    legend->SetFillStyle(0); 
    legend->AddEntry(h_obs, "observed", "lep");
    legend->AddEntry(h_air, "MC_all_particle", "l");
    legend->AddEntry(gr_band1, "#pm 1 std. deviation");
    //legend->AddEntry(gr_band2, "#pm 2 std. deviation");
    legend->AddEntry(h_mu, Form("MC_#mu (%.1f%% #pm 5.2%%)", (h_mu->Integral() * 1.0)/(h_air->Integral() * 1.0)*100.0 ), "l");
    legend->AddEntry(h_e, Form("MC_e (%.1f%% #pm 2.5%%)", (h_e->Integral() * 1.0)/(h_air->Integral() * 1.0)*100.0), "l");
    legend->AddEntry(h_other, Form("MC_other (%.1f%% #pm 0.1%%)", (h_other->Integral() * 1.0)/(h_air->Integral() * 1.0)*100.0), "l");
    // 第7行左边：空占位项
    TH1D *h_dummy = new TH1D("h_dummy", "", 1, 0, 1); // dummy histogram
    h_dummy->SetLineColor(0);
    h_dummy->SetLineWidth(0);
    legend->AddEntry(h_dummy, "", "");  // 占左栏位置
    // 第7行右边：你要的 entry
    legend->AddEntry(h_const_1, "MC 1 GeV DM-#mu scattering angle", "l");
    legend->Draw();

    /*TLegend *legend2 = new TLegend(0.57, 0.60, 0.85, 0.66);  // 手动调整位置和宽度
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.05);
    legend2->AddEntry(h_const_1, "1 GeV DM-#mu scattering angle", "l");
    legend2->Draw();*/

    pad_bot->cd();
    gStyle->SetOptStat(0);
    // 检查Bin的一致性
    if (h_obs_copy->GetNbinsX() != h_air_copy->GetNbinsX()) {
        std::cerr << "Error: Binning of h_obs_copy and h_air_copy do not match!" << std::endl;
        return 0;
    }

    // 创建一个新的直方图用于存储比值
    TH1F* h_ratio = (TH1F*)h_obs_copy->Clone("h_ratio");
    TH1F* h_ratio1 = (TH1F*)h_obs_copy->Clone("h_ratio1");
    TH1F* h_ratio2 = (TH1F*)h_obs_copy->Clone("h_ratio2");
    h_ratio->SetTitle("");
    h_ratio->SetMaximum(1.44);
    h_ratio->SetMinimum(0.66);
    h_ratio->GetXaxis()->SetTitle("#it{#theta} (rad)");
    h_ratio->GetYaxis()->SetLabelSize(0.08);     // Y轴标签
    h_ratio->GetYaxis()->SetTitleSize(0.08);     // Y轴标题
    h_ratio->GetYaxis()->SetTitle("observed/MC");
    h_ratio->GetYaxis()->SetTitleOffset(0.7);
    h_ratio->GetXaxis()->SetLabelSize(0.08);     // X轴标签
    h_ratio->GetXaxis()->SetTitleSize(0.08);     // X轴标题
    h_ratio->SetStats(0);

    // 计算比值并处理误差
    h_ratio->Divide(h_air_copy);
    h_ratio1->Divide(h_air1_copy);
    h_ratio2->Divide(h_air2_copy);

    int N_copy = h_ratio->GetNbinsX();
    TGraphAsymmErrors* gr_band_copy = new TGraphAsymmErrors(N_copy + 2);  // 多两个点

    double binWidth_copy = h_ratio->GetBinWidth(1); // 假设等宽

    for (int i = 1; i <= N_copy; ++i) {
        double x = h_ratio->GetBinCenter(i);
        double y = h_ratio->GetBinContent(i);
        double y_up = h_ratio1->GetBinContent(i);
        double y_down = h_ratio2->GetBinContent(i);

        double err_up = y_up - y;
        double err_down = y - y_down;

        gr_band_copy->SetPoint(i, x, y);
        gr_band_copy->SetPointError(i, 0, 0, err_down, err_up);
    }

    // 添加第一个点的延伸
    {
        double x = h_ratio->GetBinCenter(1) - binWidth_copy;
        double y = h_ratio->GetBinContent(1);
        double y_up = h_ratio1->GetBinContent(1);
        double y_down = h_ratio2->GetBinContent(1);
        double err_up = y_up - y;
        double err_down = y - y_down;

        gr_band_copy->SetPoint(0, x, y);
        gr_band_copy->SetPointError(0, 0, 0, err_down, err_up);
    }

    // 添加最后一个点的延伸
    {
        double x = h_ratio->GetBinCenter(N_copy) + binWidth_copy;
        double y = h_ratio->GetBinContent(N_copy);
        double y_up = h_ratio1->GetBinContent(N_copy);
        double y_down = h_ratio2->GetBinContent(N_copy);
        double err_up = y_up - y;
        double err_down = y - y_down;

        gr_band_copy->SetPoint(N_copy + 1, x, y);
        gr_band_copy->SetPointError(N_copy + 1, 0, 0, err_down, err_up);
    }

    TF1* fit_ratio = new TF1("fit_ratio", "pol0", h_ratio->GetXaxis()->GetXmin(), h_ratio->GetXaxis()->GetXmax());

    // 对 h_ratio 进行拟合
    h_ratio->Fit(fit_ratio, "R");

    // 设置绘图样式
    h_ratio->SetMarkerStyle(20);   // 数据点样式
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->SetMarkerSize(1.5);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetLineWidth(4);
    fit_ratio->SetLineColor(kRed); // 拟合线颜色
    fit_ratio->SetLineWidth(4);
    //fit_ratio->SetLineStyle(2);    // 虚线样式

    // 获取拟合常数值和误差
    double constant = fit_ratio->GetParameter(0);  // 获取常数参数值
    double constant_error = fit_ratio->GetParError(0);  // 获取常数参数的误差

    gr_band_copy->SetFillColor(kYellow); // 填充颜色
    gr_band_copy->SetFillStyle(1001);    // 半透明填充
    //gr_band->SetFillColorAlpha(kGray + 1, 0.35); // 灰色带，有透明度
    //gr_band->SetLineColor(0);                    // 去掉边框

    // 绘制图像
    h_ratio->Draw("E"); // 绘制数据点，带误差条
    // 绘制填充区域
    gr_band_copy->Draw("3 SAME");                          // "3" 表示绘制填充区域
    h_ratio->Draw("E SAME"); // 绘制数据点，带误差条
    fit_ratio->Draw("SAME"); // 在同一张图上绘制拟合曲线

    gPad->RedrawAxis();

    // 绘制拟合常数及误差
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.08);
    latex->SetTextAlign(13); // 左对齐
    latex->DrawLatex(0.67, 0.87, Form("Constant = %.3f #pm %.3f", constant, constant_error));
    //latex->DrawLatex(0.45,0.99, "observed/MC Ratio");

    // 添加图例
    TLegend* legend1 = new TLegend(0.2, 0.9, 0.9, 1.0);
    legend1->SetNColumns(3); 
    legend1->SetTextSize(0.08);
    legend1->AddEntry(h_ratio, "observed/MC Ratio", "lep");
    legend1->AddEntry(fit_ratio, "Fit (constant)", "l");
    legend1->AddEntry(gr_band_copy, "#pm 1 std. deviation");
    legend1->Draw();

    double xmin = h_obs->GetXaxis()->GetXmin();
    double xmax = h_ratio->GetXaxis()->GetXmax();
    h_ratio->GetXaxis()->SetLimits(xmin, xmax);

    return exp_count;
}
