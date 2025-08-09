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

Int_t air_hist0714(TString expROOTfilename="250212_63d_air.root",
    TString simROOTfilename="smeared_rec_CryMuAna_0_749.root")
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

    TCanvas *canvas = new TCanvas("canvas", "ang", 800, 600);

    TH1F *h_air = new TH1F("h_air", "distribution of scattering angle #theta", 50, 0.05, 0.5);
    TH1F *h_obs = new TH1F("h_obs", "distribution of scattering angle #theta(air)", 50, 0.05, 0.5);
    TH1F *h_mu = new TH1F("h_mu", "distribution of scattering angle #theta", 50, 0.05, 0.5);
    TH1F *h_e = new TH1F("h_e", "distribution of scattering angle #theta", 50, 0.05, 0.5);
    TH1F *h_other = new TH1F("h_other", "distribution of scattering angle #theta", 50, 0.05, 0.5);

    int exp_count=0;
    int sim_count=0;
    int mu_count=0;
    int e_count=0;
    int other_count=0;

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

    cout<<"exp_count="<<exp_count<<endl;
    cout<<"sim_count="<<sim_count<<endl;
    cout<<"mu="<<mu_count<<endl;
    cout<<"e="<<e_count<<endl;
    cout<<"other="<<other_count<<endl;

    double sigma_mu = (100 * 1.0) * pow((((mu_count * 1.0) * ((mu_count * 1.0) + (sim_count * 1.0))) / (pow((sim_count * 1.0) , 3))) , 0.5);
    double sigma_e = (100 * 1.0) * pow((((e_count * 1.0) * ((e_count * 1.0) + (sim_count * 1.0))) / (pow((sim_count * 1.0) , 3))) , 0.5);
    double sigma_other = (100 * 1.0) * pow((((other_count * 1.0) * ((other_count * 1.0) + (sim_count * 1.0))) / (pow((sim_count * 1.0) , 3))) , 0.5);
    cout<<"sigma_mu="<<sigma_mu<<endl;
    cout<<"sigma_e="<<sigma_e<<endl;
    cout<<"sigma_other"<<sigma_other<<endl;

    //归一化
    for (int bin = 1; bin <= 50; bin++){
        if (h_mu->GetBinContent(bin) == 0 ) h_mu->SetBinContent(bin,(h_mu->GetBinContent(bin - 1) + h_mu->GetBinContent(bin + 1)) * 1.0 / 2);
        if (h_e->GetBinContent(bin) == 0 ) h_e->SetBinContent(bin,(h_e->GetBinContent(bin - 1) + h_e->GetBinContent(bin + 1)) * 1.0 / 2);
        if (h_other->GetBinContent(bin) == 0 ) h_other->SetBinContent(bin,(h_other->GetBinContent(bin - 1) + h_other->GetBinContent(bin + 1)) * 1.0 / 2);
    }

    //h_mu->Scale(1222680./5721816);
    //h_e->Scale(1222680./5721816);
    //h_other->Scale(1222680./5721816);
    //h_obs->Scale(1179143./148832);

    /*h_obs->Scale((1.0/h_obs->Integral()));
    h_mu->Scale(((1.0)/sim_count));
    h_e->Scale(((1.0)/sim_count));
    h_other->Scale(((1.0)/sim_count));*/

    h_obs->Scale(((1.0)/h_obs->Integral()));
    h_mu->Scale(((1.0)/h_mu->Integral())*0.351485);
    h_e->Scale(((1.0)/h_e->Integral())*0.525415);
    h_other->Scale(((1.0)/h_other->Integral())*0.1231);

    h_air->Add(h_mu);
    h_air->Add(h_e);
    h_air->Add(h_other);

    /*TFile *file = new TFile("hist_of_air.root", "RECREATE");
    
    h_obs->Write();
    h_mu->Write();
    h_e->Write();
    h_other->Write();

    file->Close();*/

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
    h_obs->SetMinimum(5e-5);
    h_obs->SetMaximum(1);
	h_obs->GetYaxis()->SetTitle("relative events");
    h_obs->GetXaxis()->SetLabelSize(0);     // 隐藏X轴标签
    h_obs->GetXaxis()->SetTitleSize(0);     // 隐藏X轴标题
    h_obs->GetYaxis()->SetLabelSize(0.05);     // Y轴标签
    h_obs->GetYaxis()->SetTitleSize(0.05);     // Y轴标题
    h_obs->SetMarkerStyle(20);
    h_obs->SetMarkerSize(1.1);
	h_obs->Draw("E");
	h_obs->SetLineColor(kBlack); //黑色
    h_obs->SetLineWidth(3);

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

    TLegend *legend = new TLegend(0.2, 0.7, 0.9, 0.9);
    legend->SetNColumns(2); 
	legend->SetX1NDC(0.6);
    legend->SetY1NDC(0.70);
    legend->SetX2NDC(0.9);
    legend->SetY2NDC(0.9);
    legend->SetTextSize(0.05);
    legend->AddEntry(h_obs, "observed", "lep");
    legend->AddEntry(h_air, "MC_all_particle", "l");
    legend->AddEntry(h_mu, Form("MC_mu (%.2f%% #pm 0.49%%)", (h_mu->Integral() * 1.0)/(h_air->Integral() * 1.0)*100.0 ), "l");
    legend->AddEntry(h_e, Form("MC_e (%.2f%% #pm 0.49%%)", (h_e->Integral() * 1.0)/(h_air->Integral() * 1.0)*100.0), "l");
    legend->AddEntry(h_other, Form("MC_other (%.2f%% #pm 0.19%%)", (h_other->Integral() * 1.0)/(h_air->Integral() * 1.0)*100.0), "l");
    legend->Draw();

    pad_bot->cd();
    gStyle->SetOptStat(0);
    // 检查Bin的一致性
    h_obs->Rebin(5);
    h_air->Rebin(5);
    if (h_obs->GetNbinsX() != h_air->GetNbinsX()) {
        std::cerr << "Error: Binning of h_obs and h_air do not match!" << std::endl;
        return 0;
    }

    // 创建一个新的直方图用于存储比值
    TH1F* h_ratio = (TH1F*)h_obs->Clone("h_ratio");
    h_ratio->SetTitle("");
    //h_ratio->GetYaxis()->SetRangeUser(-1, 10); // 设置合理的 Y 轴范围
    h_ratio->SetMaximum(1.46);
    h_ratio->SetMinimum(0.46);
    h_ratio->GetXaxis()->SetTitle("#theta (rad)");
    h_ratio->GetYaxis()->SetLabelSize(0.08);     // Y轴标签
    h_ratio->GetYaxis()->SetTitleSize(0.08);     // Y轴标题
    h_ratio->GetYaxis()->SetTitle("observed/MC");
    h_ratio->GetYaxis()->SetTitleOffset(0.7);
    h_ratio->GetXaxis()->SetLabelSize(0.08);     // X轴标签
    h_ratio->GetXaxis()->SetTitleSize(0.08);     // X轴标题
    h_ratio->SetStats(0);

    // 计算比值并处理误差
    h_ratio->Divide(h_air);
    TF1* fit_ratio = new TF1("fit_ratio", "pol0", h_ratio->GetXaxis()->GetXmin(), h_ratio->GetXaxis()->GetXmax());

    // 对 h_ratio 进行拟合
    h_ratio->Fit(fit_ratio, "R");

    // 设置绘图样式
    h_ratio->SetMarkerStyle(20);   // 数据点样式
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->SetMarkerSize(1.5);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetLineWidth(3);
    fit_ratio->SetLineColor(kRed); // 拟合线颜色
    fit_ratio->SetLineWidth(3);
    //fit_ratio->SetLineStyle(2);    // 虚线样式

    // 获取拟合常数值和误差
    double constant = fit_ratio->GetParameter(0);  // 获取常数参数值
    double constant_error = fit_ratio->GetParError(0);  // 获取常数参数的误差

    // 绘制图像
    h_ratio->Draw("E"); // 绘制数据点，带误差条
    fit_ratio->Draw("SAME"); // 在同一张图上绘制拟合曲线

    // 绘制拟合常数及误差
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.08);
    latex->SetTextAlign(13); // 左对齐
    latex->DrawLatex(0.67, 0.87, Form("Constant = %.3f #pm %.3f", constant, constant_error));

    // 添加图例
    TLegend* legend1 = new TLegend(0.5, 0.9, 0.9, 1.0);
    legend1->SetNColumns(2); 
    legend1->SetTextSize(0.08);
    legend1->AddEntry(h_ratio, "observed/MC Ratio", "lep");
    legend1->AddEntry(fit_ratio, "Fit (constant)", "l");
    legend1->Draw();

    double xmin = h_obs->GetXaxis()->GetXmin();
    double xmax = h_ratio->GetXaxis()->GetXmax();
    h_ratio->GetXaxis()->SetLimits(xmin, xmax);

    return exp_count;
}