#include <math.h>
#include <fstream>
using namespace std;

#include "TStyle.h"
#include "TList.h"
#include "TH3.h"

Int_t PoCA_rec_0724(TString ROOTfilename)
{
	const Int_t Nx=80,Ny=80,Nz=80;//控制像元大小
  const Double_t Xdown=-200., Xup=200., Ydown=-200., Yup=200., Zdown=-500., Zup=500.;

  Double_t vox=(Xup-Xdown)/Nx,voy=(Yup-Ydown)/Ny,voz=(Zup-Zdown)/Nz;
  Double_t sig[Nx][Ny][Nz]={0};
  Int_t count[Nx][Ny][Nz]={0};
  Int_t num=0,fixflag;
  //---读取用于重建图像的数据------------------------------
  TFile *f = TFile::Open (ROOTfilename,"read");
  TTree *Trec = (TTree *)f->Get("Trec");
  Long64_t i,j,k;
 Long64_t nentries = Trec->GetEntries();
  Int_t time[32];
  Double_t x,y,z,xt,yt,zt,dthetaX,dthetaY,mp,nom_ang,ang,S,L,Ep=1.,d;
  Int_t u,v,w,ut,vt,wt;
  Float_t x2,y2,z2,x3,y3,z3,x1,x4,y1,y4,z1,z4;
  Trec->SetBranchAddress("x",&x);
  Trec->SetBranchAddress("y",&y);
  Trec->SetBranchAddress("z",&z);
  Trec->SetBranchAddress("x1",&x1);
  Trec->SetBranchAddress("x2",&x2);
  Trec->SetBranchAddress("x3",&x3);
  Trec->SetBranchAddress("x4",&x4);
  Trec->SetBranchAddress("y1",&y1);
  Trec->SetBranchAddress("y2",&y2);
  Trec->SetBranchAddress("y3",&y3);
  Trec->SetBranchAddress("y4",&y4);
  Trec->SetBranchAddress("z1",&z1);
  Trec->SetBranchAddress("z2",&z2);
  Trec->SetBranchAddress("z3",&z3);
  Trec->SetBranchAddress("z4",&z4);
  Trec->SetBranchAddress("time",&time);  
  Trec->SetBranchAddress("ang",&ang);
  Trec->SetBranchAddress("d",&d);
  //--------------------------------------------------

  for(i=0;i<nentries;i++)
    {
      Trec->GetEntry(i);
      L =sqrt(pow((x2-x),2)+pow((y2-y),2)+pow((z2-z),2))+sqrt(pow((x3-x),2)+pow((y3-y),2)+pow((z3-z),2));
      L=voz/10;
    if (time[1]<0 || time[10]<0 || time[11]<0 || time[12]<0 || time[13]<0 || time[14]<0 || time[15]<0 || time[16]<0 || time[17]<0 || time[18]<0 || time[19]<0 ||
        time[20]<0 || time[21]<0 || time[22]<0 || time[23]<0 || time[24]<0 || time[25]<0 || time[28]<0 || time[29]<0 || time[30]<0 || time[31]<0) 
      continue;
    if (x1 > 145 || x1 < -145 || x2 > 145 || x2 < -145 || x3 > 145 || x3 < -145 || x4 > 145 || x4 < -145 ||
        y1 > 145 || y1 < -145 || y2 > 145 || y2 < -145 || y3 > 145 || y3 < -145 || y4 > 145 || y4 < -145)
     continue;
	  if(ang>0.05)	{num++;  
      S = ang*1000;
      //给PoCA点以及路径上所像元赋值--------------------------
      u=(Int_t)((x-Xdown)/vox);
      v=(Int_t)((y-Ydown)/voy);
      w=(Int_t)((z-Zdown)/voz);
      if(u<0||u>=Nx||v<0||v>=Ny||w<0||w>=Nz) continue;
      sig[u][v][w]+=S;
      count[u][v][w]++;
    //--------------------------------------------------
	}
	}

  TH3D *h3= new TH3D("d2h3","xyz",Nx,Xdown,Xup,Ny,Ydown,Yup,Nz,Zdown,Zup);
  TH2D *h2yz= new TH2D("d2h2yz","yz",Ny,Ydown,Yup,Nz,Zdown,Zup);
  TH2D *h2xy= new TH2D("d2h2xy","xy",Nx,Xdown,Xup,Ny,Ydown,Yup);
  TH2D *h2xz= new TH2D("d2h2xz","xz",Nx,Xdown,Xup,Nz,Zdown,Zup);

  //--计算每个像元的平均散射强度-----------------------------------
  for(i=0;i<Nx;i++)
    {
      for(j=0;j<Ny;j++)
	{
	  for(k=0;k<Nz;k++)
	    {
	      if(count[i][j][k]!=0)
		sig[i][j][k]=sig[i][j][k]/count[i][j][k];
	      
	      xt=i*vox+Xdown+vox/2;
	      yt=j*voy+Ydown+voy/2;
	      zt=k*voz+Zdown+voz/2;

		  h3->Fill(xt,yt,zt,count[i][j][k]);
	      h2yz->Fill(yt,zt,count[i][j][k]);
	      h2xy->Fill(xt,yt,count[i][j][k]);
	      h2xz->Fill(xt,zt,count[i][j][k]);
	    }
	}
    }
  //--------------------------------------------------
  
  //-画图-----------------------------------------------

  TCanvas *can1=new TCanvas("d2c1","abc",300,500);
  gStyle->SetOptStat(0);
  h2yz->Draw("colz");
  h2yz->SetTitle("");
  h2yz->SetXTitle("#it{Y} (mm)");
  h2yz->SetYTitle("#it{Z} (mm)");
  h2yz->GetYaxis()->SetTitleOffset(1.5);
  h2yz->GetZaxis()->SetLabelOffset(-0.001);

    // 创建TLine对象来绘制矩形框
  TLine *line1 = new TLine(-110, -110, 110, -110);  // 底边
  TLine *line2 = new TLine(110, -110, 110, 110);    // 右边
  TLine *line3 = new TLine(110, 110, -110, 110);    // 顶边
  TLine *line4 = new TLine(-110, 110, -110, -110);  // 左边

  // 设置线的颜色为红色
  line1->SetLineColor(kRed);
  line2->SetLineColor(kRed);
  line3->SetLineColor(kRed);
  line4->SetLineColor(kRed);

  line1->SetLineWidth(6);
  line2->SetLineWidth(6);
  line3->SetLineWidth(6);
  line4->SetLineWidth(6);

  line1->SetLineStyle(2);  // 设置为虚线
  line2->SetLineStyle(2);  // 设置为虚线
  line3->SetLineStyle(2);
  line4->SetLineStyle(2);

  // 在画图中绘制这些线
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  // 绘制从矩形框到坐标轴的虚线
  /*TLine *dashedLine1 = new TLine(-110, -110, -110, Zdown);  // 从矩形左下角到Y轴
  TLine *dashedLine2 = new TLine(110, -110, 110, Zdown);
  TLine *dashedLine3 = new TLine(-110, -110, Ydown, -110);
  TLine *dashedLine4 = new TLine(-110, 110, Ydown, 110);

  // 设置虚线样式
  dashedLine1->SetLineStyle(2);  // 设置为虚线
  dashedLine2->SetLineStyle(2);  // 设置为虚线
  dashedLine3->SetLineStyle(2);
  dashedLine4->SetLineStyle(2);

  // 设置虚线颜色（例如，蓝色）
  dashedLine1->SetLineColor(kRed);
  dashedLine2->SetLineColor(kRed);
  dashedLine3->SetLineColor(kRed);
  dashedLine4->SetLineColor(kRed);

  // 设置虚线宽度
  dashedLine1->SetLineWidth(6);
  dashedLine2->SetLineWidth(6);
  dashedLine3->SetLineWidth(6);
  dashedLine4->SetLineWidth(6);

  // 绘制虚线
  dashedLine1->Draw("same");
  dashedLine2->Draw("same");
  dashedLine3->Draw("same");
  dashedLine4->Draw("same");*/

  // 使用 TLatex 在虚线末端添加标签
  /*TLatex *label1 = new TLatex(-180, Zdown+5, "z = -110");
  label1->SetTextColor(kRed);  // 设置标签颜色为蓝色
  label1->SetTextSize(0.04);  // 设置标签字体大小
  label1->Draw();

  TLatex *label2 = new TLatex(40, Zdown+5, "z = 110");
  label2->SetTextColor(kRed);  // 设置标签颜色为蓝色
  label2->SetTextSize(0.04);  // 设置标签字体大小
  label2->Draw();

  TLatex *label3 = new TLatex(Ydown+5, -150, "y = -110");
  label3->SetTextColor(kRed);  // 设置标签颜色为蓝色
  label3->SetTextSize(0.04);  // 设置标签字体大小
  label3->Draw();

  TLatex *label4 = new TLatex(Ydown+5, 70, "y = 110");
  label4->SetTextColor(kRed);  // 设置标签颜色为蓝色
  label4->SetTextSize(0.04);  // 设置标签字体大小
  label4->Draw();*/

  cout<<"num="<<num<<endl;
  return nentries;
}
