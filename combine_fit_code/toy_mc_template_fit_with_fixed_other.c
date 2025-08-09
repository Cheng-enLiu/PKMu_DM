#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include <vector>
#include <iostream>

using namespace RooFit;

void toy_mc_template_fit_with_fixed_other(int nToys = 100) {
    // 1. 读文件和直方图
    TFile* file = TFile::Open("hist_of_air.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TH1D* h_obs   = (TH1D*)file->Get("h_obs");
    TH1D* h_mu_0  = (TH1D*)file->Get("h_mu");
    TH1D* h_e_0   = (TH1D*)file->Get("h_e");
    TH1D* h_other_0 = (TH1D*)file->Get("h_other");

    if (!h_obs || !h_mu_0 || !h_e_0 || !h_other_0) {
        std::cerr << "Histograms not found!" << std::endl;
        return;
    }

    h_obs->Sumw2();
    h_mu_0->Sumw2();
    h_e_0->Sumw2();
    h_other_0->Sumw2();

    // 2. 定义 RooFit 变量
    RooRealVar x("x", "Observable", h_obs->GetXaxis()->GetXmin(), h_obs->GetXaxis()->GetXmax());

    RooDataHist dh_obs("dh_obs", "Observed Data", x, Import(*h_obs));

    // 3. 定义固定权重
    const double w_other_val = 0.1191;

    // 4. 初始化随机数发生器，固定种子保证结果稳定复现
    TRandom3 rand(2025);

    // 5. 记录拟合权重结果
    std::vector<double> w_mu_vals, w_e_vals;

    for (int i = 0; i < nToys; ++i) {
        // --- 对 h_mu 加扰动 ---
        TH1D* h_mu = (TH1D*)h_mu_0->Clone(Form("h_mu_toy_%d", i));
        h_mu->Reset();
        for (int b = 1; b <= h_mu_0->GetNbinsX(); ++b) {
            double val = h_mu_0->GetBinContent(b);
            double fluct = rand.Poisson(val);
            if (val != 0&& fluct == 0) fluct = std::max(1e-3, 1e-6 * val);
            h_mu->SetBinContent(b, fluct);
            h_mu->SetBinError(b, sqrt(fluct));
        }

        // --- 对 h_e 加扰动 ---
        TH1D* h_e = (TH1D*)h_e_0->Clone(Form("h_e_toy_%d", i));
        h_e->Reset();
        for (int b = 1; b <= h_e_0->GetNbinsX(); ++b) {
            double val = h_e_0->GetBinContent(b);
            double fluct = rand.Poisson(val);
            if (val != 0 && fluct == 0) fluct = std::max(1e-3, 1e-6 * val);
            h_e->SetBinContent(b, fluct);
            h_e->SetBinError(b, sqrt(fluct));
        }

        // --- 对 h_other 加扰动 ---
        TH1D* h_other = (TH1D*)h_other_0->Clone(Form("h_other_toy_%d", i));
        h_other->Reset();
        for (int b = 1; b <= h_other_0->GetNbinsX(); ++b) {
            double val = h_other_0->GetBinContent(b);
            double fluct = rand.Poisson(val);
            if (val != 0 && fluct == 0) fluct = std::max(1e-3, 1e-6 * val);
            h_other->SetBinContent(b, fluct);
            h_other->SetBinError(b, sqrt(fluct));
        }

        // 6. 构建对应 RooDataHist 和 RooHistPdf
        RooDataHist dh_mu(Form("dh_mu_%d", i), "Muon Template", x, Import(*h_mu));
        RooDataHist dh_e(Form("dh_e_%d", i), "Electron Template", x, Import(*h_e));
        RooDataHist dh_other(Form("dh_other_%d", i), "Other Template", x, Import(*h_other));

        RooHistPdf pdf_mu(Form("pdf_mu_%d", i), "PDF from mu", x, dh_mu);
        RooHistPdf pdf_e(Form("pdf_e_%d", i), "PDF from e", x, dh_e);
        RooHistPdf pdf_other(Form("pdf_other_%d", i), "PDF from other", x, dh_other);

        // 7. 定义可浮动参数 w_mu 和 w_e，权重总和=1，w_other固定
        RooRealVar w_mu("w_mu", "Weight mu", 0.3663, 0.0, 1.0);
        RooRealVar w_other("w_other", "Weight other", w_other_val);
        w_other.setConstant(kTRUE);
        RooFormulaVar w_e("w_e", "Weight e", "1.0 - w_mu - w_other", RooArgList(w_mu, w_other));

        // 8. 构建模型
        RooAddPdf model("model", "Combined Model",
                        RooArgList(pdf_mu, pdf_e, pdf_other),
                        RooArgList(w_mu, w_e));

        // 9. 拟合
        RooFitResult* fitResult = model.fitTo(dh_obs,
                                              Extended(kFALSE),
                                              PrintLevel(-1),
                                              Save());

        if (fitResult->status() != 0) {
            std::cerr << "Fit failed in toy " << i << std::endl;
            std::cerr << "Integral(h_mu) = " << h_mu->Integral() << ", Integral(h_e) = " << h_e->Integral() << ", Integral(h_other) = " << h_other->Integral() << std::endl;
            delete h_mu; delete h_e; delete h_other; delete fitResult;
            continue;
        }

        // 10. 记录拟合结果
        w_mu_vals.push_back(w_mu.getVal());
        w_e_vals.push_back(w_e.getVal());

        delete h_mu; delete h_e; delete h_other; delete fitResult;
    }

    // 11. 计算均值和标准差函数
    auto mean_std = [](const std::vector<double>& vec) {
        double sum = 0, sum2 = 0;
        for (double v : vec) {
            sum += v;
            sum2 += v * v;
        }
        double mean = sum / vec.size();
        double std = sqrt(sum2 / vec.size() - mean * mean);
        return std::make_pair(mean, std);
    };

    auto [mu_mean, mu_std] = mean_std(w_mu_vals);
    auto [e_mean,  e_std]  = mean_std(w_e_vals);

    // 12. 输出结果
    std::cout << "\n=== Toy MC Fit Summary over " << w_mu_vals.size() << " toys ===" << std::endl;
    std::cout << "w_mu: mean = " << mu_mean << ", std = " << mu_std << std::endl;
    std::cout << "w_e : mean = " << e_mean  << ", std = " << e_std  << std::endl;
    std::cout << "w_other: fixed = " << w_other_val << std::endl;
    std::cout << "Total weight sum = " << (mu_mean + e_mean + w_other_val) << std::endl;
}
