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

void toy_mc_template_fit_with_floating_other(int nToys = 100) {
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

    RooRealVar x("x", "Observable", h_obs->GetXaxis()->GetXmin(), h_obs->GetXaxis()->GetXmax());
    RooDataHist dh_obs("dh_obs", "Observed Data", x, Import(*h_obs));

    TRandom3 rand(2025);
    std::vector<double> w_mu_vals, w_e_vals, w_other_vals;

    for (int i = 0; i < nToys; ++i) {
        TH1D* h_mu = (TH1D*)h_mu_0->Clone(Form("h_mu_toy_%d", i));
        h_mu->Reset();
        for (int b = 1; b <= h_mu_0->GetNbinsX(); ++b) {
            double val = h_mu_0->GetBinContent(b);
            double fluct = rand.Poisson(val);
            if (val != 0&& fluct == 0) fluct = std::max(1e-3, 1e-6 * val);
            h_mu->SetBinContent(b, fluct);
            h_mu->SetBinError(b, sqrt(fluct));
        }

        TH1D* h_e = (TH1D*)h_e_0->Clone(Form("h_e_toy_%d", i));
        h_e->Reset();
        for (int b = 1; b <= h_e_0->GetNbinsX(); ++b) {
            double val = h_e_0->GetBinContent(b);
            double fluct = rand.Poisson(val);
            if (val != 0&& fluct == 0) fluct = std::max(1e-3, 1e-6 * val);
            h_e->SetBinContent(b, fluct);
            h_e->SetBinError(b, sqrt(fluct));
        }

        TH1D* h_other = (TH1D*)h_other_0->Clone(Form("h_other_toy_%d", i));
        h_other->Reset();
        for (int b = 1; b <= h_other_0->GetNbinsX(); ++b) {
            double val = h_other_0->GetBinContent(b);
            double fluct = rand.Poisson(val);
            if (val != 0&& fluct == 0) fluct = std::max(1e-3, 1e-6 * val);
            h_other->SetBinContent(b, fluct);
            h_other->SetBinError(b, sqrt(fluct));
        }

        RooDataHist dh_mu(Form("dh_mu_%d", i), "Muon Template", x, Import(*h_mu));
        RooDataHist dh_e(Form("dh_e_%d", i), "Electron Template", x, Import(*h_e));
        RooDataHist dh_other(Form("dh_other_%d", i), "Other Template", x, Import(*h_other));

        RooHistPdf pdf_mu(Form("pdf_mu_%d", i), "PDF from mu", x, dh_mu);
        RooHistPdf pdf_e(Form("pdf_e_%d", i), "PDF from e", x, dh_e);
        RooHistPdf pdf_other(Form("pdf_other_%d", i), "PDF from other", x, dh_other);

        RooRealVar w_mu("w_mu", "Weight mu", 0.3526, 0.0, 1.0);
        RooRealVar w_other("w_other", "Weight other", 0.1327, 0.0, 1.0);
        RooFormulaVar w_e("w_e", "Weight e", "1.0 - w_mu - w_other", RooArgList(w_mu, w_other));

        RooAddPdf model("model", "Combined Model",
                        RooArgList(pdf_mu, pdf_e, pdf_other),
                        RooArgList(w_mu, w_e, w_other));

        RooFitResult* fitResult = model.fitTo(dh_obs,
                                              Extended(kFALSE),
                                              PrintLevel(-1),
                                              Save());

        if (fitResult->status() != 0) {
            std::cerr << "Fit failed in toy " << i << std::endl;
            delete h_mu; delete h_e; delete h_other; delete fitResult;
            continue;
        }

        w_mu_vals.push_back(w_mu.getVal());
        w_other_vals.push_back(w_other.getVal());
        w_e_vals.push_back(w_e.getVal());

        delete h_mu; delete h_e; delete h_other; delete fitResult;
    }

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
    auto [other_mean, other_std] = mean_std(w_other_vals);

    std::cout << "\n=== Toy MC Fit Summary over " << w_mu_vals.size() << " toys ===" << std::endl;
    std::cout << "w_mu   : mean = " << mu_mean    << ", std = " << mu_std << std::endl;
    std::cout << "w_e    : mean = " << e_mean     << ", std = " << e_std  << std::endl;
    std::cout << "w_other: mean = " << other_mean << ", std = " << other_std << std::endl;
    std::cout << "Total weight sum = " << (mu_mean + e_mean + other_mean) << std::endl;
}
