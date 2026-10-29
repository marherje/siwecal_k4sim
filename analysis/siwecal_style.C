// SiW-ECAL plot style for ROOT: the colours of analysis/plot_style.py as TColors, and a basic example.
//   root -l -b -q siwecal_style.C        -> siwecal_style_root.png
#include "TCanvas.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStyle.h"

// --- palette: TColor::GetColor creates the colour once and returns its index ---------------------------------
namespace siwecal {
const Int_t kData = TColor::GetColor("#29a3dc");    // data: points with error bars
const Int_t kNoDigi = TColor::GetColor("#eda100");  // simulation without digitisation: dotted line
const Int_t kDigi = TColor::GetColor("#c0392b");    // simulation, digitised: hatched histogram
// several digitised models: light to dark reds, each with its own hatch and line style
const Int_t kReds[5] = {TColor::GetColor("#f08c76"), TColor::GetColor("#e0583c"), TColor::GetColor("#c0392b"),
                        TColor::GetColor("#8e2116"), TColor::GetColor("#5e1008")};
const Int_t kHatch[5] = {3004, 3005, 3013, 3006, 3007};  // ////, \\\\, xx, vertical, horizontal lines
const Int_t kLine[5] = {1, 2, 4, 3, 5};                  // solid, dashed, dash-dot, dotted, long dash-dot

void SetStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetGridColor(TColor::GetColor("#e4e3df"));
  gStyle->SetGridStyle(1);
  gStyle->SetFrameLineColor(TColor::GetColor("#52514e"));
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetHatchesLineWidth(2);
  gStyle->SetHatchesSpacing(1.2);
  gStyle->SetEndErrorSize(0);
}
void Data(TH1* h) {  // points with error bars (and the bin width as a horizontal bar)
  h->SetMarkerStyle(20); h->SetMarkerSize(1.0);
  h->SetMarkerColor(kData); h->SetLineColor(kData); h->SetLineWidth(2);
}
void NoDigi(TH1* h) { h->SetLineColor(kNoDigi); h->SetLineStyle(3); h->SetLineWidth(3); h->SetFillStyle(0); }
void Digi(TH1* h, int i = 2) {  // i-th digitised model, 0 (light) .. 4 (dark); 2 = the reference red
  h->SetLineColor(kReds[i]); h->SetLineStyle(kLine[i % 5]); h->SetLineWidth(2);
  h->SetFillColor(kReds[i]); h->SetFillStyle(kHatch[i % 5]);
}
}  // namespace siwecal

// --- basic example: an energy distribution, data against undigitised and two digitised simulations ------------
void siwecal_style() {
  siwecal::SetStyle();
  TRandom3 rng(1);
  auto make = [&](const char* name, int n, double mu, double sigma) {
    auto* h = new TH1D(name, ";Energy [MIP];Events (unit area)", 40, 1400, 2200);
    for (int i = 0; i < n; ++i) h->Fill(rng.Gaus(mu, sigma));
    h->Scale(1.0 / h->Integral("width"));
    return h;
  };
  TH1D* data = make("data", 4000, 1800, 95);
  TH1D* nodigi = make("nodigi", 40000, 1830, 80);
  TH1D* v7 = make("v7", 40000, 1840, 92);
  TH1D* v8 = make("v8", 40000, 1805, 94);
  siwecal::Data(data);
  siwecal::NoDigi(nodigi);
  siwecal::Digi(v7, 1);  // lighter red, \\\\ hatch, dashed
  siwecal::Digi(v8, 3);  // dark red, vertical lines

  auto* c = new TCanvas("c", "", 800, 600);
  c->SetLeftMargin(0.14); c->SetBottomMargin(0.12);
  v7->SetMaximum(1.6 * nodigi->GetMaximum());
  v7->Draw("hist");
  v8->Draw("hist same");
  nodigi->Draw("hist same");
  data->Draw("E1 X0 same");  // points with vertical error bars; drop X0 to show the bin width
  auto* leg = new TLegend(0.16, 0.70, 0.55, 0.88);
  leg->AddEntry(data, "Data", "pe");
  leg->AddEntry(nodigi, "Sim (no digitisation)", "l");
  leg->AddEntry(v7, "Sim (digitised), model 1", "f");
  leg->AddEntry(v8, "Sim (digitised), model 2", "f");
  leg->Draw();
  c->SaveAs("siwecal_style_root.png");
}
