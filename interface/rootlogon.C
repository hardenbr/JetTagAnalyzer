#include "TStyle.h"

class defaultStyle{

public:
  void setStyle() {
    Int_t font = 42;
    gStyle->SetLegendFillColor(kWhite);

    gStyle->SetOptTitle(00000);
    gStyle->SetTitleX(0.3);
    gStyle->SetTitleW(0.4);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasDefH(600); //Height of canvas
    gStyle->SetCanvasDefW(800); //Width of canvas
    gStyle->SetCanvasDefX(0);   //POsition on screen
    gStyle->SetCanvasDefY(0);

    gStyle->SetPadBorderMode(0);
    // gStyle->SetPadBorderSize(Width_t size = 1);
    gStyle->SetPadColor(kWhite);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetGridColor(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    //for the frame
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetFrameLineWidth(1);

    gStyle->SetPaperSize(20,26);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetPadLeftMargin(0.15);

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFont(font,"xyz");  // set the all 3 axes title font
    gStyle->SetTitleFont(font," ");    // set the pad title font
    gStyle->SetTitleSize(0.07,"xyz"); // set the 3 axes title size
    gStyle->SetTitleSize(0.07," ");   // set the pad title size
    gStyle->SetLabelFont(font,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelColor(1,"xyz");
    gStyle->SetTextFont(font);
    gStyle->SetTextSize(0.08);
    gStyle->SetStatFont(font);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleW(0.4);

    //tick marks
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

  }
  //  gSystem->Load("libFWCoreFWLite.so");
  //AutoLibraryLoader::enable();
};
