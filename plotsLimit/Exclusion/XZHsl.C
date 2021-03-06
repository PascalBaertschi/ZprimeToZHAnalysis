void XZHsl()
{
//=========Macro generated from canvas: c1/Exclusion Limits
//=========  (Sun Feb 16 18:00:17 2020) by ROOT version6.10/09
   TCanvas *c1 = new TCanvas("c1", "Exclusion Limits",0,0,800,600);
   gStyle->SetOptStat(0);
   c1->SetHighLightColor(2);
   c1->Range(171.1133,-2.678449,5271.836,4.106039);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.12);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.06);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   Double_t _fx3001[41] = {
   800,
   900,
   1000,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy3001[41] = {
   25.9375,
   13.9375,
   15.1875,
   4.5938,
   3.5312,
   2.9297,
   2.3906,
   2.0078,
   1.6797,
   1.4258,
   1.2109,
   1.0312,
   0.8945,
   0.7812,
   0.6914,
   0.6172,
   0.5527,
   0.5,
   0.4551,
   0.418,
   0.3857,
   0.3574,
   0.335,
   0.3164,
   0.2998,
   0.2871,
   0.2754,
   0.2656,
   0.2578,
   0.251,
   0.2441,
   0.2373,
   0.2314,
   0.2266,
   0.2217,
   0.2178,
   0.2129,
   0.21,
   0.208,
   0.207,
   0.207};
   Double_t _felx3001[41] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3001[41] = {
   5.0659,
   5.4988,
   2.907,
   2.3149,
   1.8345,
   1.5564,
   1.2793,
   1.0823,
   0.9055,
   0.7686,
   0.6527,
   0.5559,
   0.4857,
   0.4272,
   0.3808,
   0.3448,
   0.3109,
   0.2852,
   0.2631,
   0.2466,
   0.2305,
   0.2178,
   0.2068,
   0.1977,
   0.1897,
   0.1839,
   0.1786,
   0.1743,
   0.1702,
   0.1667,
   0.164,
   0.1594,
   0.1573,
   0.154,
   0.1516,
   0.1497,
   0.1464,
   0.1452,
   0.1438,
   0.1439,
   0.1439};
   Double_t _fehx3001[41] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3001[41] = {
   15.0361,
   11.9863,
   8.1927,
   5.285,
   4.2509,
   3.6563,
   3.0357,
   2.6088,
   2.2053,
   1.8692,
   1.6134,
   1.4017,
   1.2397,
   1.1238,
   1.0214,
   0.9362,
   0.8661,
   0.8133,
   0.7679,
   0.7254,
   0.695,
   0.667,
   0.6464,
   0.6238,
   0.6089,
   0.5832,
   0.5594,
   0.5396,
   0.5238,
   0.5099,
   0.4961,
   0.4822,
   0.4703,
   0.4603,
   0.4504,
   0.4425,
   0.4326,
   0.4266,
   0.4227,
   0.4207,
   0.4207};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(41,_fx3001,_fy3001,_felx3001,_fehx3001,_fely3001,_fehy3001);
   grae->SetName("");
   grae->SetTitle("");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ffcc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   
   TH1F *Graph_Graph3001 = new TH1F("Graph_Graph3001","",100,380,5420);
   Graph_Graph3001->SetMinimum(0.01);
   Graph_Graph3001->SetMaximum(5000);
   Graph_Graph3001->SetDirectory(0);
   Graph_Graph3001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3001->SetLineColor(ci);
   Graph_Graph3001->GetXaxis()->SetTitle("m_{Z'} (GeV)");
   Graph_Graph3001->GetXaxis()->SetRange(9,92);
   Graph_Graph3001->GetXaxis()->SetMoreLogLabels();
   Graph_Graph3001->GetXaxis()->SetNoExponent();
   Graph_Graph3001->GetXaxis()->SetLabelFont(42);
   Graph_Graph3001->GetXaxis()->SetLabelSize(0.045);
   Graph_Graph3001->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph3001->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph3001->GetXaxis()->SetTitleFont(42);
   Graph_Graph3001->GetYaxis()->SetTitle("#sigma(Z') #bf{#it{#Beta}}(Z' #rightarrow Zh) (fb)");
   Graph_Graph3001->GetYaxis()->SetMoreLogLabels();
   Graph_Graph3001->GetYaxis()->SetNoExponent();
   Graph_Graph3001->GetYaxis()->SetLabelFont(42);
   Graph_Graph3001->GetYaxis()->SetLabelSize(0.045);
   Graph_Graph3001->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph3001->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph3001->GetYaxis()->SetTitleFont(42);
   Graph_Graph3001->GetZaxis()->SetLabelFont(42);
   Graph_Graph3001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3001);
   
   grae->Draw("a3");
   
   Double_t _fx3002[41] = {
   800,
   900,
   1000,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy3002[41] = {
   25.9375,
   13.9375,
   15.1875,
   4.5938,
   3.5312,
   2.9297,
   2.3906,
   2.0078,
   1.6797,
   1.4258,
   1.2109,
   1.0312,
   0.8945,
   0.7812,
   0.6914,
   0.6172,
   0.5527,
   0.5,
   0.4551,
   0.418,
   0.3857,
   0.3574,
   0.335,
   0.3164,
   0.2998,
   0.2871,
   0.2754,
   0.2656,
   0.2578,
   0.251,
   0.2441,
   0.2373,
   0.2314,
   0.2266,
   0.2217,
   0.2178,
   0.2129,
   0.21,
   0.208,
   0.207,
   0.207};
   Double_t _felx3002[41] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3002[41] = {
   3.4037,
   3.4367,
   1.9531,
   1.4468,
   1.1466,
   0.9728,
   0.8096,
   0.6764,
   0.5659,
   0.4864,
   0.413,
   0.3517,
   0.3111,
   0.2737,
   0.244,
   0.2209,
   0.2016,
   0.1849,
   0.1727,
   0.1618,
   0.1513,
   0.1429,
   0.1373,
   0.1313,
   0.126,
   0.1236,
   0.12,
   0.1171,
   0.1143,
   0.112,
   0.1102,
   0.1077,
   0.1056,
   0.1047,
   0.103,
   0.1018,
   0.0995,
   0.0981,
   0.0977,
   0.0978,
   0.0978};
   Double_t _fehx3002[41] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3002[41] = {
   6.5134,
   5.1667,
   3.3296,
   2.2522,
   1.8017,
   1.5415,
   1.2769,
   1.0965,
   0.9239,
   0.7786,
   0.671,
   0.5796,
   0.5099,
   0.4578,
   0.4134,
   0.3788,
   0.3481,
   0.3229,
   0.3047,
   0.2899,
   0.2737,
   0.2593,
   0.251,
   0.2447,
   0.2366,
   0.2312,
   0.2261,
   0.2224,
   0.2199,
   0.2141,
   0.2122,
   0.21,
   0.2049,
   0.2023,
   0.1997,
   0.1962,
   0.1952,
   0.1925,
   0.1907,
   0.1898,
   0.1915};
   grae = new TGraphAsymmErrors(41,_fx3002,_fy3002,_felx3002,_fehx3002,_fely3002,_fehy3002);
   grae->SetName("");
   grae->SetTitle("");

   ci = TColor::GetColor("#00cc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#00cc00");
   grae->SetLineColor(ci);
   
   TH1F *Graph_Graph3002 = new TH1F("Graph_Graph3002","",100,380,5420);
   Graph_Graph3002->SetMinimum(0.09828);
   Graph_Graph3002->SetMaximum(35.68507);
   Graph_Graph3002->SetDirectory(0);
   Graph_Graph3002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3002->SetLineColor(ci);
   Graph_Graph3002->GetXaxis()->SetLabelFont(42);
   Graph_Graph3002->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3002->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3002->GetXaxis()->SetTitleFont(42);
   Graph_Graph3002->GetYaxis()->SetLabelFont(42);
   Graph_Graph3002->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3002->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3002->GetYaxis()->SetTitleOffset(0);
   Graph_Graph3002->GetYaxis()->SetTitleFont(42);
   Graph_Graph3002->GetZaxis()->SetLabelFont(42);
   Graph_Graph3002->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3002->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3002);
   
   grae->Draw(", 3");
   
   Double_t _fx3003[43] = {
   800,
   900,
   1000,
   1100,
   1200,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy3003[43] = {
   578.0999,
   358.3873,
   231.6927,
   154.8033,
   106.2204,
   74.50991,
   53.23298,
   38.62051,
   28.39542,
   21.11701,
   15.85439,
   12.00795,
   9.163542,
   7.039392,
   5.439621,
   4.225485,
   3.297668,
   2.584068,
   2.032477,
   1.604015,
   1.269687,
   1.009061,
   0.8013809,
   0.6388036,
   0.5101927,
   0.4081523,
   0.3269965,
   0.2623107,
   0.2106453,
   0.1693069,
   0.1361855,
   0.109623,
   0.08828075,
   0.07113544,
   0.05733185,
   0.04620448,
   0.03724663,
   0.03002855,
   0.02420639,
   0.01950782,
   0.01571446,
   0.01265338,
   0.01018433};
   Double_t _felx3003[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3003[43] = {
   39.09341,
   26.12542,
   18.07946,
   12.62288,
   9.051203,
   6.765647,
   5.163873,
   3.902392,
   2.984897,
   2.32527,
   1.835073,
   1.440604,
   1.138088,
   0.9148231,
   0.738287,
   0.5978874,
   0.4836005,
   0.3938987,
   0.3264199,
   0.2685276,
   0.2230208,
   0.1841806,
   0.1522582,
   0.1284647,
   0.1085944,
   0.09113777,
   0.07673766,
   0.06445908,
   0.05474387,
   0.04656555,
   0.0394106,
   0.03330537,
   0.02813664,
   0.02420247,
   0.02077107,
   0.01776471,
   0.01515063,
   0.01287661,
   0.01038,
   0.008365192,
   0.006738552,
   0.005425925,
   0.004367165};
   Double_t _fehx3003[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3003[43] = {
   39.09341,
   26.12542,
   18.07946,
   12.62288,
   9.051203,
   6.765647,
   5.163873,
   3.902392,
   2.984897,
   2.32527,
   1.835073,
   1.440604,
   1.138088,
   0.9148231,
   0.738287,
   0.5978874,
   0.4836005,
   0.3938987,
   0.3264199,
   0.2685276,
   0.2230208,
   0.1841806,
   0.1522582,
   0.1284647,
   0.1085944,
   0.09113777,
   0.07673766,
   0.06445908,
   0.05474387,
   0.04656555,
   0.0394106,
   0.03330537,
   0.02813664,
   0.02420247,
   0.02077107,
   0.01776471,
   0.01515063,
   0.01287661,
   0.01038,
   0.008365192,
   0.006738552,
   0.005425925,
   0.004367165};
   grae = new TGraphAsymmErrors(43,_fx3003,_fy3003,_felx3003,_fehx3003,_fely3003,_fehy3003);
   grae->SetName("");
   grae->SetTitle("");

   ci = TColor::GetColor("#ff66ff");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3013);

   ci = TColor::GetColor("#cc33cc");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   
   TH1F *Graph_Graph3003 = new TH1F("Graph_Graph3003","",100,380,5420);
   Graph_Graph3003->SetMinimum(0.005235446);
   Graph_Graph3003->SetMaximum(678.9121);
   Graph_Graph3003->SetDirectory(0);
   Graph_Graph3003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3003->SetLineColor(ci);
   Graph_Graph3003->GetXaxis()->SetLabelFont(42);
   Graph_Graph3003->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3003->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3003->GetXaxis()->SetTitleFont(42);
   Graph_Graph3003->GetYaxis()->SetLabelFont(42);
   Graph_Graph3003->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3003->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3003->GetYaxis()->SetTitleOffset(0);
   Graph_Graph3003->GetYaxis()->SetTitleFont(42);
   Graph_Graph3003->GetZaxis()->SetLabelFont(42);
   Graph_Graph3003->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3003->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3003);
   
   grae->Draw(", l3");
   
   Double_t _fx3004[43] = {
   800,
   900,
   1000,
   1100,
   1200,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy3004[43] = {
   578.0999,
   358.3873,
   231.6927,
   154.8033,
   106.2204,
   74.50991,
   53.23298,
   38.62051,
   28.39542,
   21.11701,
   15.85439,
   12.00795,
   9.163542,
   7.039392,
   5.439621,
   4.225485,
   3.297668,
   2.584068,
   2.032477,
   1.604015,
   1.269687,
   1.009061,
   0.8013809,
   0.6388036,
   0.5101927,
   0.4081523,
   0.3269965,
   0.2623107,
   0.2106453,
   0.1693069,
   0.1361855,
   0.109623,
   0.08828075,
   0.07113544,
   0.05733185,
   0.04620448,
   0.03724663,
   0.03002855,
   0.02420639,
   0.01950782,
   0.01571446,
   0.01265338,
   0.01018433};
   Double_t _felx3004[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3004[43] = {
   39.09341,
   26.12542,
   18.07946,
   12.62288,
   9.051203,
   6.765647,
   5.163873,
   3.902392,
   2.984897,
   2.32527,
   1.835073,
   1.440604,
   1.138088,
   0.9148231,
   0.738287,
   0.5978874,
   0.4836005,
   0.3938987,
   0.3264199,
   0.2685276,
   0.2230208,
   0.1841806,
   0.1522582,
   0.1284647,
   0.1085944,
   0.09113777,
   0.07673766,
   0.06445908,
   0.05474387,
   0.04656555,
   0.0394106,
   0.03330537,
   0.02813664,
   0.02420247,
   0.02077107,
   0.01776471,
   0.01515063,
   0.01287661,
   0.01038,
   0.008365192,
   0.006738552,
   0.005425925,
   0.004367165};
   Double_t _fehx3004[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3004[43] = {
   39.09341,
   26.12542,
   18.07946,
   12.62288,
   9.051203,
   6.765647,
   5.163873,
   3.902392,
   2.984897,
   2.32527,
   1.835073,
   1.440604,
   1.138088,
   0.9148231,
   0.738287,
   0.5978874,
   0.4836005,
   0.3938987,
   0.3264199,
   0.2685276,
   0.2230208,
   0.1841806,
   0.1522582,
   0.1284647,
   0.1085944,
   0.09113777,
   0.07673766,
   0.06445908,
   0.05474387,
   0.04656555,
   0.0394106,
   0.03330537,
   0.02813664,
   0.02420247,
   0.02077107,
   0.01776471,
   0.01515063,
   0.01287661,
   0.01038,
   0.008365192,
   0.006738552,
   0.005425925,
   0.004367165};
   grae = new TGraphAsymmErrors(43,_fx3004,_fy3004,_felx3004,_fehx3004,_fely3004,_fehy3004);
   grae->SetName("");
   grae->SetTitle("");

   ci = TColor::GetColor("#ff66ff");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3013);

   ci = TColor::GetColor("#cc33cc");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   
   TH1F *Graph_Graph_Graph30033004 = new TH1F("Graph_Graph_Graph30033004","",100,380,5420);
   Graph_Graph_Graph30033004->SetMinimum(0.005235446);
   Graph_Graph_Graph30033004->SetMaximum(678.9121);
   Graph_Graph_Graph30033004->SetDirectory(0);
   Graph_Graph_Graph30033004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph30033004->SetLineColor(ci);
   Graph_Graph_Graph30033004->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph30033004->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph30033004->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph30033004->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph30033004->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph30033004->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph30033004->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph30033004->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_Graph30033004->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph30033004->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph30033004->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph30033004->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph30033004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph30033004);
   
   grae->Draw(", l3x0y0");
   
   Double_t _fx3005[43] = {
   800,
   900,
   1000,
   1100,
   1200,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy3005[43] = {
   485.1621,
   367.2001,
   263.8817,
   188.6818,
   135.7415,
   98.58745,
   72.32583,
   53.57184,
   40.04648,
   30.18608,
   22.91748,
   17.51997,
   13.47566,
   10.42173,
   8.100062,
   6.323773,
   4.956944,
   3.899242,
   3.077397,
   2.435993,
   1.93346,
   1.540306,
   1.225945,
   0.9791569,
   0.7834139,
   0.627741,
   0.5036664,
   0.4045769,
   0.3252906,
   0.2617513,
   0.2107641,
   0.1698193,
   0.1368799,
   0.1103879,
   0.08903561,
   0.07180615,
   0.0579232,
   0.0467275,
   0.03768939,
   0.03039015,
   0.02449322,
   0.01973149,
   0.01588843};
   Double_t _felx3005[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3005[43] = {
   32.80858,
   26.76786,
   20.59123,
   15.38538,
   11.56674,
   8.951936,
   7.015978,
   5.413142,
   4.209645,
   3.323897,
   2.652594,
   2.101886,
   1.673642,
   1.354384,
   1.099373,
   0.8947859,
   0.7269321,
   0.5943755,
   0.4942362,
   0.4078087,
   0.3396127,
   0.2811469,
   0.2329231,
   0.1969104,
   0.1667494,
   0.1401705,
   0.1181975,
   0.09941895,
   0.08453863,
   0.07199112,
   0.06099281,
   0.05159406,
   0.04362604,
   0.03755738,
   0.0322572,
   0.02760804,
   0.02356115,
   0.02003733,
   0.01616168,
   0.01303167,
   0.010503,
   0.008461106,
   0.006813155};
   Double_t _fehx3005[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3005[43] = {
   32.80858,
   26.76786,
   20.59123,
   15.38538,
   11.56674,
   8.951936,
   7.015978,
   5.413142,
   4.209645,
   3.323897,
   2.652594,
   2.101886,
   1.673642,
   1.354384,
   1.099373,
   0.8947859,
   0.7269321,
   0.5943755,
   0.4942362,
   0.4078087,
   0.3396127,
   0.2811469,
   0.2329231,
   0.1969104,
   0.1667494,
   0.1401705,
   0.1181975,
   0.09941895,
   0.08453863,
   0.07199112,
   0.06099281,
   0.05159406,
   0.04362604,
   0.03755738,
   0.0322572,
   0.02760804,
   0.02356115,
   0.02003733,
   0.01616168,
   0.01303167,
   0.010503,
   0.008461106,
   0.006813155};
   grae = new TGraphAsymmErrors(43,_fx3005,_fy3005,_felx3005,_fehx3005,_fely3005,_fehy3005);
   grae->SetName("");
   grae->SetTitle("");

   ci = TColor::GetColor("#ff6666");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3002);

   ci = TColor::GetColor("#cc3333");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   
   TH1F *Graph_Graph3005 = new TH1F("Graph_Graph3005","",100,380,5420);
   Graph_Graph3005->SetMinimum(0.008167749);
   Graph_Graph3005->SetMaximum(569.7668);
   Graph_Graph3005->SetDirectory(0);
   Graph_Graph3005->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3005->SetLineColor(ci);
   Graph_Graph3005->GetXaxis()->SetLabelFont(42);
   Graph_Graph3005->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3005->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3005->GetXaxis()->SetTitleFont(42);
   Graph_Graph3005->GetYaxis()->SetLabelFont(42);
   Graph_Graph3005->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3005->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3005->GetYaxis()->SetTitleOffset(0);
   Graph_Graph3005->GetYaxis()->SetTitleFont(42);
   Graph_Graph3005->GetZaxis()->SetLabelFont(42);
   Graph_Graph3005->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3005->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3005->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3005);
   
   grae->Draw(", l3");
   
   Double_t _fx3006[43] = {
   800,
   900,
   1000,
   1100,
   1200,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy3006[43] = {
   485.1621,
   367.2001,
   263.8817,
   188.6818,
   135.7415,
   98.58745,
   72.32583,
   53.57184,
   40.04648,
   30.18608,
   22.91748,
   17.51997,
   13.47566,
   10.42173,
   8.100062,
   6.323773,
   4.956944,
   3.899242,
   3.077397,
   2.435993,
   1.93346,
   1.540306,
   1.225945,
   0.9791569,
   0.7834139,
   0.627741,
   0.5036664,
   0.4045769,
   0.3252906,
   0.2617513,
   0.2107641,
   0.1698193,
   0.1368799,
   0.1103879,
   0.08903561,
   0.07180615,
   0.0579232,
   0.0467275,
   0.03768939,
   0.03039015,
   0.02449322,
   0.01973149,
   0.01588843};
   Double_t _felx3006[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3006[43] = {
   32.80858,
   26.76786,
   20.59123,
   15.38538,
   11.56674,
   8.951936,
   7.015978,
   5.413142,
   4.209645,
   3.323897,
   2.652594,
   2.101886,
   1.673642,
   1.354384,
   1.099373,
   0.8947859,
   0.7269321,
   0.5943755,
   0.4942362,
   0.4078087,
   0.3396127,
   0.2811469,
   0.2329231,
   0.1969104,
   0.1667494,
   0.1401705,
   0.1181975,
   0.09941895,
   0.08453863,
   0.07199112,
   0.06099281,
   0.05159406,
   0.04362604,
   0.03755738,
   0.0322572,
   0.02760804,
   0.02356115,
   0.02003733,
   0.01616168,
   0.01303167,
   0.010503,
   0.008461106,
   0.006813155};
   Double_t _fehx3006[43] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3006[43] = {
   32.80858,
   26.76786,
   20.59123,
   15.38538,
   11.56674,
   8.951936,
   7.015978,
   5.413142,
   4.209645,
   3.323897,
   2.652594,
   2.101886,
   1.673642,
   1.354384,
   1.099373,
   0.8947859,
   0.7269321,
   0.5943755,
   0.4942362,
   0.4078087,
   0.3396127,
   0.2811469,
   0.2329231,
   0.1969104,
   0.1667494,
   0.1401705,
   0.1181975,
   0.09941895,
   0.08453863,
   0.07199112,
   0.06099281,
   0.05159406,
   0.04362604,
   0.03755738,
   0.0322572,
   0.02760804,
   0.02356115,
   0.02003733,
   0.01616168,
   0.01303167,
   0.010503,
   0.008461106,
   0.006813155};
   grae = new TGraphAsymmErrors(43,_fx3006,_fy3006,_felx3006,_fehx3006,_fely3006,_fehy3006);
   grae->SetName("");
   grae->SetTitle("");

   ci = TColor::GetColor("#ff6666");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3002);

   ci = TColor::GetColor("#cc3333");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   
   TH1F *Graph_Graph_Graph30053006 = new TH1F("Graph_Graph_Graph30053006","",100,380,5420);
   Graph_Graph_Graph30053006->SetMinimum(0.008167749);
   Graph_Graph_Graph30053006->SetMaximum(569.7668);
   Graph_Graph_Graph30053006->SetDirectory(0);
   Graph_Graph_Graph30053006->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph30053006->SetLineColor(ci);
   Graph_Graph_Graph30053006->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph30053006->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph30053006->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph30053006->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph30053006->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph30053006->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph30053006->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph30053006->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_Graph30053006->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph30053006->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph30053006->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph30053006->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph30053006->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph30053006);
   
   grae->Draw(", l3x0y0");
   
   Double_t _fx1[41] = {
   800,
   900,
   1000,
   1300,
   1400,
   1500,
   1600,
   1700,
   1800,
   1900,
   2000,
   2100,
   2200,
   2300,
   2400,
   2500,
   2600,
   2700,
   2800,
   2900,
   3000,
   3100,
   3200,
   3300,
   3400,
   3500,
   3600,
   3700,
   3800,
   3900,
   4000,
   4100,
   4200,
   4300,
   4400,
   4500,
   4600,
   4700,
   4800,
   4900,
   5000};
   Double_t _fy1[41] = {
   25.9375,
   13.9375,
   15.1875,
   4.5938,
   3.5312,
   2.9297,
   2.3906,
   2.0078,
   1.6797,
   1.4258,
   1.2109,
   1.0312,
   0.8945,
   0.7812,
   0.6914,
   0.6172,
   0.5527,
   0.5,
   0.4551,
   0.418,
   0.3857,
   0.3574,
   0.335,
   0.3164,
   0.2998,
   0.2871,
   0.2754,
   0.2656,
   0.2578,
   0.251,
   0.2441,
   0.2373,
   0.2314,
   0.2266,
   0.2217,
   0.2178,
   0.2129,
   0.21,
   0.208,
   0.207,
   0.207};
   TGraph *graph = new TGraph(41,_fx1,_fy1);
   graph->SetName("");
   graph->SetTitle("");
   graph->SetFillColor(1);
   graph->SetLineStyle(2);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,380,5420);
   Graph_Graph1->SetMinimum(0.1863);
   Graph_Graph1->SetMaximum(28.51055);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw(", l");
   TLatex *   tex = new TLatex(0.15,0.95,"");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.79,"0l, 2l categories, no VBF");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.95,0.99,"RunII, 137.2 fb^{-1}  (13 TeV)");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(42);
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.88,"CMS");
tex->SetNDC();
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.84,"Preliminary");
tex->SetNDC();
   tex->SetTextFont(52);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.55,0.54,0.98,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","95% CL upper limits","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","Observed","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","Expected","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","#pm 1 std. deviation","f");

   ci = TColor::GetColor("#00cc00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#00cc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","#pm 2 std. deviation","f");

   ci = TColor::GetColor("#ffcc00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","HVT model A","fl");

   ci = TColor::GetColor("#ff66ff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3013);

   ci = TColor::GetColor("#cc33cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","HVT model B","fl");

   ci = TColor::GetColor("#ff6666");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3002);

   ci = TColor::GetColor("#cc3333");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TH1F *Graph_copy = new TH1F("Graph_copy","",100,380,5420);
   Graph_copy->SetMinimum(0.01);
   Graph_copy->SetMaximum(5000);
   Graph_copy->SetDirectory(0);
   Graph_copy->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_copy->SetLineColor(ci);
   Graph_copy->GetXaxis()->SetTitle("m_{Z'} (GeV)");
   Graph_copy->GetXaxis()->SetRange(9,92);
   Graph_copy->GetXaxis()->SetMoreLogLabels();
   Graph_copy->GetXaxis()->SetNoExponent();
   Graph_copy->GetXaxis()->SetLabelFont(42);
   Graph_copy->GetXaxis()->SetLabelSize(0.045);
   Graph_copy->GetXaxis()->SetTitleSize(0.05);
   Graph_copy->GetXaxis()->SetTitleOffset(0.9);
   Graph_copy->GetXaxis()->SetTitleFont(42);
   Graph_copy->GetYaxis()->SetTitle("#sigma(Z') #bf{#it{#Beta}}(Z' #rightarrow Zh) (fb)");
   Graph_copy->GetYaxis()->SetMoreLogLabels();
   Graph_copy->GetYaxis()->SetNoExponent();
   Graph_copy->GetYaxis()->SetLabelFont(42);
   Graph_copy->GetYaxis()->SetLabelSize(0.045);
   Graph_copy->GetYaxis()->SetTitleSize(0.05);
   Graph_copy->GetYaxis()->SetTitleOffset(1.25);
   Graph_copy->GetYaxis()->SetTitleFont(42);
   Graph_copy->GetZaxis()->SetLabelFont(42);
   Graph_copy->GetZaxis()->SetLabelSize(0.035);
   Graph_copy->GetZaxis()->SetTitleSize(0.035);
   Graph_copy->GetZaxis()->SetTitleFont(42);
   Graph_copy->Draw("sameaxis");
   
   leg = new TLegend(0.12,0.125,0.65,0.225,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
