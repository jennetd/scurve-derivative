#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <cmath>
#include <iostream>

TRandom3 myDice;
TGraph *myGraph = nullptr;
TH1D *myHisto = nullptr;

int nTrials = 1e3;

// S-curve goes from minCurve to maxCurve
const double minCurve = 0;
const double maxCurve = 100;

// number of steps
const int nThresholds = 40;

int nEvents = 100;

class Estimates {
public:
  double fitMean;
  double fitNoise;
  double diffMean;
  double diffNoise;
  double mean;
  double noise;
};

// Makes a realistic S-curve random mean/noise
void generateRandomParameters(double &mean, double &noise) {
  mean = (minCurve + maxCurve) / 2. + (maxCurve - minCurve) / 6 * myDice.Rndm();
  noise = (maxCurve - minCurve) * (myDice.Rndm() + 1) / 20.;
}

// Function to fit S-curve with
TF1 *myFitFunction = new TF1(
			     "myFunc", "(TMath::Erf((x - [0])/sqrt(2)/[1])+1)*0.5", minCurve, maxCurve);

void printResult(Estimates est) {
  std::cout << "fit method:" << std::endl;
  std::cout << "p0 = " << est.fitMean << " / " << est.mean << std::endl;
  std::cout << "p1 = " << est.fitNoise << " / " << est.noise << std::endl;
  std::cout << "histogram method:" << std::endl;
  std::cout << "p0 = " << est.diffMean << " / " << est.mean << std::endl;
  std::cout << "p1 = " << est.diffNoise << " / " << est.noise << std::endl;
}

// makes S-curve points
Estimates makeScurve(bool randomParameters = false, bool debugPrint = false,
                     bool draw = true) {
  double mean;
  double noise;
  
  if (randomParameters) {
    generateRandomParameters(mean, noise);
  } 
  // value of mean and noise are specified here
  else {
    mean = (maxCurve + minCurve) / 2;
    noise = (maxCurve - minCurve) / 10;
  }

  // make a graph object
  if (myGraph)
    delete myGraph;
  myGraph = new TGraph();

  const double aStep = (maxCurve - minCurve) / nThresholds;
  double aValue;

  // a small simulation to get "data"
  for (double threshold = minCurve; threshold < maxCurve; threshold += aStep) {
    double efficiency = 0;
    for (int i = 0; i < nEvents; ++i) {
      aValue = myDice.Gaus(mean, noise);
      if (aValue < threshold)
        efficiency++;
    }
    efficiency /= nEvents;
    myGraph->SetPoint(myGraph->GetN(), threshold, efficiency);
  }

  if (draw) {
    myGraph->SetMarkerStyle(8);
    myGraph->Draw("ap");
  }

  // Fit method
  myFitFunction->SetParameter(0, mean);
  myFitFunction->SetParameter(1, noise);
  myGraph->Fit(myFitFunction, "Q");

  // Differential method
  if (myHisto)
    delete myHisto;
  myHisto =
    new TH1D("myHisto", "my histogram", nThresholds - 1, minCurve, maxCurve);
  double x1, x2, y1, y2;
  for (int i = 0; i < nThresholds - 1; ++i) {
    myGraph->GetPoint(i, x1, y1);
    myGraph->GetPoint(i + 1, x2, y2);
    double x = aStep * (0.5 + i);
    myHisto->Fill(x, y2 - y1);
  }

  Estimates result;
  result.fitMean = myFitFunction->GetParameter(0);
  result.fitNoise = myFitFunction->GetParameter(1);
  result.diffMean = myHisto->GetMean();
  result.diffNoise = myHisto->GetRMS();
  result.mean = mean;
  result.noise = noise;

  if (debugPrint)
    printResult(result);

  return result;
}

TH1D *fitMeanHisto = nullptr;
TH1D *fitNoiseHisto = nullptr;
TH1D *diffMeanHisto = nullptr;
TH1D *diffNoiseHisto = nullptr;

void methodComparison() {
  Estimates est;

  if (fitMeanHisto)
    delete fitMeanHisto;
  if (fitNoiseHisto)
    delete fitNoiseHisto;
  if (diffMeanHisto)
    delete diffMeanHisto;
  if (diffNoiseHisto)
    delete diffNoiseHisto;

  const int nBins = 100;
  const double minNoise = 0;
  const double maxNoise = 20;
  const double minMean = 40;
  const double maxMean = 60;

  fitMeanHisto =
    new TH1D("fitMeanHisto", "fitMeanHisto", nBins, minMean, maxMean);
  fitNoiseHisto =
    new TH1D("fitNoiseHisto", "fitNoiseHisto", nBins, minNoise, maxNoise);

  diffMeanHisto =
    new TH1D("diffMeanHisto", "diffMeanHisto", nBins, minMean, maxMean);
  diffNoiseHisto =
    new TH1D("diffNoiseHisto", "diffNoiseHisto", nBins, minNoise, maxNoise);

  for (int i = 0; i < nTrials; ++i) {
    est = makeScurve();
    fitMeanHisto->Fill(est.fitMean);
    diffMeanHisto->Fill(est.diffMean);
    fitNoiseHisto->Fill(est.fitNoise);
    diffNoiseHisto->Fill(est.diffNoise);
  }

  std::cout << "MEAN COMPARISON" << std::endl;
  std::cout << "diff: " << diffMeanHisto->GetMean()
            << " with σ= " << diffMeanHisto->GetRMS() << std::endl;
  std::cout << "fit : " << fitMeanHisto->GetMean()
            << " with σ= " << fitMeanHisto->GetRMS() << std::endl;

  std::cout << "NOISE COMPARISON" << std::endl;
  std::cout << "diff: " << diffNoiseHisto->GetMean()
            << " with σ= " << diffNoiseHisto->GetRMS() << std::endl;
  std::cout << "fit : " << fitNoiseHisto->GetMean()
            << " with σ= " << fitNoiseHisto->GetRMS() << std::endl;
}
