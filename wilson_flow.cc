#include "Utils.h"
#include "WilsonFlow.h"
#include "parameters.h"
#include "TopologicalCharge.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);
  GridLogLayout();

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());

  WilsonFlow_para WF_para;
  WF_init(argc, argv, WF_para);

  MyWilsonFlow<PeriodicGimplR> WF(WF_para.step_size, WF_para.adaptiveErrorTolerance);

  std::string inFile;
  std::string outFile(WF_para.topoChargeOutFile);
  LatticeGaugeField U(grid);
  LatticeGaugeField Uflow(U._grid);
  for(int i=WF_para.StartTrajectory; i<=WF_para.EndTrajectory; i+=WF_para.TrajectoryInterval)
  {
    inFile = WF_para.inFilePrefix + "." + std::to_string(i);
    readField(U, inFile);
    if(WF_para.doSmear) {
      WF.smear_adaptive(Uflow, U);
      if(WF_para.saveSmearField) writeField(Uflow, WF_para.smearFieldFilePrefix + "." + std::to_string(i));
    }
    else {
      Uflow = U;
    }
    if(WF_para.calculateTopoCharge) {
      std::vector<double> topoCharge = timeSliceTopologicalCharge(Uflow);
      writeVector(topoCharge, outFile, U._grid->ThisRank());
    }
  }

  Grid_finalize();

}
