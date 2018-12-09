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

  std::string filename;
  std::string outFile(WF_para.outFile);
  LatticeGaugeField U(grid);
  LatticeGaugeField Uflow(U._grid);
  for(int i=WF_para.StartTrajectory; i<=WF_para.EndTrajectory; i+=WF_para.TrajectoryInterval)
  {
    filename = "ckpoint_lat." + std::to_string(i);
    readField(U, filename);
    WF.smear_adaptive(Uflow, U);
    std::vector<double> topoCharge = timeSliceTopologicalCharge(Uflow);
    writeVector(topoCharge, outFile, U._grid->ThisRank());
  }

  Grid_finalize();

}
