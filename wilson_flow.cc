#include "Utils.h"
#include "WilsonFlow.h"
#include "parameters.h"
#include "TopologicalCharge.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  LatticeGaugeField U(grid);
  std::string filename ("ckpoint_lat.1000");
  readField(U, filename);

  WilsonFlow_para WF_para;
  WF_init(argc, argv, WF_para);

  MyWilsonFlow<PeriodicGimplR> WF(WF_para.steps, WF_para.step_size, WF_para.adaptiveErrorTolerance, WF_para.measure_interval);

  LatticeGaugeField Uflow(U._grid);
  WF.smear_adaptive(Uflow, U);

  // cout << WilsonLoops<PeriodicGimplR>::TopologicalCharge(U) << endl;
  // // cout << topologicalCharge(U) << endl;
  // std::vector<double> ret = timeSliceTopologicalCharge(U);
  // cout << ret << endl;
  //
  // double sum=0;
  // for(auto x: ret) sum += x;
  // cout << sum << endl;


  Grid_finalize();

}
