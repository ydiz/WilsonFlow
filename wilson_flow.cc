#include "WilsonFlow.h"

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

  // LatticeComplex c(grid);
  // c = 2.0;

  // cout << norm2(c) << endl;

  // GridParallelRNG           pRNG(grid);   pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  LatticeGaugeField U(grid);
  U = 1.0;

  Lattice<iVector<iScalar<iScalar<vRealD>>, Nd> > U_norm2(grid);

  parallel_for(int ss=0;ss<U_norm2._grid->oSites();ss++){
    for(int mu=0; mu<Nd; mu++) {
      U_norm2[ss](mu)()() = 0.;
      for(int c1=0;c1<Nc;c1++)
        for(int c2=0;c2<Nc;c2++)
            U_norm2[ss](mu)()() += toReal(U[ss](mu)()(c1, c2) * U[ss](mu)()(c1, c2));
    }
  }

  typedef typename decltype(U_norm2)::scalar_type scalar_type;

  // U_norm2[10](2)()() = 9.0;

  double max_val = 0.;
  #pragma omp parallel for reduction(max : max_val)
  for(int ss=0;ss<U_norm2._grid->oSites();ss++)
  {
    for(int mu=0; mu<Nd; mu++) {
      scalar_type *sobj = (scalar_type *)& U_norm2[ss](mu)()();
      for(int idx=0; idx<U_norm2._grid->Nsimd(); ++idx)
        if( *(sobj + idx) > max_val)
            max_val = *(sobj + idx);
    }

  }

  #ifndef GRID_COMMS_NONE
  MPI_Allreduce(MPI_IN_PLACE, &max_val, 1, MPI_DOUBLE, MPI_MAX, U_norm2.communicator)
  #endif

  cout << max_val << endl;

  Grid_finalize();

}
