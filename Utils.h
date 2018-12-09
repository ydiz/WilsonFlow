#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

template<class T>
void writeVector(const std::vector<T> &vec, const std::string &filename) {

#ifdef GRID_COMMS_MPI3
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
#endif

    std::ofstream outFile(filename, std::ofstream::app);
    for(const T& x : vec) outFile << x << " ";
    outFile << "\n";
    outFile.close();

#ifdef GRID_COMMS_MPI3
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}


void readField(LatticeGaugeField &U, const std::string &filename)
{
		FieldMetaData header;
			NerscIO::readConfiguration(U,header,filename);
}

void writeField(LatticeGaugeField &U, const std::string &filename)
{
		NerscIO::writeConfiguration(U,filename,0,0);
}



}
}
