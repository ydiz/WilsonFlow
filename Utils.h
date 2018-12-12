#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

std::vector<double> operator+(const std::vector<double> &x1, const std::vector<double> &x2) {
  std::vector<double> ret(x1.size());
  for(int i=0; i<ret.size(); ++i) ret[i] = x1[i] + x2[i];
  return ret;
}

std::vector<double> operator*(double x, const std::vector<double> &vec) {
  std::vector<double> ret(vec.size());
  for(int i=0; i<ret.size(); ++i) ret[i] = x * vec[i];
  return ret;
}


template<class T>
void writeVector(const std::vector<T> &vec, const std::string &filename, int rank) {
  if(rank == 0) {
    std::ofstream outFile(filename, std::ofstream::app);
    for(const T& x : vec) outFile << x << " ";
    outFile << "\n";
    outFile.close();
  }
#ifdef GRID_COMMS_MPI3
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
