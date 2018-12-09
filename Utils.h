#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

  
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
