#include <stdlib.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace Grid {
namespace QCD {

void WF_init(int argc, char **argv, WilsonFlow_para &WF_para)
{
  po::options_description desc("GFFA options");
  desc.add_options()("help", "help message")
                    ("steps", po::value<int>(&WF_para.steps)->default_value(50), "")
                    ("step_size", po::value<double>(&WF_para.step_size)->default_value(1.0), "parameter for smearing")
                    ("adaptiveErrorTolerance", po::value<double>(&WF_para.adaptiveErrorTolerance)->default_value(2e-6), "Wilson flow integration steps for calculating topological charge")
                    ("measure_interval", po::value<int>(&WF_para.measure_interval)->default_value(50), "Wilson flow integration steps for calculating topological charge")
                    ("maxTau", po::value<double>(&WF_para.maxTau)->default_value(2.0), "parameter for smearing")
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options
  po::store(po::parse_config_file<char>("WF.ini", desc), vm);
  po::notify(vm);

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

}

}
}
