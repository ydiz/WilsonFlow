#include <stdlib.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace Grid {
namespace QCD {

void WF_init(int argc, char **argv, WilsonFlow_para &WF_para)
{
  po::options_description desc("GFFA options");
  desc.add_options()("help", "help message")
                    // ("steps", po::value<int>(&WF_para.steps)->default_value(50), "")
                    ("step_size", po::value<double>(&WF_para.step_size)->default_value(1.0), "")
                    ("adaptiveErrorTolerance", po::value<double>(&WF_para.adaptiveErrorTolerance)->default_value(2e-6), "")
                    ("StartTrajectory", po::value<int>(&WF_para.StartTrajectory)->default_value(0))
                    ("EndTrajectory", po::value<int>(&WF_para.EndTrajectory)->default_value(0))
                    ("TrajectoryInterval", po::value<int>(&WF_para.TrajectoryInterval)->default_value(1))
                    ("outFile", po::value<std::string>(&WF_para.outFile)->default_value("topoCharge.txt"))
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
