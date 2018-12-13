#include <stdlib.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace Grid {
namespace QCD {

void WF_init(int argc, char **argv, WilsonFlow_para &WF_para)
{
  po::options_description desc("GFFA options");
  desc.add_options()("help", "help message")
                    ("lat", po::value<std::string>()->default_value("8.8.8.16"))
                    ("step_size", po::value<double>(&WF_para.step_size)->default_value(1.0), "")
                    ("adaptiveErrorTolerance", po::value<double>(&WF_para.adaptiveErrorTolerance)->default_value(2e-6), "")
                    ("Nstep", po::value<int>(&WF_para.Nstep)->default_value(0))
                    ("StartTrajectory", po::value<int>(&WF_para.StartTrajectory)->default_value(0))
                    ("EndTrajectory", po::value<int>(&WF_para.EndTrajectory)->default_value(0))
                    ("TrajectoryInterval", po::value<int>(&WF_para.TrajectoryInterval)->default_value(1))
                    ("inFilePrefix", po::value<std::string>(&WF_para.inFilePrefix)->default_value("ckpoint_lat"))
                    ("doSmear", po::value<bool>(&WF_para.doSmear)->default_value(true))
                    ("saveSmearField", po::value<bool>(&WF_para.saveSmearField)->default_value(false))
                    ("smearFieldFilePrefix", po::value<std::string>(&WF_para.smearFieldFilePrefix)->default_value("ckpoint_lat_smear"))
                    ("calculateTopoCharge", po::value<bool>(&WF_para.calculateTopoCharge)->default_value(true))
                    ("topoChargeOutFile", po::value<std::string>(&WF_para.topoChargeOutFile)->default_value("topoCharge.txt"))
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options
  po::store(po::parse_config_file<char>("WF.ini", desc), vm);
  po::notify(vm);

  WF_para.lat.resize(0);
  std::stringstream lat_ss(vm["lat"].as<std::string>());
  std::string segment;
  while(std::getline(lat_ss, segment, '.'))
  {
     WF_para.lat.push_back(std::stoi(segment));
  }
  assert(WF_para.lat.size() == 4);

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

}

}
}
