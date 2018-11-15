// g++ GF_para.cc -lboost_program_options
#include <stdlib.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace Grid {

void init(int argc, char **argv, HMC_PARA &hmc_para)
{
  po::options_description desc("GFFA options");
  desc.add_options()("help", "help message")
                    ("StartingType", po::value<std::string>(&hmc_para.StartingType)->default_value("ColdStart"), "Stariing configuration. It can be HotStart, ColdStart, TepidStart, or CheckpointStart.")
                    ("StartingTrajectory", po::value<int>(&hmc_para.StartingTrajectory), "If StartingType is CheckpointStart, ckpoint_lat.xx and ckpoint_rng.xx corresponding to StartingTrajectory will be loaded.")
                    ("Thermalizations", po::value<int>(&hmc_para.Thermalizations)->default_value(200), "Number of trajectories without Metropolis test.")
                    ("Trajectories", po::value<int>(&hmc_para.Trajectories)->default_value(0), "Number of trajectories after those without Metropolis test. p.s. At the moment Metropolis is disabled")
                    ("mdSteps", po::value<int>(&hmc_para.mdSteps)->default_value(20), "Number of MD steps within each trajectory.")
                    ("trajL", po::value<double>(&hmc_para.trajL)->default_value(1.0), "Trajectory length.")
                    ("saveInterval", po::value<int>(&hmc_para.saveInterval)->default_value(50), "Save interval for checker pointers.")
                    ("newAction", po::value<bool>(&hmc_para.newAction)->default_value(true), "Determine whether use gauge-fixing action.")
                    ("newHp", po::value<bool>(&hmc_para.newHp)->default_value(false), "Determine whether use Fourier-accelerated kinetic energy term.")
                    ("beta", po::value<double>(&hmc_para.beta)->default_value(5.6), "beta")
                    ("M", po::value<double>(&hmc_para.M)->default_value(1.0), "M for soft gauge-fixing.")
                    ("epsilon", po::value<double>(&hmc_para.epsilon)->default_value(0.2), "Infrared regulator")
                    ("hb_offset", po::value<int>(&hmc_para.hb_offset)->default_value(100), "number of heatbath sweeps to reach equilibrium")
                    ("innerMC_N", po::value<int>(&hmc_para.innerMC_N)->default_value(100), "number of heatbath sweeps to calculated inner Monte Carlo")
                    // ("hb_nsweeps", po::value<int>(&hmc_para.hb_nsweeps)->default_value(1), "number of heatbath sweeps for within each innerMC_N")
                    ("hb_multi_hit", po::value<int>(&hmc_para.hb_multi_hit)->default_value(1),"heatbath multiple hits")
                    // ("UInitEquil", po::value<bool>(),"") // delete this parameter, replace it with !UFileName.empty();
                    ("UFileName", po::value<std::string>(&hmc_para.UFile), "If starting type if not CheckpointStart and UFileName is not empty, gauge configuration with corresponding filename with be loaded")
                    ;

  po::variables_map vm;
  // po::store(po::parse_command_line(argc, argv, desc), vm); // command line options have higher priority
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options
  po::store(po::parse_config_file<char>("GFFA.ini", desc), vm);
  po::notify(vm);

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

  hmc_para.betaMM = hmc_para.beta * hmc_para.M * hmc_para.M;
}

}
