// g++ GF_para.cc -lboost_program_options
#include <stdlib.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;


template<typename T>
void str2vector(const std::string &str, std::vector<T> &vec)
{
  vec.clear();
  std::stringstream ss(str);
  int i;
  while (ss >> i){
    vec.push_back(i);
    if(std::ispunct(ss.peek()))
      ss.ignore();
  }
  return;
}


namespace Grid {
namespace QCD {

void init(int argc, char **argv, HMC_PARA &hmc_para)
{
  std::string measure_A_coor_str, meas_taus_str;
  po::options_description desc("GFFA options");
  desc.add_options()("help", "help message")
                    ("StartingType", po::value<std::string>(&hmc_para.StartingType)->default_value("ColdStart"), "Stariing configuration. It can be HotStart, ColdStart, TepidStart, or CheckpointStart.")
                    ("StartingTrajectory", po::value<int>(&hmc_para.StartingTrajectory), "If StartingType is CheckpointStart, ckpoint_lat.xx and ckpoint_rng.xx corresponding to StartingTrajectory will be loaded.")
                    ("Thermalizations", po::value<int>(&hmc_para.Thermalizations)->default_value(200), "Number of trajectories without Metropolis test.")
                    ("Trajectories", po::value<int>(&hmc_para.Trajectories)->default_value(0), "Number of trajectories after those without Metropolis test. p.s. At the moment Metropolis is disabled")
                    ("mdSteps", po::value<int>(&hmc_para.mdSteps)->default_value(20), "Number of MD steps within each trajectory.")
                    ("trajL", po::value<double>(&hmc_para.trajL)->default_value(1.0), "Trajectory length.")
                    ("saveInterval", po::value<int>(&hmc_para.saveInterval)->default_value(50), "Save interval for checker pointers.")
                    ("action", po::value<std::string>(&hmc_para.action)->default_value("Wilson"), "Action name; available choices: Wilson, GF_Wilson, DBW2, GF_DBW2")
                    ("newHp", po::value<bool>(&hmc_para.newHp)->default_value(false), "Determine whether use Fourier-accelerated kinetic energy term.")
                    ("beta", po::value<double>(&hmc_para.beta)->default_value(5.6), "beta")
                    ("M", po::value<double>(&hmc_para.M)->default_value(1.0), "M for soft gauge-fixing.")
                    ("epsilon", po::value<double>(&hmc_para.epsilon)->default_value(0.2), "Infrared regulator")
                    ("hb_offset", po::value<int>(&hmc_para.hb_offset)->default_value(100), "number of heatbath sweeps to reach equilibrium")
                    ("innerMC_N", po::value<int>(&hmc_para.innerMC_N)->default_value(100), "number of heatbath sweeps to calculated inner Monte Carlo")
                    ("table_path", po::value<std::string>(&hmc_para.table_path)->default_value("."),"path of integral look up table")
                    ("UFileName", po::value<std::string>(&hmc_para.UFile), "If starting type if not CheckpointStart and UFileName is not empty, gauge configuration with corresponding filename with be loaded")

                    ("measure_A", po::value<bool>(&hmc_para.measure_A)->default_value(false))
                    ("fixed_P_k", po::value<double>(&hmc_para.fixed_P_k)->default_value(0.5))
                    ("measure_A_coors", po::value<std::string>(&measure_A_coor_str)->default_value(""))

                    ("TC.type", po::value<std::string>(&hmc_para.tc_para.type)->default_value("fixed_taus"), "")
                    ("TC.step_size", po::value<double>(&hmc_para.tc_para.step_size)->default_value(1.0), "")
                    ("TC.adaptiveErrorTolerance", po::value<double>(&hmc_para.tc_para.adaptiveErrorTolerance)->default_value(2e-6), "")
                    ("TC.meas_taus", po::value<std::string>(&meas_taus_str)->default_value(""), "")
                    ("TC.TrajectoryStart", po::value<int>(&hmc_para.tc_para.TrajectoryStart)->default_value(20))
                    ("TC.TrajectoryInterval", po::value<int>(&hmc_para.tc_para.TrajectoryInterval)->default_value(1))
                    // ("TC.topoChargeOutFile", po::value<std::string>(&hmc_para.tc_para.topoChargeOutFile)->default_value("topoCharge.txt"))
                    ;

  po::variables_map vm;
  // po::store(po::parse_command_line(argc, argv, desc), vm); // command line options have higher priority
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options
  po::store(po::parse_config_file<char>("GFFA.ini", desc), vm);
  po::notify(vm);

  // Wison flow times to measure energy
  str2vector(meas_taus_str, hmc_para.tc_para.meas_taus);

  // meausre_A coors
  std::stringstream ss(measure_A_coor_str);
  std::string tmp;
  while(ss >> tmp) {
    std::vector<int> vec;
    GridCmdOptionIntVector(tmp, vec);
    assert(vec.size()==4);
    hmc_para.measure_A_coors.push_back(vec);
    hmc_para.gm_para.coors.push_back(vec);
  }

  if(hmc_para.measure_A) {
    hmc_para.Thermalizations = 1; // If measure_A, only run one trajectory
    hmc_para.tc_para.TrajectoryStart = 999999; // If measure_A, do not peroform Wilson flow and measure Topological charge
  }

  hmc_para.isGFFA = hmc_para.action.substr(0,2) == "GF";

  // std::cout << "zyd Warning: there is a discrepancy in evolution time between cps and old version Grid." << std::endl;
  // std::cout << "In terms of cps, your trajectory length is " <<  hmc_para.trajL << std::endl;
  //
  // hmc_para.trajL = hmc_para.trajL * std::sqrt(2);
  // std::cout << "In terms of old version Grid, your trajectory length is " <<  hmc_para.trajL << std::endl;

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

  hmc_para.betaMM = hmc_para.beta * hmc_para.M * hmc_para.M;
}

}}
