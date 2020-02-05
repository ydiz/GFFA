namespace Grid{
namespace QCD {

struct HMC_PARA{
  std::string StartingType;
  int StartingTrajectory;
  int Thermalizations;
  int Trajectories;
  int mdSteps;
  double trajL;
  int saveInterval;

  bool newHp;
  std::string action;
  double beta;
  double M;
  double betaMM;
  double epsilon;
  int innerMC_N;
  int hb_offset;
  std::string table_path;
  std::string UFile;

  // MyTC_para tc_para;

  // // Topological Charge
  // int TC_interval;
  // bool TC_do_smearing;
  // // int TC_Smearing_steps; // does not matter
  // double TC_Smearing_step_size;
  // int TC_Smearing_meas_interval;
  // double TC_Smearing_maxTau;

  // HMC_PARA(): beta(5.6), M(1.0), epsilon(0.2), hb_offset(1000), innerMC_N(10000), hb_nsweeps(1),
  //             newHp(false), newAction(true), SDGF(false), UInitEquil(true), hb_multi_hit(1) {} //betaMM is initialzed in GF_init()
};

std::ostream& operator<<(std::ostream& out, const HMC_PARA &HMC_para)
{
  out << "StartingType: " << HMC_para.StartingType << std::endl;
  out << "StartingTrajectory: " << HMC_para.StartingTrajectory << std::endl;
  out << "Thermalizations: " << HMC_para.Thermalizations << std::endl;
  out << "Trajectories: " << HMC_para.Trajectories << std::endl;
  out << "mdSteps: " << HMC_para.mdSteps << std::endl;
  out << "trajL: " << HMC_para.trajL << std::endl;
  out << "saveInterval: " << HMC_para.saveInterval << std::endl;
  out << "===================================================" << std::endl;

  out << "newHp: " << std::boolalpha << HMC_para.newHp << std::endl;
  out << "action: " << HMC_para.action << std::endl;
  out << "beta: " << HMC_para.beta << std::endl;
  out << "M: " << HMC_para.M << std::endl;
  out << "betaMM: " << HMC_para.betaMM << std::endl;
  out << "epsilon: " << HMC_para.epsilon << std::endl;
  out << "innerMC_N: " << HMC_para.innerMC_N << std::endl;
  out << "hb_offset: " << HMC_para.hb_offset << std::endl;
  out << "table_path: " << HMC_para.table_path << std::endl;

  // out << "================ Topological Charge ===============" << std::endl;
  // out << "TC_interval: " << HMC_para.TC_interval << std::endl;
  // out << "TC_do_smearing: " << HMC_para.TC_do_smearing << std::endl;
  // // out << "TC_Smearing_steps: " << HMC_para.TC_Smearing_steps << std::endl;
  // out << "TC_Smearing_step_size: " << HMC_para.TC_Smearing_step_size << std::endl;
  // out << "TC_Smearing_meas_interval: " << HMC_para.TC_Smearing_meas_interval << std::endl;
  // out << "TC_Smearing_maxTau: " << HMC_para.TC_Smearing_maxTau << std::endl;

  return out;
}



}}
