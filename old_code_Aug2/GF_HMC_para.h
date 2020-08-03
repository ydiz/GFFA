#pragma once

namespace Grid{
namespace QCD {

// class MyTC_para {
// public:
//   std::string type;
//   double step_size;
//   double adaptiveErrorTolerance;
//   double maxTau;
//
//   int TrajectoryStart;
//   int TrajectoryInterval;
//
//   bool saveSmearField;
//   std::string smearFieldFilePrefix;
//   std::string topoChargeOutFile;
// };


struct HMC_PARA{
  // std::string StartingType;
  // int StartingTrajectory;
  // int Thermalizations;
  // int Trajectories;
  // int mdSteps;
  // double trajL;
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

  // Measuring A
  bool measure_A;
  double fixed_P_k;
  std::vector<std::vector<int>> measure_A_coors;
  // std::vector<int> measure_A_coorsss;
  
  bool isGFFA;

  // MyTC_para tc_para;
  // GaugeModes_para gm_para;


};

std::ostream& operator<<(std::ostream& out, const HMC_PARA &HMC_para)
{
  // out << "StartingType: " << HMC_para.StartingType << std::endl;
  // out << "StartingTrajectory: " << HMC_para.StartingTrajectory << std::endl;
  // out << "Thermalizations: " << HMC_para.Thermalizations << std::endl;
  // out << "Trajectories: " << HMC_para.Trajectories << std::endl;
  // out << "mdSteps: " << HMC_para.mdSteps << std::endl;
  // out << "trajL: " << HMC_para.trajL << std::endl;
  out << "saveInterval: " << HMC_para.saveInterval << std::endl;
  out << "===================================================" << std::endl;

  out << "newHp: " << std::boolalpha << HMC_para.newHp << std::endl;
  out << "action: " << HMC_para.action << std::endl;
  out << "isGFFA: " << HMC_para.isGFFA << std::endl;
  out << "beta: " << HMC_para.beta << std::endl;
  out << "M: " << HMC_para.M << std::endl;
  out << "betaMM: " << HMC_para.betaMM << std::endl;
  out << "epsilon: " << HMC_para.epsilon << std::endl;
  out << "innerMC_N: " << HMC_para.innerMC_N << std::endl;
  out << "hb_offset: " << HMC_para.hb_offset << std::endl;
  out << "table_path: " << HMC_para.table_path << std::endl;

  out << "================ Measure A ===============" << std::endl;
  out << "measure_A: " << HMC_para.measure_A << std::endl;
  out << "fixed_P_k: " << HMC_para.fixed_P_k << std::endl;
  out << "measure_A_coors: " << HMC_para.measure_A_coors << std::endl;

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