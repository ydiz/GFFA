namespace Grid{

struct HMC_PARA{
  std::string StartingType;
  int StartingTrajectory;
  int Thermalizations;
  int Trajectories;
  int mdSteps;
  double trajL;
  int saveInterval;

  bool newHp;
  bool newAction;
  double beta;
  double M;
  double betaMM;
  double epsilon;
  int innerMC_N;
  int hb_offset;
  // int hb_nsweeps;
  int hb_multi_hit;
  // bool SDGF; // SteepestDescentGaugeFix
  // bool UInitEquil; //use equilibrium U
  std::string UFile;

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
  out << "newAction: " << HMC_para.newAction << std::endl;
  out << "beta: " << HMC_para.beta << std::endl;
  out << "M: " << HMC_para.M << std::endl;
  out << "betaMM: " << HMC_para.betaMM << std::endl;
  out << "epsilon: " << HMC_para.epsilon << std::endl;
  out << "innerMC_N: " << HMC_para.innerMC_N << std::endl;
  out << "hb_offset: " << HMC_para.hb_offset << std::endl;
  // out << "hb_nsweeps: " << HMC_para.hb_nsweeps << std::endl; //inteval
  // out << "SDGF: " << HMC_para.SDGF << std::endl;
  // out << "UInitEquil: " << HMC_para.UInitEquil << std::endl;
  out << "hb_multi_hit: " << HMC_para.hb_multi_hit << std::endl;
  return out;
}

}
