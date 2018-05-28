namespace Grid{

struct HMC_PARA{
  Real beta;
  Real M;
  Real betaMM;
  Real epsilon;
  int innerMC_N;
  int hb_offset;
  int hb_nsweeps;
  bool newHp;
  bool newAction;
  bool SDGF; // SteepestDescentGaugeFix
  bool UInitEquil; //use equilibrium U

  HMC_PARA(): beta(5.6), M(3.0), epsilon(0.2), hb_offset(300), innerMC_N(500), hb_nsweeps(1),
              newHp(false), newAction(true), SDGF(true), UInitEquil(true)
  {
    betaMM = beta * M *M;
  }
};

std::ostream& operator<<(std::ostream& out, const HMC_PARA &HMC_para)
{
  out << "beta: " << HMC_para.beta << std::endl;
  out << "M: " << HMC_para.M << std::endl;
  out << "epsilon: " << HMC_para.epsilon << std::endl;
  out << "innerMC_N: " << HMC_para.innerMC_N << std::endl;
  out << "hb_offset: " << HMC_para.hb_offset << std::endl;
  out << "hb_nsweeps: " << HMC_para.hb_nsweeps << std::endl; //inteval
  out << "newHp: " << HMC_para.newHp << std::endl;
  out << "newAction: " << HMC_para.newAction << std::endl;
  out << "SDGF: " << HMC_para.SDGF << std::endl;
  out << "UInitEquil: " << HMC_para.UInitEquil << std::endl;
  return out;
}

}
