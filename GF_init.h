namespace Grid{

// "--noMetro, --traj(n.o. of trajectories), --mdSteps, --trajL, --M, --epsilon, --innerMC_N, --hb_nsweeps, --newHp"

void GF_init(int argc, char **argv, int &noMetro, int &traj, int &mdSteps, Real &trajL, HMC_PARA &HMC_para)
{
  std::string arg;
  std::stringstream ss;

  if(GridCmdOptionExists(argv, argv+argc, "--noMetro")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--noMetro");
    ss.clear();
    ss.str(arg);
    ss >> noMetro;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--traj")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    ss.clear();
    ss.str(arg);
    ss >> traj;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--mdSteps")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--mdSteps");
    ss.clear();
    ss.str(arg);
    ss >> mdSteps;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--trajL")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--trajL");
    ss.clear();
    ss.str(arg);
    ss >> trajL;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--M")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--M");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.M;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--epsilon")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--epsilon");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.epsilon;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--hb_offset")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--hb_offset");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.hb_offset;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--hb_nsweeps")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--hb_nsweeps");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.hb_nsweeps;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--innerMC_N")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--innerMC_N");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.innerMC_N;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--newHp")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--newHp");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.newHp;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--newAction")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--newAction");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.newAction;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--SDGF")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--SDGF");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.SDGF;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--UInitEquil")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--UInitEquil");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.UInitEquil;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--hb_multi_hit")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--hb_multi_hit");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.hb_multi_hit;
  }
  if(GridCmdOptionExists(argv, argv+argc, "--UFile")){
    arg = GridCmdOptionPayload(argv, argv+argc, "--UFile");
    ss.clear();
    ss.str(arg);
    ss >> HMC_para.UFile;
  }


  HMC_para.betaMM = HMC_para.beta * HMC_para.M * HMC_para.M;
}

}
