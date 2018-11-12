# ./Test_hmc --M 1.0 --threads 64 --hb_offset 0 --innerMC_N 10
# ./Test_hmc --M 1.0 --threads 64 --UFile ckpoint/ckpoint_lat.120 --innerMC_N 1 > 10.24_M1.0_traj120-180.txt
# ./Test_hmc --M 1.0 --threads 64 --startingType CheckpointStart --StartTrajectory 60  --saveInterval 1000  --innerMC_N 1 > 10.24_M1.0_traj120-180.txt
#./Test_hmc --M 1.0 --threads 64 --startingType CheckpointStart --StartingTrajectory 120 --saveInterval 1000  --innerMC_N 1000 --hb_offset 100 --traj 999 > 10.25_M1.0_innerMC1000_traj120-xx.txt

./Test_hmc --M 1.0 --startingType ColdStart --saveInterval 1000  --innerMC_N 3 --hb_offset 0 --traj 0 --UInitEquil 0 --noMetro 3
