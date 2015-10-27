#log 		  	Si_sw.lammps
log			Si_sw_wfnho.lammps

units			metal
atom_style      	atomic

boundary		p p p
read_data       	Si100v.out
mass 			1 28.0855
	
#pair_style		sw
#pair_coeff		* * /homer2/lmhale/Desktop/lammps-5Jun10/potentials/Si.sw Si 
pair_style		sw/wfnho
pair_coeff		* * /homer2/lmhale/Desktop/lammps-5Jun10/potentials/SiO.sw_wfnho Si(a)

compute			peatom all pe/atom
pair_modify		shift yes

timestep        	0.001
fix 			boing all nvt temp 300.0 300.0 10.0 
thermo_style    	custom lx ly lz pxx pyy pzz pe 
thermo          	1	
thermo_modify   	flush yes
#dump 			dumpy all custom 1 atomsw.txt id type x y z c_peatom	
dump 			dumpy all custom 1 atomsw_wfnho.txt id type x y z c_peatom	

run			1000
