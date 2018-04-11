# this makefile is written by Hongguan Wu in March, 2018
# makefile should be arranged according to your main program

objects = global_module_constants.o global_module_parameters.o \
      global_module_energy.o local_module_mfmt.o local_module_parameters.o \
      main.o sub_system_def.o sub_comp_external_pot.o sub_miu_att.o\
	  sub_comp_bulk_pot.o sub_comp_hs_pot.o sub_comp_att_pot.o


cc = gfortran

edit:$(objects)
	$(cc) -o edit $(objects)


global_module_constants.o:global_module_constants.f90
	$(cc) -c global_module_constants.f90

global_module_parameters.o:global_module_constants.o global_module_parameters.f90
	$(cc) -c global_module_parameters.f90

global_module_energy.o:global_module_parameters.o global_module_energy.f90
	$(cc) -c global_module_energy.f90

local_module_parameters.o:local_module_parameters.f90
	$(cc) -c local_module_parameters.f90

local_module_mfmt.o:local_module_mfmt.f90
	$(cc) -c local_module_mfmt.f90

main.o:global_module_energy.o main.f90
	$(cc) -c main.f90

sub_system_def.o:global_module_energy.o sub_system_def.f90
	$(cc) -c sub_system_def.f90

sub_comp_external_pot.o:global_module_energy.o sub_comp_external_pot.f90
	$(cc) -c sub_comp_external_pot.f90

sub_comp_bulk_pot.o:global_module_energy.o local_module_parameters.o sub_comp_bulk_pot.f90
	$(cc) -c sub_comp_bulk_pot.f90

sub_miu_att.o:global_module_energy.o sub_miu_att.f90
	$(cc) -c sub_miu_att.f90

sub_comp_hs_pot.o:global_module_energy.o local_module_mfmt.o sub_comp_hs_pot.f90
	$(cc) -c sub_comp_hs_pot.f90

sub_comp_att_pot.o:global_module_energy.o sub_comp_att_pot.f90
	$(cc) -c sub_comp_att_pot.f90

clean:
	rm  $(objects) edit
