#
# Customized Makefile for RRTMG_SW GCM package
# This builds the package into a shared library that can be linked into a coupled system
#

INSTALLDIR=/scratch4/NCEPDEV/ocean/noscrub/Fei.Liu/RRTMG/rrtmg_sw_v4.02/gcm_model/install
MACHINE=theia

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif
include         $(ESMFMKFILE)
ESMF_INC        = $(ESMF_F90COMPILEPATHS)
ESMF_LIB        = $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)

VERSION=4.0.2
FC=mpiifort
FCFLAGS=-O2 -g -real-size 64 -integer-size 64 -fPIC $(ESMF_INC)
LINKFLAGS=-shared
AR=ar
ARFLAGS=-r

TARGET=librrtmg_sw.$(VERSION).so
SOLIB=librrtmg_sw.so
INCLUDE=-I. -Imod -I$(NETCDF)/include
LIB=-L$(NETCDF)/lib -lnetcdf
LIBRARY=librrtmg_sw.a

MODULES_DIR=modules
SRC_DIR=src
OBJ_DIR=obj

TARGET_OBJ=$(OBJ_DIR)/parrrsw.o	 $(OBJ_DIR)/rrsw_kg16.o	$(OBJ_DIR)/rrsw_kg20.o  $(OBJ_DIR)/rrsw_kg24.o  $(OBJ_DIR)/rrsw_kg28.o   $(OBJ_DIR)/rrsw_tbl.o $(OBJ_DIR)/rrsw_aer.o  $(OBJ_DIR)/rrsw_kg17.o	$(OBJ_DIR)/rrsw_kg21.o  $(OBJ_DIR)/rrsw_kg25.o  $(OBJ_DIR)/rrsw_kg29.o   $(OBJ_DIR)/rrsw_vsn.o $(OBJ_DIR)/rrsw_cld.o  $(OBJ_DIR)/rrsw_kg18.o	$(OBJ_DIR)/rrsw_kg22.o  $(OBJ_DIR)/rrsw_kg26.o  $(OBJ_DIR)/rrsw_ncpar.o  $(OBJ_DIR)/rrsw_wvn.o $(OBJ_DIR)/parkind.o		   $(OBJ_DIR)/rrsw_con.o  $(OBJ_DIR)/rrsw_kg19.o	$(OBJ_DIR)/rrsw_kg23.o  $(OBJ_DIR)/rrsw_kg27.o  $(OBJ_DIR)/rrsw_ref.o $(OBJ_DIR)/rrtmg_sw_cldprmc.o	$(OBJ_DIR)/rrtmg_sw_k_g.o	  	$(OBJ_DIR)/rrtmg_sw_spcvmc.o  $(OBJ_DIR)/rrtmg_sw_vrtqdr.o $(OBJ_DIR)/mcica_random_numbers.o  $(OBJ_DIR)/rrtmg_sw_cldprop.o	$(OBJ_DIR)/rrtmg_sw_rad.o	  $(OBJ_DIR)/rrtmg_sw_reftra.o	$(OBJ_DIR)/rrtmg_sw_spcvrt.o $(OBJ_DIR)/mcica_subcol_gen_sw.o   $(OBJ_DIR)/rrtmg_sw_init.o	  $(OBJ_DIR)/rrtmg_sw_setcoef.o	$(OBJ_DIR)/rrtmg_sw_taumol.o $(OBJ_DIR)/rrtmg_cap.o

#TARGET with nc read
#TARGET_NC_OBJ=$(OBJ_DIR)/parrrsw.o	 $(OBJ_DIR)/rrsw_kg16.o	$(OBJ_DIR)/rrsw_kg20.o  $(OBJ_DIR)/rrsw_kg24.o  $(OBJ_DIR)/rrsw_kg28.o   $(OBJ_DIR)/rrsw_tbl.o $(OBJ_DIR)/rrsw_aer.o  $(OBJ_DIR)/rrsw_kg17.o	$(OBJ_DIR)/rrsw_kg21.o  $(OBJ_DIR)/rrsw_kg25.o  $(OBJ_DIR)/rrsw_kg29.o   $(OBJ_DIR)/rrsw_vsn.o $(OBJ_DIR)/rrsw_cld.o  $(OBJ_DIR)/rrsw_kg18.o	$(OBJ_DIR)/rrsw_kg22.o  $(OBJ_DIR)/rrsw_kg26.o  $(OBJ_DIR)/rrsw_ncpar.o  $(OBJ_DIR)/rrsw_wvn.o $(OBJ_DIR)/parkind.o		   $(OBJ_DIR)/rrsw_con.o  $(OBJ_DIR)/rrsw_kg19.o	$(OBJ_DIR)/rrsw_kg23.o  $(OBJ_DIR)/rrsw_kg27.o  $(OBJ_DIR)/rrsw_ref.o $(OBJ_DIR)/rrtmg_sw_cldprmc.o	$(OBJ_DIR)/rrtmg_sw_k_g.o	  $(OBJ_DIR)/rrtmg_sw_read_nc.o	$(OBJ_DIR)/rrtmg_sw_spcvmc.o  $(OBJ_DIR)/rrtmg_sw_vrtqdr.o $(OBJ_DIR)/mcica_random_numbers.o  $(OBJ_DIR)/rrtmg_sw_cldprop.o	$(OBJ_DIR)/rrtmg_sw_rad.o	  $(OBJ_DIR)/rrtmg_sw_reftra.o	$(OBJ_DIR)/rrtmg_sw_spcvrt.o $(OBJ_DIR)/mcica_subcol_gen_sw.o   $(OBJ_DIR)/rrtmg_sw_init.o	  $(OBJ_DIR)/rrtmg_sw_setcoef.o	$(OBJ_DIR)/rrtmg_sw_taumol.o 

all: $(TARGET) $(LIBRARY)
$(TARGET): $(TARGET_OBJ)
	$(FC) $(LINKFLAGS) -o $(TARGET) $(TARGET_OBJ)
	ln -fs $(TARGET) $(SOLIB)
$(LIBRARY): $(TARGET_OBJ)
	$(AR) $(ARFLAGS) $@ $?

install: $(LIBRARY)
	rm -f rrtmg_sw.mk.install
	@echo "# ESMF self-describing build dependency makefile fragment" > rrtmg_sw.mk.install
	@echo "# src location $(MACHINE): $pwd" >> rrtmg_sw.mk.install
	@echo  >> rrtmg_sw.mk.install
	@echo "ESMF_DEP_FRONT     = rrtmg_sw_cap_mod" >> rrtmg_sw.mk.install
	@echo "ESMF_DEP_INCPATH   = $(INSTALLDIR)" >> rrtmg_sw.mk.install
	@echo "ESMF_DEP_CMPL_OBJS = " >> rrtmg_sw.mk.install
	@echo "ESMF_DEP_LINK_OBJS = $(INSTALLDIR)/librrtmg_sw.a " >> rrtmg_sw.mk.install
	mkdir -p $(INSTALLDIR)
	mv *.mod $(INSTALLDIR) 
	cp -f librrtmg_sw.a mod/*.mod $(INSTALLDIR) 
	cp -f rrtmg_sw.mk.install $(INSTALLDIR)/rrtmg_sw.mk

clean:
	rm -rf *.o *.so *.mod obj/*.o mod/*.mod
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
# The rest is generated from makedepf90
# makedepf90 -I modules/:src/ -b obj -r '$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<' src/*
# they can be included through a depend file.

$(OBJ_DIR)/parkind.o : $(MODULES_DIR)/parkind.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/parrrsw.o : $(MODULES_DIR)/parrrsw.f90 $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_aer.o : $(MODULES_DIR)/rrsw_aer.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_cld.o : $(MODULES_DIR)/rrsw_cld.f90 $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_con.o : $(MODULES_DIR)/rrsw_con.f90 $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg16.o : $(MODULES_DIR)/rrsw_kg16.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg17.o : $(MODULES_DIR)/rrsw_kg17.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg18.o : $(MODULES_DIR)/rrsw_kg18.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg19.o : $(MODULES_DIR)/rrsw_kg19.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg20.o : $(MODULES_DIR)/rrsw_kg20.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg21.o : $(MODULES_DIR)/rrsw_kg21.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg22.o : $(MODULES_DIR)/rrsw_kg22.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg23.o : $(MODULES_DIR)/rrsw_kg23.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg24.o : $(MODULES_DIR)/rrsw_kg24.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg25.o : $(MODULES_DIR)/rrsw_kg25.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg26.o : $(MODULES_DIR)/rrsw_kg26.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg27.o : $(MODULES_DIR)/rrsw_kg27.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg28.o : $(MODULES_DIR)/rrsw_kg28.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_kg29.o : $(MODULES_DIR)/rrsw_kg29.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_ncpar.o : $(MODULES_DIR)/rrsw_ncpar.f90 $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_ref.o : $(MODULES_DIR)/rrsw_ref.f90 $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_tbl.o : $(MODULES_DIR)/rrsw_tbl.f90 $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_vsn.o : $(MODULES_DIR)/rrsw_vsn.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrsw_wvn.o : $(MODULES_DIR)/rrsw_wvn.f90 $(OBJ_DIR)/parrrsw.o $(OBJ_DIR)/parkind.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/mcica_random_numbers.o : $(SRC_DIR)/mcica_random_numbers.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/mcica_subcol_gen_sw.o : $(SRC_DIR)/mcica_subcol_gen_sw.f90 $(OBJ_DIR)/mcica_random_numbers.o $(OBJ_DIR)/mcica_random_numbers.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_cldprmc.o : $(SRC_DIR)/rrtmg_sw_cldprmc.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_cldprop.o : $(SRC_DIR)/rrtmg_sw_cldprop.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_init.o : $(SRC_DIR)/rrtmg_sw_init.f90 $(OBJ_DIR)/rrtmg_sw_setcoef.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_k_g.o : $(SRC_DIR)/rrtmg_sw_k_g.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_rad.o : $(SRC_DIR)/rrtmg_sw_rad.f90 $(OBJ_DIR)/rrtmg_sw_spcvmc.o $(OBJ_DIR)/rrtmg_sw_setcoef.o $(OBJ_DIR)/rrtmg_sw_cldprmc.o $(OBJ_DIR)/mcica_subcol_gen_sw.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_rad.nomcica.o : $(SRC_DIR)/rrtmg_sw_rad.nomcica.f90 $(OBJ_DIR)/rrtmg_sw_spcvrt.o $(OBJ_DIR)/rrtmg_sw_setcoef.o $(OBJ_DIR)/rrtmg_sw_cldprop.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_read_nc.o : $(SRC_DIR)/rrtmg_sw_read_nc.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_reftra.o : $(SRC_DIR)/rrtmg_sw_reftra.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_setcoef.o : $(SRC_DIR)/rrtmg_sw_setcoef.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_spcvmc.o : $(SRC_DIR)/rrtmg_sw_spcvmc.f90 $(OBJ_DIR)/rrtmg_sw_vrtqdr.o $(OBJ_DIR)/rrtmg_sw_taumol.o $(OBJ_DIR)/rrtmg_sw_reftra.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_spcvrt.o : $(SRC_DIR)/rrtmg_sw_spcvrt.f90 $(OBJ_DIR)/rrtmg_sw_vrtqdr.o $(OBJ_DIR)/rrtmg_sw_taumol.o $(OBJ_DIR)/rrtmg_sw_reftra.o 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_taumol.o : $(SRC_DIR)/rrtmg_sw_taumol.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
$(OBJ_DIR)/rrtmg_sw_vrtqdr.o : $(SRC_DIR)/rrtmg_sw_vrtqdr.f90 
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ -c $<
obj/rrtmg_cap.o : nuopc/rrtmg_cap.F90 
	$(FC) $(INCLUDE) $(ESMF_F90COMPILEOPTS) $(UTILINCS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -o $@ -c $<

