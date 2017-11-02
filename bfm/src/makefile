#include ${PETSC_DIR}/conf/variables
FC        = gfortran

#    IFORT
#FFLAGS    = -traceback -O1 -g -assume byterecl -convert big_endian -fp-model precise
#FFLAGS    = -traceback -check all -g -warn interfaces,nouncalled -assume byterecl -convert big_endian
#FFLAGS    = -traceback -O1 -g -cpp ${PETSC_FC_INCLUDES} -assume byterecl -convert big_endian
#FFLAGS    = -traceback -check all -g -warn interfaces,nouncalled -cpp ${PETSC_FC_INCLUDES} -assume byterecl -convert big_endian -fp-stack-check -fp-model precise -fpe0

#    GFORTRAN
#FFLAGS    = -g -ffree-line-length-0 -fbacktrace -Wall -cpp ${PETSC_FC_INCLUDES}
#FFLAGS    = -g -O3 -ffree-line-length-0 -fbacktrace -Wall -cpp ${PETSC_FC_INCLUDES}
#FFLAGS    = -Wall -O0 -g -ggdb -ffree-line-length-none -fcheck=all -fsanitize=thread -fbacktrace -fstack-protector-all -ffpe-trap=invalid,zero,overflow
#FFLAGS    = -Wall -O0 -g -ggdb -ffree-line-length-none -fcheck=all -fsanitize=address -fsanitize=undefined   -fbacktrace -fstack-protector-all -ffpe-trap=invalid,zero,overflow
 FFLAGS    = -Wall -O3 -g -ggdb -ffree-line-length-none -ffpe-trap=invalid,zero,overflow


BINDIR = exe

SCRATCH_DIR = Build
EXE_DIR     = exe

RULEBOOK   = RuleBook

SCRIPT_RULE   := ./scripts/BFM_CreateRelevantObjectList
SFMAKEDEPEND  := ./scripts/sfmakedepend
SimpleDEPEND  := ./scripts/WW3_GetSimpleDepend
SCRIPT_OBJECT := ./scripts/BFM_ObjectFiles
RULE_BINARY   := ./scripts/BFM_GetEXEC
          AR := ar
     ARFLAGS := r
      RANLIB := ranlib
        TEST := test
       MKDIR := mkdir -p
MDEPFLAGS += --silent --moddir $(SCRATCH_DIR)

   FFLAGS += -cpp -I$(NETCDF_INCDIR) -I../Includes
#     LIBS := -L$(NETCDF_LIBDIR) -lnetcdf ${PETSC_LIB}
#     LIBS := $(NETCDF_FORTRAN_LINK) ${PETSC_LIB}
     LIBS := $(NETCDF_FORTRAN_LINK)

.PHONY: clean
clean:
	$(QUIET)rm -rf Build/* 
	$(QUIET)rm -rf exe/* 

SRC1 = $(notdir $(wildcard BFM/*/*.F90))
OBJ1 = $(addprefix $(SCRATCH_DIR)/, $(patsubst %.F90,%.o, $(notdir $(SRC1))))

SRC2 = $(notdir $(wildcard nemo/*.F90))
OBJ2 = $(addprefix $(SCRATCH_DIR)/, $(patsubst %.F90,%.o, $(notdir $(SRC2))))

SRC3 = $(notdir $(wildcard share/*.F90))
OBJ3 = $(addprefix $(SCRATCH_DIR)/, $(patsubst %.F90,%.o, $(notdir $(SRC3))))

SRC4 = $(notdir $(wildcard standalone/*.F90))
OBJ4 = $(addprefix $(SCRATCH_DIR)/, $(patsubst %.F90,%.o, $(notdir $(SRC4))))


OBJtot := $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)

-include $(SCRATCH_DIR)/OBJsrc
-include $(SCRATCH_DIR)/MakeDepend
-include $(SCRATCH_DIR)/MakeSimpleObject
-include $(SCRATCH_DIR)/MakeSimpleBINARIES


MDEPFLAGS := --file=- --silent --moddir $(SCRATCH_DIR) --objdir=$(SCRATCH_DIR)

$(SCRATCH_DIR)/OBJsrc: $(SCRATCH_DIR)
	$(SCRIPT_RULE) BFM/*/*.F90 nemo/*.F90 share/*.F90 standalone/*.F90 > $(SCRATCH_DIR)/OBJsrc

$(SCRATCH_DIR)/MakeDepend: $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(SWITCHFILE) BFM/*/*.F90 nemo/*.F90 share/*.F90 standalone/*.F90 srcEXEC/*.F90 > $(SCRATCH_DIR)/MakeDepend

$(SCRATCH_DIR)/MakeSimpleObject: $(SCRATCH_DIR)
	$(SCRIPT_OBJECT) BFM/*/*.F90 nemo/*.F90 share/*.F90 standalone/*.F90 srcEXEC/*.F90 > $(SCRATCH_DIR)/MakeSimpleObject

$(SCRATCH_DIR)/MakeSimpleBINARIES: $(SCRATCH_DIR) $(BINDIR)
	$(RULE_BINARY) srcEXEC/*.F90 > $(SCRATCH_DIR)/MakeSimpleBINARIES

$(SCRATCH_DIR)/libtotal.a: $(SCRATCH_DIR)/OBJsrc $(SCRATCH_DIR)/MakeDepend $(SCRATCH_DIR)/MakeSimpleBINARIES $(OBJsrc)
	$(AR) $(ARFLAGS) $(SCRATCH_DIR)/libtotal.a $(OBJsrc)
	$(RANLIB) $(SCRATCH_DIR)/libtotal.a

$(SCRATCH_DIR):
	$(shell $(TEST) -d $(SCRATCH_DIR) || $(MKDIR) $(SCRATCH_DIR) )

$(BINDIR):
	$(shell $(TEST) -d $(BINDIR) || $(MKDIR) $(BINDIR) )
# DO NOT DELETE THIS LINE - used by make depend
