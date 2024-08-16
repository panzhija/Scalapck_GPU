#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS         = -DAdd_


#
#  The fortran and C compilers, loaders, and their flags
#

FC            = mpiifort
CC            = mpiicc 
NOOPT         = -O0
FCFLAGS       = -O3
CCFLAGS       = -O3
FCLOADER      = $(FC)
CCLOADER      = $(CC)
FCLOADFLAGS   = $(FCFLAGS)
CCLOADFLAGS   = $(CCFLAGS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#

ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = /usr/bin/ranlib

#
#  The name of the ScaLAPACK library to be created
#

SCALAPACKLIB  = libhipscalapack.a

#
#  BLAS, LAPACK (and possibly other) libraries needed for linking test programs
#

BLASLIB       = -L/usr/lib64/ -lblas
LAPACKLIB     = -L/usr/lib64/ -llapack

HIPDIR       = /public/software/compiler/dtk/21.04
HIPINC       = -I$(HIPDIR)/include -D__HIP_PLATFORM_HCC__
HIPLIBS      = -L/public/software/compiler/dtk/21.04/hipblas/lib/ -lhipblas
HIPLIBS      += -L/public/software/compiler/dtk/21.04/lib/ -lamdhip64
CCFLAGS      += $(HIPINC)

#HIPDIR       = /public/software/compiler/dtk/dtk-22.10
#HIPINC       = -I$(HIPDIR)/include
#HIPLIBS      = -L/public/software/compiler/dtk/dtk-22.10/hipblas/lib/ -lhipblas
#HIPLIBS      += -L/public/software/compiler/dtk/dtk-22.10/lib/ -lamdhip64
#CCFLAGS      += $(HIPINC)

MAGMADIR     = /public/software/mathlib/magma/magma-2.7.1
MAGMAINC     = -I$(MAGMADIR)/include
MAGMALIBS    = -L/public/software/mathlib/magma/magma-2.7.1/lib -lmagma
CCFLAGS      += $(MAGMAINC)

#MAGMADIR     = /public/home/aicao/panzhijie/test/magma
#MAGMAINC     = -I$(MAGMADIR)/include
#MAGMALIBS    = -L$(MAGMADIR)/lib -lmagma
#CCFLAGS      += $(MAGMAINC)

LIBS          = $(LAPACKLIB) $(BLASLIB) $(HIPLIBS) $(MAGMALIBS)
#LIBS          = $(LAPACKLIB) $(BLASLIB) $(HIPLIBS)