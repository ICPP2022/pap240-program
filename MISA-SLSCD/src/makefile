# -------------------Begin-------------------------
#$Header:
#
#SHELL = /bin/sh
# -------------------End---------------------------

#Compiler and Compile options
FC = mpif90
FF = -g
CFLAGS = -O2

#Source code files, Object files and Target(Executable) file defeine
SRCS = $(wildcard src_*.f90)
OBJS = $(SRCS:.f90=.o)
TARGET = misascd

MODSRCS =   mod_structuretype.f90\
            mod_constants.f90\
            mod_globalvariables.f90\
            mod_randdp.f90\
            mod_update_reactions.f90
#MODSRCS =   $(wildcard mod_*.f90)   #make twice
MODOBJS = $(MODSRCS:.f90=.o)
MODS = $(MODSRCS:.f90=.mod)

#Generate Target file using Object files
$(TARGET): $(OBJS) $(MODOBJS)
	$(FC) $(CFLAGS) $^ -o $@

#Generate Object files using Source files
%.o: %.f90 $(MODS)
	$(FC) $(CFLAGS) -c $<

$(MODS) $(MODOBJS): $(MODSRCS)
	$(FC) $(CFLAGS) -c $^

# Clean Object files
clean:
	rm $(TARGET) $(OBJS) $(MODOBJS) $(MODS)