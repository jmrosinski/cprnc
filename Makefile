FC = ifort
FFLAGS = -g -O3 -I$(NETCDF)/include
LDFLAGS = -L$(NETCDF)/lib -lnetcdf
OBJS = utils.o stats.o cprnc.o

%.o: %.F90
	$(FC) -c $(FFLAGS) $<

cprnc: $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) *.mod cprnc

utils.o:
stats.o: utils.o
cprnc.o: utils.o
