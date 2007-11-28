.SUFFIXES :
.SUFFIXES : .f90 .o

FC = g95
FFLAGS = -O3 -I/usr/local/include 
LD = f90
LDFLAGS = -flat_namespace -undefined suppress -L/usr/local/lib
LDLIBS = -lfftw3
MODEL = semilag

SRC = constant_module.f90 parameter_module.f90 \
	fft_module.f90 \
	glatwgt_module.f90 alf_module.f90 \
	legendre_transform_module.f90 init_module.f90 \
	upstream_module.f90 interpolate_module.f90 \
	polint_module.f90 cubicspline_module.f90 bicubic_module.f90 \
	$(MODEL)_module.f90 sphere_module.f90 io_module.f90 \
	main_$(MODEL).f90
OBJ = ${SRC:.f90=.o}
TARGET=adv

$(TARGET) : $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LDLIBS) -o $@

fft_module.o : constant_module.o
glatwgt_module.o : constant_module.o
alf_module.o : constant_module.o
parameter_module.o : constant_module.o
legendre_transform_module.o : glatwgt_module.o alf_module.o parameter_module.o constant_module.o
init_module.o : io_module.o sphere_module.o legendre_transform_module.o parameter_module.o constant_module.o
upstream_module.o : constant_module.o sphere_module.o glatwgt_module.o interpolate_module.o
interpolate_module.o : constant_module.o glatwgt_module.o sphere_module.o bicubic_module.o polint_module.o cubicspline_module.o
bicubic_module.o : constant_module.o
sphere_module.o : constant_module.o
io_module.o : constant_module.o
polint_module.o : constant_module.o
cubicspline_module.o : constant_module.o
eulerian_module.o : constant_module.o parameter_module.o io_module.o glatwgt_module.o legendre_transform_module.o init_module.o
semilag_module.o : constant_module.o parameter_module.o io_module.o legendre_transform_module.o init_module.o upstream_module.o
nisl_module.o : constant_module.o parameter_module.o io_module.o legendre_transform_module.o init_module.o upstream_module.o sphere_module.o
main_$(MODEL).o : init_module.o $(MODEL)_module.o parameter_module.o constant_module.o

clean :
	rm -f *.o *.mod *.dat *.log $(TARGET)

.f90.o :
	$(FC) $(FFLAGS) $< -c
