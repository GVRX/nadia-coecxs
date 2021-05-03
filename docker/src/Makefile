export LD_RUN_PATH=:/usr/include/hdf4/lib/:/home/ics/sw/NADIA/lib

C_SOURCES =$(wildcard *.c)


EXEC=$(C_SOURCES:.c=.exe)

all: $(EXEC)

%.exe: %.c
	g++ -I/usr/local/hdf4/include -O3 -fPIC -ffast-math  -DDOUBLE_PRECISION -I/opt/nadia-coecxs/NADIA/include $< -o $@ \
	  -L/opt/nadia-coecxs/NADIA/lib -lNADIAd \
	  -L/usr/local/hdf4/lib/  -llapack -ltiff -lmfhdf -ldf -lz -ljpeg -lpthread -lm -lfftw3 -lstdc++ 


clean:                                                     
	rm -f *.exe* *~
