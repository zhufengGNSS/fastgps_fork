SRC = fastgps.cpp \
      code_table.cpp \
      acquire.cpp \
      correlator.cpp \
      tracking.cpp \
      nav.cpp \
      lin_alg.cpp \
      ephemerides.cpp \
      datetime.cpp \
      intrpsp3c.cpp \
      gnss_utils.cpp

CSRC = kiss_fft.c kiss_fftr.c

OBJ  = $(SRC:.cpp=.o)
COBJ = $(CSRC:.c=.o)

OUT = libfastgps.a

CC   =  gcc
CCC  =  g++
LIBS = -lm 

# flags for Intel Core 2
CFLAGS   = -g -O2 -march=prescott -msse2 -fomit-frame-pointer -Wall
CPPFLAGS = -g -O2 -march=prescott -msse2 -fomit-frame-pointer -Wall

default: .depend $(OUT)

$(OUT): $(OBJ) $(COBJ)
	ar rcs $(OUT) $(OBJ) $(COBJ)
.cpp.o:
	$(CCC) $(CPPFLAGS) -c $< -o $@
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.depend: 
	$(CC) $(INCLUDES) -MM $(SRC) $(CSRC) >.depend
	
clean:
	-rm $(OBJ) $(COBJ) $(OUT)
	-rm .depend

ifeq (,$(findstring clean,$(MAKECMDGOALS)))
-include .depend
endif

.PHONY : clean
