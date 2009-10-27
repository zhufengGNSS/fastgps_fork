SRC  = main.cpp
OBJ  = $(SRC:.cpp=.o)
OUT  = ../../fastgps
CCC  =  g++
LIBS = -lm -L../libfastgps -lfastgps

# flags for Intel Core Duo
INCLUDES = -I../libfastgps
CPPFLAGS = -Wall -Werror -g -O2 -msse2 -fomit-frame-pointer $(INCLUDES)
#CPPFLAGS= -g -O2 $(INCLUDES)

default: .depend $(OUT)

$(OUT): $(OBJ) $(COBJ) ../libfastgps/libfastgps.a
	$(CCC) $(OBJ) $(LIBS) -o $(OUT)
.cpp.o:
	$(CCC) $(CPPFLAGS) -c $< -o $@
.depend: 
	$(CC) $(INCLUDES) -MM $(SRC) >.depend
clean:
	-rm $(OBJ) $(OUT)
	-rm .depend

ifeq (,$(findstring clean,$(MAKECMDGOALS)))
-include .depend
endif

.PHONY : clean

