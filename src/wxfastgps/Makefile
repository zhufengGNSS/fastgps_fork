# source files
SRC  = fastgps_frame.cpp fastgps_app.cpp
OBJ  = $(SRC:.cpp=.o)
OUT  = ../../wxfastgps
LIBS = -lm -L../libfastgps -lfastgps `wx-config --libs`
CCC  = g++

INCLUDES = -I../libfastgps
CPPFLAGS = -g -O2 -msse2 -fomit-frame-pointer $(INCLUDES) `wx-config --cxxflags`
#CPPFLAGS = -g -O2 $(INCLUDES) `wx-config --cxxflags`

default: $(OUT)

.cpp.o:
	$(CCC) $(CPPFLAGS) -c $< -o $@
$(OUT): $(OBJ) ../libfastgps/libfastgps.a
	$(CCC) $(OBJ) $(LIBS) -o $(OUT)
clean:
	-rm $(OUT) *.o
	#rm *.o
	#cp fastgps_GUI ../
	#cp fastgps_GUI ../fastgps_GUI.app/Contents/MacOS

