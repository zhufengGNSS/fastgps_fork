default:
	cd src/libfastgps && make
	cd src/fastgps && make
	cd src/wxfastgps && make

clean:
	cd src/libfastgps && make clean
	cd src/fastgps && make clean
	cd src/wxfastgps && make clean

.PHONY : clean
