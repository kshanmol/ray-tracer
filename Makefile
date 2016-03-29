all:
	g++ -O3 -mssse3 -malign-double -c raytracer.cpp `libpng-config --cflags` -I $PREFIX/include
	g++  -O3 -mssse3 -malign-double -o ray_tracer -lpng12 raytracer.o `libpng-config --ldflags`
	rm -f raytracer.o


simple:
	g++ -o ray_tracer  -lpng12 raytracer.cpp

opt:
	g++ -O3 -o ray_tracer -lpng12 raytracer.cpp

clean:
	rm -f ray_tracer raytracer.o
