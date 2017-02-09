all: pdistort

pdistort: pdistort.o
	gcc pdistort.o -o pdistort `pkg-config --libs libpng libperlin`

pdistort.o: pdistort.c
	gcc -c pdistort.c `pkg-config --cflags libpng libperlin` 
