CC=g++
CFLAGS = -std=c++14

test: Exec
	./Exec

Exec: test.o fft.o
	$(CC) $(CFLAGS) -o Exec test.o fft.o

test.o: test.cpp fft.h
	$(CC) $(CFLAGS) -c test.cpp

fft.o: fft.cpp fft.h
	$(CC) $(CFLAGS) -c fft.cpp

clean:
	rm -rf *.o Exec
