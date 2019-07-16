sudo apt-get install r-base-core
gcc   -g -c bamcat.c 
g++  -g -c samToErrorRate.c
g++   -o samToErrorRate -lz -lpthread bamcat.o samToErrorRate.o
