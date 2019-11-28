#!/bin/bash
# # NB - the next line is needed on Mac OS X because of the gcc version...
#   export PATH=/usr/local/bin:$PATH
x86_64-w64-mingw32-gcc --version
i686-w64-mingw32-gcc --version

echo "compiling serial..."
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/TBmain.cc -o  temp/TB.o -O3
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/utilities.cc -o temp/utilities.o -O3
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/parameters.cc -o temp/parameters.o -O3
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/results.cc -o temp/results.o -O3
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/person.cc -o temp/person.o -O3
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/household.cc -o temp/household.o -O3
# x86_64-w64-mingw32-gcc -Wall -I/usr/local/include -c src/population.cc -o temp/population.o -O3
# x86_64-w64-mingw32-g++ temp/TB.o temp/utilities.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o  -lgsl -lgslcblas  -o TB.exe -O3


# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/TBmain.cc -o  temp/TB.o -O3
# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/utilities.cc -o temp/utilities.o -O3
# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/parameters.cc -o temp/parameters.o -O3
# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/results.cc -o temp/results.o -O3
# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/person.cc -o temp/person.o -O3
# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/household.cc -o temp/household.o -O3
# i686-w64-mingw32-gcc -Wall -I/usr/local/include -c src/population.cc -o temp/population.o -O3
# i686-w64-mingw32-g++ temp/TB.o temp/utilities.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o  -lgsl -lgslcblas  -o TB.exe -O3


/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/TBmain.cc -o  temp/TB.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/utilities.cc -o temp/utilities.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/parameters.cc -o temp/parameters.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/results.cc -o temp/results.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/person.cc -o temp/person.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/household.cc -o temp/household.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-gcc -Wall -I/usr/local/include -c src/population.cc -o temp/population.o -O3
/home/pjd/Downloads/mxe/usr/bin/i686-w64-mingw32.static-g++ temp/TB.o temp/utilities.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o  -lgsl -lgslcblas  -o TB.exe -O3
