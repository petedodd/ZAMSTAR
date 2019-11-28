#!/bin/bash
# # NB - the next line is needed on Mac OS X because of the gcc version...
#   export PATH=/usr/local/bin:$PATH
gcc --version


echo option: $1

#BUILDING...
PAR=$1
if [ $PAR  =  3 ]
then
    echo "compiling parallel..."
    # gcc -Wall -c src/TBmain.cc -o  temp/TBbd.o -O3 -DPARA
    # # gcc -Wall -c -fopenmp src/IBM.cc -o temp/IBM.o -O3 -DPARA
    # gcc -Wall -c src/utils.cc -o temp/utils.o -O3
    # gcc -Wall -c src/parameters.cc -o temp/parameters.o -O3
    # gcc -Wall -c src/results.cc -o temp/results.o -O3
    # gcc -Wall -c src/person.cc -o temp/person.o -O3
    # gcc -Wall -c src/household.cc -o temp/household.o -O3
    # gcc -Wall -c -fopenmp src/population.cc -o temp/population.o -O3 -DPARA
    # gcc -Wall -c -fopenmp src/metapopulation.cc -o temp/metapopulation.o -O3 -DPARA
    # g++ -fopenmp temp/TBbd.o  temp/utils.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o temp/metapopulation.o -lgsl -lgslcblas -o TBbd

elif [ $PAR = 0 ]
then
    echo "no compilation..."
elif [ $PAR = 2 ]
then
    echo "compiling serial, no optimization..."
    # gcc -Wall -c src/TBmain.cc -o  temp/TBbd.o
    # # gcc -Wall -c src/IBM.cc -o temp/IBM.o
    # gcc -Wall -c src/utils.cc -o temp/utils.o
    # gcc -Wall -c src/parameters.cc -o temp/parameters.o
    # gcc -Wall -c src/results.cc -o temp/results.o
    # gcc -Wall -c src/person.cc -o temp/person.o
    # gcc -Wall -c src/household.cc -o temp/household.o
    # gcc -Wall -c src/population.cc -o temp/population.o
    # gcc -Wall -c src/metapopulation.cc -o temp/metapopulation.o
    # g++ temp/TBbd.o  temp/utils.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o temp/metapopulation.o -lgsl -lgslcblas -o TBbd

elif [ $PAR = 1 ]
then
    echo "compiling serial..."
    gcc -Wall -c src/TBmain.cc -o  temp/TB.o -O3
    gcc -Wall -c src/utilities.cc -o temp/utilities.o -O3
    gcc -Wall -c src/parameters.cc -o temp/parameters.o -O3
    gcc -Wall -c src/results.cc -o temp/results.o -O3
    gcc -Wall -c src/person.cc -o temp/person.o -O3
    gcc -Wall -c src/household.cc -o temp/household.o -O3
    gcc -Wall -c src/population.cc -o temp/population.o -O3
    g++ temp/TB.o temp/utilities.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o  -lgsl -lgslcblas  -o TB -O3
    # -std=c++11

elif [ $PAR = 5 ]
then
    echo "<<<VALGRIND!!>>>"
    # gcc -Wall -c src/TBmain.cc -o  temp/TBbd.o -O1 -g -fno-inline
    # # gcc -Wall -c src/IBM.cc -o temp/IBM.o -O1 -g -fno-inline
    # gcc -Wall -c src/utils.cc -o temp/utils.o -O1 -g -fno-inline
    # gcc -Wall -c src/parameters.cc -o temp/parameters.o -O1 -g -fno-inline
    # gcc -Wall -c src/results.cc -o temp/results.o -O1 -g -fno-inline
    # gcc -Wall -c src/person.cc -o temp/person.o -O1 -g -fno-inline
    # gcc -Wall -c src/household.cc -o temp/household.o -O1 -g -fno-inline
    # gcc -Wall -c src/population.cc -o temp/population.o -O1 -g -fno-inline
    # g++ temp/TBbd.o  temp/utils.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o  -lgsl -lgslcblas -o TBbd
    # valgrind -v  --log-file=vout.txt --leak-check=full --show-reachable=yes --dsymutil=yes ./TBbd data/

elif [ $PAR = 4 ]
then
    echo "DEBUG!"
    # gcc -Wall -c src/TBmain.cc -o  temp/TBbd.o -g
    # # gcc -Wall -c src/IBM.cc -o temp/IBM.o -g
    # gcc -Wall -c src/utils.cc -o temp/utils.o -g
    # gcc -Wall -c src/parameters.cc -o temp/parameters.o -g
    # gcc -Wall -c src/results.cc -o temp/results.o -g
    # gcc -Wall -c src/person.cc -o temp/person.o -g
    # gcc -Wall -c src/household.cc -o temp/household.o -g
    # gcc -Wall -c src/population.cc -o temp/population.o -g
    # g++ temp/TBbd.o  temp/utils.o temp/parameters.o temp/results.o temp/person.o temp/household.o temp/population.o  -lgsl -lgslcblas -o TBbd

else
    echo "wrong option..."
fi


if [ $2 = 1 ]
then
    echo "EXECUTING...."
    #EXECUTING...
    # data/use.sh S			# using this or that parm set
    # data/use.sh 413
    # time ./TB
    ./TB data/ 1
else
    echo "no execution ... "
fi


if [ $3 = 1 ]
then
      echo "PLOTTING..."
      #PLOTTING...

      # R plots
      R --vanilla --slave < src/hhplot.R #
      R --vanilla --slave < src/sazplot.R
      R --vanilla --slave < src/zmzplot.R
      R --vanilla --slave < src/snapshot.R  #


      # gnuplots
      cd src
      for i in *.gp;
      do
      	  gnuplot "$i"
      done
      cd ..

      # file conversion and tidy
      cd plots
      for i in *.ps;
      do
	  ps2pdf "$i" "${i/.ps}".pdf;
	  rm "$i";
      done

      #extra and bundle
      if [ -e Rplots.pdf ]
      then
	  rm Rplots.pdf 		# if R scripts used
      fi

      rm ALL.pdf
      pdftk *.pdf cat output ALL.pdf

       # make cd4 movie
       # cd ..
       # R --vanilla --slave < src/cd4movie.R

else
    echo "no plotting..."
fi
