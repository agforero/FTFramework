#cd $1
rm -f test.bats # reset it 
echo "# Generic Makefile test" >> test.bats
echo "@test \"generic Makefile test for $1\" {" >> test.bats
echo "  unset FC" >> test.bats
echo "  unset FFLAGS" >> test.bats
if [ $# -ne 1 ]; then # if arguments include more than just the directory
    echo "  export FC=$2" >> test.bats # make the new compiler the value for FC
    echo "  export FFLAGS=${@:3}" >> test.bats
else # if there was only the directory as an argument
    echo "  export FC=gfortran" >> test.bats
    echo "  export FFLAGS=" >> test.bats
fi
if [ -f env.txt ]; then
    echo cat env.txt >> test.bats
fi
echo "  make clean" >> test.bats
echo "  make" >> test.bats
echo "}" >> test.bats