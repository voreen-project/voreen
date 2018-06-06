# Check for proper number of command line args.

EXPECTED_ARGS=4
E_BADARGS=65

if [ $# -ne $EXPECTED_ARGS ]
then
echo "Usage: `basename $0` SIZE RADIUS PATH MAXCWS"
exit $E_BADARGS
fi

./voltool synth $1 $3synth$1
./voltool grad $3synth$1.dat $3synth$1_32bit
./voltool histogram $3synth$1.dat $3synth$1.hist
./voltool daosphere $2 $3synth$1
./voltool vqpack $3synth$1
./voltool vqtrain $3synth$1 $4 
./voltool vq $3synth$1 $4 
./voltool vqunpack $3synth$1 $4 
