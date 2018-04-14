# Check for proper number of command line args.

EXPECTED_ARGS=4
E_BADARGS=65

if [ $# -ne $EXPECTED_ARGS ]
then
echo "Usage: `basename $0` SIZE RADIUS PATH MAXCWS"
exit $E_BADARGS
fi

./voltool cornell $1 $3cornell$1.dat
./voltool grad $3cornell$1.dat $3cornell$1_32bit.dat
./voltool histogram $3cornell$1.dat $3cornell$1.hist
./voltool daosphere $2 $3cornell$1
./voltool vqpack $3cornell$1
./voltool vqtrain $3cornell$1 $4 5
./voltool vq $3cornell$1 $4 
./voltool vqunpack $3cornell$1 $4 
