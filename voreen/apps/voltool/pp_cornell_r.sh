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
./voltool filtergrad $3cornell$1_32bit.dat $3cornell$1_32bit-tf.dat
./voltool histogram $3cornell$1.dat $3cornell$1.hist
./voltool r_daosphere $2 $3cornell$1
./voltool vqpack $3cornell$1
./voltool r_vqtrain $3cornell$1 $4 
./voltool r_vq $3cornell$1 $4 
./voltool vqunpack $3cornell$1 $4 
./voltool dao32216 $3cornell$1_dao.dat 0
cp $3cornell$1_dao.cb $3cornell$1_dao16-0.cb
mv $3 /share/voreen/voreen-dao