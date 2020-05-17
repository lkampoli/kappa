for dir in `ls -d [abimp]*`

do
   echo $dir
   cd $dir
   make
   ./runall.sh
   pwd
   cd ..
   pwd
done


