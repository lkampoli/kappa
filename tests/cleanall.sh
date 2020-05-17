for dir in `ls -d [abimp]*`

do
   echo $dir
   cd $dir
   make clean
   pwd
   cd ..
   pwd
done


