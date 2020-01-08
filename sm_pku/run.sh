#!/bin/bash
#SBATCH -A hubsh
#SBATCH -o out%j_SM.dat
#SBATCH --partition=compute_B
#SBATCH --qos=low
#SBATCH -J SM
#SBATCH --get-user-env
#SBATCH --nodes=10
#SBATCH --nodelist=node[01,04,05,06,08,09,10,11,12,13]
##SBATCH --nodelist=node[18,19]
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=end
#SBATCH --time=480:00:00

srun hostname -s | sort -u > slurm.hosts

cd $PWD

path="./operators_input/PKU_format/"
outfile="sd-shell_n3lo_im3b_-0.2_0.71_1.05_N12_hw24_O16A16_HF"
 sp="sd-shell_n3lo_im3b_-0.2_0.71_1.05_N12_hw24_O16A16_HF_pku.sp"
int="sd-shell_n3lo_im3b_-0.2_0.71_1.05_N12_hw24_O16A16_HF_pku.int"
nucleus=(O F Ne Na Mg Al Si P S)
coreA=16
EM=E2
Jm2=0
Par=0
OBorTB="1"

for vZ in 4 5 6 7 8 
#for vZ in 0 2 4 6 8
#for vZ in 1 3 5 7
#for vZ in 0 1 2 3 4 5 6 7 8
do
#for vN in 7
#for vN in 1 3 5 7 9 11
#for vN in 0 2 4 6 8 10 12
for vN in 0 1 2 3 4 5 6 7 8
do

name=${nucleus[${vZ}]}
#name=S
vA=`expr ${vN} + ${vZ}`
A=`expr ${coreA} + ${vA}`

JM2=`echo "scale=0; $A%2"|bc`
if [ ${JM2} -lt 1 ]
then
JM2=${Jm2}
fi
#if [ ${Jm2} -lt 2 ]
#then
#JM2=`echo "scale=0; $A%2"|bc`
#else
#JM2=${Jm2}
#fi
PAR=`echo "scale=0; ${Par}*(${JM2}%2)"|bc`
#PAR=${Par}
echo ${outfile}
echo ${name}$A"   "${EM}"    JM2:"${JM2}"    PAR:"${PAR}

echo ${path}${sp} > input4.txt
echo ${path}${int} >> input4.txt
if [ ${JM2} -eq 2 ]
then
echo "./results/IMSRG/"${name}${A}_JM2_${EM}_${outfile}".dat" >> input4.txt
else
echo "./results/IMSRG/"${name}${A}_${EM}_${outfile}".dat" >> input4.txt
fi
echo ${vN}>> input4.txt
echo ${vZ}>> input4.txt
echo "10000  0  0">> input4.txt
echo ${JM2}>> input4.txt
echo ${PAR}>> input4.txt
echo "10">> input4.txt
echo "0">> input4.txt
if [ ${OBorTB} == "0" ] 
then
echo "0 "${path}${outfile}"_"${EM}"_1b.op" >> input4.txt
else
echo "1 "${path}${outfile}"_"${EM}"_1b.op  "${path}${outfile}"_"${EM}"_2b.op" >> input4.txt
fi

echo >> input4.txt
echo >> input4.txt
echo >> input4.txt
echo >> input4.txt
echo "#############################################
spfile          #single particle file
vfile           #interaction file
outputfile      #output file
vN              #valence neutron number
vZ              #valence proton number
Emax, restriction for n,p   #restriction on configuration energy, total scatering state
MM              #all configurations have angular momentum projection MM=2M
                 #for even A, MM=0; for odd A, MM=1
Par             #parity, 0 for +, 1 for -
num of states   #number of states wanted
printOccuNums   #whether print occupation number, 1 for yes, 0 for no
calOper   1BOperfile  2BOperfile      #whether cal. other operators, 0 for no, 1 for yes.
                                      #if yes, the following file will be used
" >> input4.txt

#./test.exe < input4.txt
mpiexec -n 10 -machinefile slurm.hosts ./test.exe < input4.txt
#/opt/gopenmpi/bin/mpiexec -n 1 -machinefile slurm.hosts ./test.exe < input4.txt

done
done
