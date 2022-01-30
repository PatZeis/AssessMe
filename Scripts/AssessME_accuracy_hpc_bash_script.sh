#! /bin/bash
#SBATCH -p bioinfo
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --export=ALL
#SBATCH --mail-user=zeis@ie-freiburg.mpg.de
#SBATCH --mail-type=ALL


assess=''
object=''
assess_ind=1
iterator=1
cross=50
samp=70
tree=200
logreg=false
rawdata=false

print_usage() {
  printf "Usage: ..."
}


while getopts 'a:j:t:i:c:p:e:r:g' flag; do  ### : behind flag = flag can have an argument, if not flag = T, no flag = F
    case "${flag}" in 
        a) assess="${OPTARG}" ;;
	j) object="${OPTARG}" ;;
	t) iterator="${OPTARG}" ;;
	i) assess_ind="${OPTARG}" ;;
	c) cross="${OPTARG}" ;;
        p) samp="${OPTARG}" ;;
	e) tree="${OPTARG}" ;;
	g) logreg='true' ;;
	r) rawdata="${OPTARG}" ;;
	*) print_usage
 	   exit 1 ;;
    esac
done
if [ "$logreg" = false ]
then
    Rscript --vanilla /home/zeis/AssessME_accuracy_hpc_script.R -a $assess -j $object -i $assess_ind -t $iterator -c $cross -p $samp -e $tree -r $rawdata
else 
    Rscript --vanilla /home/zeis/AssessME_accuracy_hpc_script.R -a $assess -j $object -i $assess_ind -t $iterator -c $cross -p $samp -e $tree -g $logreg -r $rawdata
fi
