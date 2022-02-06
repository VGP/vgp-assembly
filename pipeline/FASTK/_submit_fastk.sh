#!/bin/bash

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF
	
	Usage: '$0 -l filelist -o outdir -p partition -c cpu [-a 10x]'
	
	-a use 10x to trim barcodes 

EOF

exit 0

fi

while getopts "l:o:p:c:a:k:" opt; do

	case $opt in
		l)
			filelist="$OPTARG"
			;;
		o)
			outdir="$OPTARG"
			;;
        p)
        	partition="$OPTARG"
            ;;
        c)
        	cpu="$OPTARG"
            ;;
        a)
        	optarg="$OPTARG"
            ;;
	k)
          	klen="$OPTARG"
            ;;
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac

done

file_N=$(wc -l $filelist | cut -f1 -d' ')

mkdir -p $outdir logs

rm -f fastk_count.jid

log=logs/$outdir.count.%A_%a.log
script=fastk_array.sh

echo "\
sbatch -p $partition -c $cpu --array=1-${file_N} --job-name=$outdir.count --error=$log --output=$log $script $filelist $outdir $cpu $klen $optarg | tail -n 1 | cut -d' ' -f4 >> fastk_count.jid"
sbatch -p $partition -c $cpu --array=1-${file_N} --job-name=$outdir.count --error=$log --output=$log $script $filelist $outdir $cpu $klen $optarg | tail -n 1 | cut -d' ' -f4 >> fastk_count.jid

WAIT="afterok:"

job_nums=`wc -l fastk_count.jid | awk '{print $1}'`

if [[ $job_nums -eq 1 ]]; then
   jid=`cat fastk_count.jid`
   WAIT=$WAIT$jid
else
    for jid in $(cat fastk_count.jid)
    do
        WAIT=$WAIT$jid","
    done
fi

log=logs/$outdir.sum.%A.log
script=fastk_union.sh

echo "\
sbatch -p $partition -c $cpu --job-name=$outdir.sum --error=$log --output=$log --dependency=$WAIT $script $filelist $outdir $cpu"
sbatch -p $partition -c $cpu --job-name=$outdir.sum --error=$log --output=$log --dependency=$WAIT $script $filelist $outdir $cpu
