# Simple script to run all simulations at once

for i in ./*; do
	if [[ -d $i ]]; then
		cd $i
		echo $i
		for i in ./*; do
			if [[ $i == './EC' ]] || [[ $i == './noEC' ]]; then
				echo $i
				cd $i
				rm nohup.out
				nohup python -B sim_launcher.py &
				echo "Process ID: ${!}" >> aaa_info_dataset
				echo `cat aaa_info_dataset`
				echo ""
				cd ..
			fi
			cd ..
		done
	fi
done
