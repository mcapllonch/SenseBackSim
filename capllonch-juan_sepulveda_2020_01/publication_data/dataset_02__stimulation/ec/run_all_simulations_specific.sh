# Simple script to run specific simulations
# The folders listed in the conditional block in this script were left just as an example;
# the user should list his/her own choices

for i in ./*; do
	if [[ $i == './current_08000nA' ]] || [[ $i == './current_08500nA' ]] || [[ $i == './current_09000nA' ]] || [[ $i == './current_09500nA' ]] || [[ $i == './current_10000nA' ]]; then
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
