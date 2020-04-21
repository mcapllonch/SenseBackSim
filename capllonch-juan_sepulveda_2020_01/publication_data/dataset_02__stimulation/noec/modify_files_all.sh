# Simple script to modify some files for all simulations
# Use this script in case you want to change a file in all or some of the folders (modify it at will for this)
# Note that the origin files presented in this script won't necessarily work; they are presented just as an example. The user must type in the paths to his/her desired origin and destination files, as well as the commands he/she sees fit.

for i in ./*; do
	if [[ -d $i ]]; then
		echo $i
		cd $i
		for i in ./*; do
			if [[ -d $i ]]; then
				if [[ $i == './EC' ]] || [[ $i == './noEC' ]]; then
					echo $i
					cd $i
					rm ./*py
					rm ./settings/tissues_JA1.json
					cp ../../../../../codeversion/revisions_round_2/*py .
					cp ../../../../../codeversion/revisions_round_2/settings/tissues.json ./settings/
#					echo "UPDATE: 2ND ROUND OF REVISIONS:" >> aaa_info_dataset
#					echo "This version of the code does NOT adjust rho_endo_longt of the fascicles as I did originally to make rho_endo_longt_eq match the nominal value 175 Ohm*cm. Instead, rho_endo is a isotropic resistivity with value 1211 Ohm*cm." >> aaa_info_dataset
#					echo "The only adjustment is to xraxial, by lumping  endoneurium and epineurium for those fibers overlapping with the epineurium. But this is something I did after the first round of revisions." >> aaa_info_dataset
#					echo ""
#					echo `cat aaa_info_dataset`
					cd ..
				fi
			fi
		done
		cd ..
	fi
done
