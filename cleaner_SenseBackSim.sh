#!/bin/bash
function explore {
		pwd
		for item in *; do 
			if [[ -d $item ]]; then
				# $item is a folder

				# Remove unwanted folders
				# if [[ $item == 'results' ]] || [[ $item == 'outputs' ]]; then
				# 	rm -rf $item

				# Create certain folder trees
				if [[ $item == 'data' ]]; then
					cd data
					mkdir results
					cd results
					mkdir images
					cd images
					touch .gitkeep
					cd ../../../
				else
					# If $item is not a folder we want to remove, we go inside and keep exploring its possible subfolders
					cd $item
					explore
				fi
			else
				# $item is not a folder, so it may be a file
				if [[ -f $item ]]; then
					# $item is a file
					if [[ $item == *'.xls'* ]]; then
						# rm *.xls
						rm $item
					fi

					# if [[ $item == *'.xml'* ]]; then
					# 	# rm *.xml
					# 	rm $item
					# fi

					if [[ $item == *'.pyc'* ]]; then
						rm $item
					fi

					if [[ $item == *'.blg'* ]]; then
						rm $item
					fi

					if [[ $item == *'.aux'* ]]; then
						rm $item
					fi

					if [[ $item == *'.log'* ]]; then
						rm $item
					fi

					if [[ $item == *'.mp4'* ]]; then
						rm $item
					fi

					if [[ $item == 'AXNODE.mod' ]]; then
						rm $item
					fi

					# if [[ $item == *'.msh'* ]] && [[ $item != *'axon1.msh'* ]]; then
					if [[ $item == *'.msh'* ]]; then
						if [[ $item != *'axon0.msh'* ]]; then
							# rm *.msh
							rm $item
						fi
					fi

				fi
			fi
		done
		cd ..
		pwd
}

# Call the function. This is all the main code does
explore