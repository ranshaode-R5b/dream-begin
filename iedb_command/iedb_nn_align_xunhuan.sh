#! /bin/bash
path_fasta="/home/renshida/mhcii-data/IEDB_data_length"
path_code="/home/renshida/mhcii-data/code_bench"
path_result="/home/renshida/mhcii-data/code_result"
tool_consensus="/home/renshida/apps/IEDB_3.1.12"
	result_consensus="$path_result/iedb_nn_align"
	cd $tool_consensus ||{ echo "Failed to change directory to $tool_netmhccons"; exit 1; }
	com_file="$path_code/iedb_nn_align.sh"
	for file_path in "$path_fasta"/*.fa; do
		file_name2=$(basename "$file_path")
		file_name1=$(head -n 1 "$file_path")
		file_name=${file_name1:1}
		cleaned_name=$(echo "$file_name" | tr -d '[:punct:]')
		number_part="${file_name2##*_}"
		number="${number_part%.fa}"
		echo "python $tool_consensus/mhc_II_binding.py nn_align $file_name $file_path $number > $result_consensus/$cleaned_name.txt" >> "$com_file"
	done
	sed -i '1i #! /bin/bash' "$com_file"
