################################ netmhcpan_ba-4.1
    #https://downloads.iedb.org/tools/mhci/LATEST/README
    tool_netmhcpan_ba="/data1/renshida/data/mhc_benchmark/tools/mhc_i"
    result_netmhcpan_ba="$path_result/netmhcpan_ba"
    #run configure
    cd "$tool_netmhcpan_ba" || { echo "Failed to change directory to $tool_netmhcpan_ba"; exit 1; }
    ./configure || { echo "Configuration failed"; exit 1; }
    #produce command lines
    rm "$path_code/netmhcpan_ba_command.sh"
    while IFS= read -r file_path; do
        file_name=$(basename "$file_path")
        hla_part="${file_name%%_*}"
        number_part="${file_name##*_}"
        number="${number_part%.fasta}" #remove appendix
        echo "$tool_netmhcpan_ba/src/predict_binding.py netmhcpan_ba \"$hla_part\" $number $file_path > $result_netmhcpan_ba/${hla_part}_$number.txt" >> "$path_code/netmhcpan_ba_command.sh"
    done < "$filelist_mhci_fasta"
    #give #!/bin/bash
    com_file="$path_code/netmhcpan_ba_command.sh"
    sed -i '1i #!/bin/bash' "$com_file"
    #add x right
    chmod +x "$com_file"
    # /home/renshida/biocode/mhc_benchmark/netmhcpan_ba_command.sh > /home/renshida/biocode/mhc_benchmark/mhcflurry.log 2>&1 &
    # jobs
    # disown -h %1 #(%1 means %[num],num is job number , can change)
    #successfully working if python process is running
################################ netmhcpan_el-4.1
    #https://downloads.iedb.org/tools/mhci/LATEST/README
    tool_netmhcpan_el="/data1/renshida/data/mhc_benchmark/tools/mhc_i"
    result_netmhcpan_el="$path_result/netmhcpan_el"
    #run configure
    cd "$tool_netmhcpan_el" || { echo "Failed to change directory to $tool_netmhcpan_el"; exit 1; }
    ./configure || { echo "Configuration failed"; exit 1; }
    #produce command lines
    rm "$path_code/netmhcpan_el_command.sh"
    while IFS= read -r file_path; do
        file_name=$(basename "$file_path")
        hla_part="${file_name%%_*}"
        number_part="${file_name##*_}"
        number="${number_part%.fasta}" #remove appendix
        echo "$tool_netmhcpan_el/src/predict_binding.py netmhcpan_el \"$hla_part\" $number $file_path > $result_netmhcpan_el/${hla_part}_$number.txt" >> "$path_code/netmhcpan_el_command.sh"
    done < "$filelist_mhci_fasta"
    #give #!/bin/bash
    com_file="$path_code/netmhcpan_el_command.sh"
    sed -i '1i #!/bin/bash' "$com_file"
    #add x right
    chmod +x "$com_file"
    # /home/renshida/biocode/mhc_benchmark/netmhcpan_el_command.sh > /home/renshida/biocode/mhc_benchmark/mhcflurry.log 2>&1 &
    # jobs
    # disown -h %1 #(%1 means %[num],num is job number , can change)
    #successfully working if python process is running