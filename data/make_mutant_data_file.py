original_data_file = open("mutants_data.txt")
reformatted_data_file = open("kras_mutants_data.csv",'w')

reformatted_data_file.write("mutant_name,percent_folded\n")

for line in original_data_file:
	if len(line)>1:
		entries = line.split()
		percent_folded = 1.0 - float(entries[1])/100.
		reformatted_data_file.write(entries[0] + "," + str(percent_folded) +"\n")

reformatted_data_file.close()
original_data_file.close()
