fasta_file = open("msa_kras_mutants.fasta")
output_file = open("kras_mutant_sequences.csv",'w')

output_file.write("mutant_name,sequence\n")

for line in fasta_file:
	if len(line)<10:
		output_file.write(line[1:-1] + ",")
	else:
		output_file.write(line)

output_file.close()
fasta_file.close()
