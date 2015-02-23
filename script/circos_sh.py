#Copyright Stephen Fletcher
###run from root directory

"""
Input format - tab sepated:

mutant WT
"""

def circos_sh(in_list, out_file, species = 'Arabidopsis'):

	sample_names=[]
	with open(in_list, 'rU') as loaded_list:


		for line in loaded_list:

			sample_names.append((line.strip().split('\t')))
	print sample_names
	with open(out_file, 'w') as write_file:


		for comparison in sample_names:
			conf_out = comparison[0]+'_'+comparison[1]

	    	
			line =  'python ./script/bin_to_circos.py ./bin_files/{0}  ./bin_files/{1}  ./bin_files/{2}  ./bin_files/{3} ./bin_files/{4}  ./bin_files/{5} {6} -org {7} \n'.format(comparison[0]+'_24_bin_chr', comparison[1]+'_24_bin_chr', comparison[0]+'_22_bin_chr', comparison[1]+'_22_bin_chr', comparison[0]+'_21_bin_chr', comparison[1]+'_21_bin_chr', conf_out, species)
			write_file.write(line)



def dm_sh(seq_file_list, out_file, species = 'Arabidopsis', window = 100000):
    arab_chr = ['01','02','03','04','05']
    ref_arab_chr = ['TAIR10_chr01.fas', 'TAIR10_chr02.fas', 'TAIR10_chr03.fas', 'TAIR10_chr04.fas','TAIR10_chr05.fas']
    nt = ['21', '22', '24']


    seq_names=[]

    with open(seq_file_list, 'rU') as seq_list:

        for line in seq_list:
        	seq_names.append(line.strip())
    print seq_names
    with open(out_file, 'w') as write_file:


	    for seq_name in seq_names:
	    	chromo = 0
	    	while chromo < len(arab_chr):
	    		for sRNA_len in nt:
	    			out_file_dm = './bin_files/{0}_{1}_{2}.csv'.format(seq_name,sRNA_len,arab_chr[chromo])
	    			out_file_bin_dm = './bin_files/{0}_{1}_bin_chr{2}.csv'.format(seq_name,sRNA_len,arab_chr[chromo])
	    			line_1 = 'python ~/Software/sRNAmap_cur/sRNAmapper.py -ana dm -s1 {0} -r {1} -nt1 {2} -o {3} -nosplit\n'.format('./seq/{0}'.format(seq_name), './ref/{0}'.format(ref_arab_chr[chromo]), sRNA_len, out_file_dm)
	    			line_2 = 'python ~/Software/sRNA_mapped_to_bin/sRNAmapped_to_bin_zero.py -i {0} -o {1}  -b {2}\n'.format(out_file_dm, out_file_bin_dm, window)
	    			write_file.write(line_1)
	    			write_file.write(line_2)

	    		chromo+=1

    	
	    	



dm_sh('in_seq.txt', 'dm_bin.sh')
circos_sh('comparison.txt','circos.sh')

