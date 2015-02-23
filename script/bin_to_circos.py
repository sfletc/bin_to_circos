#Copyright Stephen Fletcher

import math
import csv
import os
import argparse

circos_loc = '~/circos-0.67-5/bin/circos'

"""
Take output from chromosome/bin files from /Users/steve/Software/sRNA_mapped_to_bin/sRNAmapped_to_bin_zero.py
and produce circos plot

Input example:
1. 24nt chromosome bin WT
2. 24nt chromosome bin mutant
3. 22nt chromosome bin WT
4. 22nt chromosome bin mutant
5. 21nt chromosome bin WT
6. 21nt chromosome bin mutant
7. No of chromosomes

***Chromsome number must always be final two numbers in bin file

Circos plot schema (outer to inner track):
1. Ideogram
2. 24nt mutant rel to WT - LFC (log10)
3. 22nt mutant rel to WT - LFC (log10)
4. 22nt mutant rel to WT - LFC (log10)
5. WT 24nt RPMR
6. WT 22nt RPMR
7. WT 21nt RPMR

"""

def load_bin(input_file):
    line_count = 0
    output_bin = []	
    
    with open(input_file, 'rU') as loaded_bin:

        for line in loaded_bin:
        	if line_count ==0 or line_count ==1:
        		pass
        	else:
        		clean_line = line.strip().split(',')
        		output_bin.append((clean_line[1],clean_line[2],clean_line[3]))
        	line_count +=1
    return output_bin[:-1]    	



def mod_circos_conf_template(circos_out, organism_species):
	line_count = 1
	if organism_species == 'Arabidopsis':
		with open('./template/circos_template_arab.conf', 'rU') as loaded_conf:
			with open('./template/{0}.conf'.format(circos_out), 'w') as conf_out:
				for line in loaded_conf:
					if line_count == 7:
						line = '	file* = {0}.png\n'.format(circos_out)
					elif line_count == 50:
						line = '		file    = data/{0}_24.csv\n'.format(circos_out)
					elif line_count == 78:
						line = '	    file    = data/{0}_22.csv\n'.format(circos_out)
					elif line_count == 107:
						line = '		file    = data/{0}_21.csv\n'.format(circos_out)
					elif line_count == 133:
						line = '		file    = data/{0}_RPMR_24.csv\n'.format(circos_out)
					conf_out.write(line)
					line_count+=1
	elif organism_species == 'tomato':
		with open('./template/circos_template_tomato.conf', 'rU') as loaded_conf:
			with open('./template/{0}.conf'.format(circos_out), 'w') as conf_out:
				for line in loaded_conf:
					if line_count == 7:
						line = '	file* = {0}.png\n'.format(circos_out)
					elif line_count == 50:
						line = '		file    = data/{0}_24.csv\n'.format(circos_out)
					elif line_count == 78:
						line = '	    file    = data/{0}_22.csv\n'.format(circos_out)
					elif line_count == 107:
						line = '		file    = data/{0}_21.csv\n'.format(circos_out)
					elif line_count == 133:
						line = '		file    = data/{0}_RPMR_24.csv\n'.format(circos_out)
					conf_out.write(line)
					line_count+=1


def circos_plot(nt24_a, nt24_b, nt22_a, nt22_b, nt21_a, nt21_b, circos_out, organism_species = 'Arabidopsis'):
	arab_chr = ['01','02','03','04','05']
	tomato_chr = ['01','02','03','04','05','06','07','08','09','10','11','12']
	circos_data_24=[]
	circos_data_22=[]
	circos_data_21=[]
	RPMR_data_24=[]
	circos_conf = './template/{0}.conf'.format(circos_out)
	
	if organism_species =='Arabidopsis': species_chr = arab_chr
	elif organism_species =='tomato': species_chr = tomato_chr

	chr_no = len(species_chr) #set for species

	for chromo in range(chr_no):
		bin_count_24 = 0
		bin_count_22 = 0
		bin_count_21 = 0

		bin_file_nt24_a = nt24_a+species_chr[chromo]+'.csv'
		bin_file_nt24_b = nt24_b+species_chr[chromo]+'.csv'
		bin_file_nt22_a = nt22_a+species_chr[chromo]+'.csv'
		bin_file_nt22_b = nt22_b+species_chr[chromo]+'.csv'			
		bin_file_nt21_a = nt21_a+species_chr[chromo]+'.csv'
		bin_file_nt21_b = nt21_b+species_chr[chromo]+'.csv'



		###24nt

		circos_nt24_a_chr = load_bin(bin_file_nt24_a)
		circos_nt24_b_chr = load_bin(bin_file_nt24_b)

		while bin_count_24 < len(circos_nt24_a_chr):
			try:
				
				circos_data_24.append(['chr'+species_chr[chromo], circos_nt24_a_chr[bin_count_24][0], circos_nt24_a_chr[bin_count_24][1],  "%.2f" %math.log10(float(circos_nt24_a_chr[bin_count_24][2])/float(circos_nt24_b_chr[bin_count_24][2]))])
				
			except:
				circos_data_24.append(['chr'+species_chr[chromo], circos_nt24_a_chr[bin_count_24][0], circos_nt24_a_chr[bin_count_24][1],  '0.00'])
			
			RPMR_data_24.append(['chr'+species_chr[chromo], circos_nt24_a_chr[bin_count_24][0], circos_nt24_a_chr[bin_count_24][1],  "%.2f" %float(circos_nt24_a_chr[bin_count_24][2])])
			bin_count_24+=1	
	

		with open('./data/'+circos_out+'_24.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter='\t')
			for row in circos_data_24:
				mycsv.writerow(row)
		with open('./data/'+circos_out+'_RPMR_24.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter='\t')
			for row in RPMR_data_24:
				mycsv.writerow(row)
		
		###22 nt

		circos_nt22_a_chr = load_bin(bin_file_nt22_a)
		circos_nt22_b_chr = load_bin(bin_file_nt22_b)

		while bin_count_22 < len(circos_nt22_a_chr):
			try:
				
				circos_data_22.append(['chr'+species_chr[chromo], circos_nt22_a_chr[bin_count_22][0], circos_nt22_a_chr[bin_count_22][1],  "%.2f" %math.log10(float(circos_nt22_a_chr[bin_count_22][2])/float(circos_nt22_b_chr[bin_count_22][2]))])
			except:
				circos_data_22.append(['chr'+species_chr[chromo], circos_nt22_a_chr[bin_count_22][0], circos_nt22_a_chr[bin_count_22][1],  '0.00'])
			bin_count_22+=1	
		with open('./data/'+circos_out+'_22.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter='\t')
			for row in circos_data_22:
				mycsv.writerow(row)

		###21 nt

		circos_nt21_a_chr = load_bin(bin_file_nt21_a)
		circos_nt21_b_chr = load_bin(bin_file_nt21_b)

		while bin_count_21 < len(circos_nt21_a_chr):
			try:
				
				circos_data_21.append(['chr'+species_chr[chromo], circos_nt21_a_chr[bin_count_21][0], circos_nt21_a_chr[bin_count_21][1],  "%.2f" %math.log10(float(circos_nt21_a_chr[bin_count_21][2])/float(circos_nt21_b_chr[bin_count_21][2]))])
			except:
				circos_data_21.append(['chr'+species_chr[chromo], circos_nt21_a_chr[bin_count_21][0], circos_nt21_a_chr[bin_count_21][1],  '0.00'])
			bin_count_21+=1
		
		with open('./data/'+circos_out+'_21.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter='\t')
			for row in circos_data_21:
				mycsv.writerow(row)		

	mod_circos_conf_template(circos_out, organism_species)
	circos_command = 'perl {0} -conf {1}'.format(circos_loc,circos_conf)
	os.system(circos_command)


def comline():
	"""
	Command line parser - add agruments here if needed
	"""
	parser = argparse.ArgumentParser("bin_to_circos.py")
	parser.add_argument('bin_file_a_24nt', type = str)
	parser.add_argument('bin_file_b_24nt', type = str)
	parser.add_argument('bin_file_a_22nt', type = str)
	parser.add_argument('bin_file_b_22nt', type = str)
	parser.add_argument('bin_file_a_21nt', type = str)
	parser.add_argument('bin_file_b_21nt', type = str)
	parser.add_argument('circos_out', type = str)
	parser.add_argument('-org','--organism_species', type=str)
	args = parser.parse_args()


	return args

def run_ana():
	"""
	Run the analysis.  GQ_cutoff and SNP filter defaults can be modified here
	"""

	a = comline()
	
	if a.organism_species == None:
		print 'No species chosen; default = Arabidopsis'
		circos_plot(a.bin_file_a_24nt,a.bin_file_b_24nt,a.bin_file_a_22nt,a.bin_file_b_22nt,a.bin_file_a_21nt,a.bin_file_b_21nt,a.circos_out)
	else:
		print 'Species = {0}'.format(a.organism_species)
		circos_plot(a.bin_file_a_24nt,a.bin_file_b_24nt,a.bin_file_a_22nt,a.bin_file_b_22nt,a.bin_file_a_21nt,a.bin_file_b_21nt,a.circos_out,a.organism_species)


run_ana()	

