### update by Mei-Shiue KUO
### date : 20181129

## query.py
it takes mfasta file and returns pssm matrix with length x 21 columns (20 amine acids + 1 gap)
	# agrv 0: query.py
	# agrv 1: mfasta file (input)
	# agrv 2: pssm file (output)

## query_ss2.py
it takes mfasta file and a ss2 file returned by psipred_runpsipred_single to returns pssm matrix with length x 24 columns (20 amine acids + 1 gap + 3 predictions of secondary structure : C, H and E)
* Must perform runpsipred_single to have secondary structure before the execution of query_ss2
	# agrv 0: query_ss2.py
	# agrv 1: mfasta file (input)
	# agrv 2: ss2 file (input)
	# agrv 3: pssm_ss2 file (output)

## template.py
it takes map file and a ss2 file returned by psipred_runpsipred_single to returns pssm matrix with length x 24 columns (20 amine acids + 1 gap + 3 predictions of secondary structure : C, H and E)
* Must perform runpsipred_single to have secondary structure before the execution of template_ss2
	# agrv 0: query_ss2.py
	# agrv 1: map file (input)
	# agrv 2: template_ss2 file (input)
	# agrv 3: pssm_ss2 file (output)

 
