# update by Mei-Shiue KUO
# date : 20181125

### query.py
it takes mfasta file and returns pssm matrix with length x 21 columns (20 amine acids + 1 gap)
	# agrv 0: query.py
	# agrv 1: mfasta file (input)
	# agrv 2: pssm file (output)

### query_ss2.py
it takes mfasta file and a ss2 file returned by psipred_runpsipred_single to returns pssm matrix with length x 24 columns (20 amine acids + 1 gap + 3 predictions of secondary structure : C, H and E)

	# agrv 0: query.py
	# agrv 1: mfasta file (input)
	# agrv 2: ss2 file (input)
	# agrv 3: pssm_ss2 file (output)

# Things have to be understood and performed
  1. can't run psipred using a given mfasta file (makemat step), I just used input fasta to perform psipred_single. However, it's impossible to use this for template.
  2. query_ss2 doesn't take into account the neither seq weight nor pseudo-count

