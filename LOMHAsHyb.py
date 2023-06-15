# Loci With Multiple Haplotypes (LOMHAs) Detection From A Population of A Hybrid Species

# Get libraries
import pandas as pd
import glob

# import the FASTA files of the hybrid indvs. with .fa extension
h0_list = glob.glob("*.fa")
h0 = pd.DataFrame()
for h_file in h0_list:
    h0_pr = pd.read_csv(h_file, sep="\t", header=None)
    h0 = h0.append(h0_pr)


    # Make the data frame of fasta a two-column data frame (i.e. Col.1:The fasta tag, Col.2:The sequence)
    rowIndex = h0.index[:]
    idx = rowIndex.values.tolist()

    even_idx = [x for x in idx if x%2 == 0]
    odd_idx = [x for x in idx if x%2 != 0]

    h0DF = pd.DataFrame(columns=['SeqNames', 'Seq'])

    h0DF['SeqNames'] = list(h0.loc[even_idx, 0])
    h0DF['Seq'] = list(h0.loc[odd_idx, 0])

# Change the SeqNames in the 1st column as a locus number in the hybrid data frame
for seqName in range(0, len(h0DF)):
  h0DF.at[seqName, 'SeqNames'] = h0DF.iloc[seqName, 0].split('_')[5]


# Filter the rows including 'N' in Seq column
h0DFNo_N = h0DF[~h0DF.Seq.str.contains('N')]

# Filter rows contain duplicate seqs in the Seq column
h0DFNo_Dup = h0DFNo_N.drop_duplicates(subset=['Seq'])

# Sort the h0DFNo_Dup file as to locus numbers
h0DFNo_DupSort = h0DFNo_Dup.sort_values(by=['SeqNames'], ascending=False)


# Detect tri or more haplotypic states through the hybrid indvs.
# and save them in the "h0DF*" output object
for line in h0DFNo_DupSort['SeqNames'].unique():
    lidx = h0DFNo_DupSort.index[h0DFNo_DupSort['SeqNames'] == str(line)].tolist()
    if len(lidx) < 3 :
        h0DFNo_DupSort.drop(lidx, inplace=True)

# Export the result file
h0DFNo_DupSort.to_csv('h0DFMultiHaplotypes', sep = '\t', index=False, header=None)
