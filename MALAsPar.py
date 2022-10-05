# MULTI ALLELLIC LOCUS ALLELES (MALAs) DETECTION IN THE PARENTAL POPULATIONS OF A HYBRID SPECIES

# Get the libraries
import pandas as pd
import sys

# Take the maternal and paternal fasta files including the all data of all individuals
with open(sys.argv[1], "r") as p1:
    with open(sys.argv[2], "r") as m1:

        # Load the data
        # Create Pandas object from input fasta files
        p0 = pd.read_csv(p1, sep="\t", header=None)
        m0 = pd.read_csv(m1, sep="\t", header=None)

        # Convert the single-columns fasta format to two-column data frame
        SeqNamesP = []
        SeqP = []

        for line in range(0, len(p0), 2):
            SeqNamesP.append(p0.iloc[line,0])
            SeqP.append(p0.iloc[line+1,0])

        p0DF = pd.DataFrame(zip(SeqNamesP, SeqP), columns=['SeqNames', 'Seq'])

        SeqNamesM = []
        SeqM = []

        for line in range(0, len(m0), 2):
            SeqNamesM.append(m0.iloc[line,0])
            SeqM.append(m0.iloc[line+1,0])

        m0DF = pd.DataFrame(zip(SeqNamesM, SeqM), columns=['SeqNames', 'Seq'])

        # Change the seq names as Locus Numbers in 1st column of
        # the data frame of parentals
        for seqName in range(0, len(m0DF)):
            m0DF.at[seqName, 'SeqNames'] = m0DF.iloc[seqName, 0].split('_')[5]

        for seqName in range(0, len(p0DF)):
            p0DF.at[seqName, 'SeqNames'] = p0DF.iloc[seqName, 0].split('_')[5]

        # Export the maternal and paternal files with two-column format
        m0DF.to_csv('m0DF_twoCols', sep = '\t', index=False, header=None)
        p0DF.to_csv('p0DF_twoCols', sep = '\t', index=False, header=None)



## Hybrid code part
# Get the hybrid MALAs file from the cureent working directory
import glob
for h_file in glob.glob("h0DFMultiAlleles"):
    h0DF = pd.read_csv(h_file, sep="\t", header=None)
    h0DF.rename(columns={0: 'SeqNames', 1:'Seq'}, inplace=True)


# Detect if each multi-allele of each MALAs locus in the hybrids is present in parental populations
# Create output files PatDF for MALAs in paternal population and MatDF for MALAs in maternal population 
# and Mat_PatDF for MALAs in both maternal and paternal populations
PatDF = pd.DataFrame(columns=['SeqNames', 'Seq'])
MatDF = pd.DataFrame(columns=['SeqNames', 'Seq'])
Mat_PatDF = pd.DataFrame(columns=['SeqNames', 'Seq'])


for sq in range(0, len(h0DF)):
    # Inspect if the hybrid allele is in paternal side 
    inspectorPh = pd.DataFrame([h0DF.loc[sq,'Seq']]).isin(list(p0DF.loc[:, 'Seq'])).any().any()
    # Inspect if the hybrid allle is in maternal side
    inspectorMh = pd.DataFrame([h0DF.loc[sq,'Seq']]).isin(list(m0DF.loc[:, 'Seq'])).any().any()

    # Save the MALAs as maternal and paternal (if any) in the output files
    if inspectorPh==True:
        PatDF = PatDF.append(pd.DataFrame(h0DF.loc[sq,:]).T, ignore_index=True)

    if inspectorMh==True:
        MatDF = MatDF.append(pd.DataFrame(h0DF.loc[sq,:]).T, ignore_index=True)

    if inspectorMh==True and inspectorPh==True:
        Mat_PatDF = Mat_PatDF.append(pd.DataFrame(h0DF.loc[sq,:]).T, ignore_index=True)

# Export the paternal and maternal MALAs files of the hybrid species
PatDF.to_csv("./P_h0DFMultiAlleles", sep = '\t', index=False, header=None)
MatDF.to_csv("./M_h0DFMultiAlleles", sep = '\t', index=False, header=None)
Mat_PatDF.to_csv("./MP_h0DFMultiAlleles", sep = '\t', index=False, header=None)
