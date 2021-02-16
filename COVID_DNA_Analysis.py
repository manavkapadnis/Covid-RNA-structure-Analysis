##########################################################################################################################
#This python script does the DNA analysis of COVID virus. It also does the phylogenetic analysis of the COVID-19 DNA #####
#It uses genebank to read COVID DNA files and does the analysis								 #
#Written by Jitendra Moon												 #
##########################################################################################################################

from __future__ import division
#Import Necessary Libraries linear algebra
import numpy as num_lib 
#Import Necessary Libraries linear algebra
#data processing, CSV file I/O (e.g. pd.read_csv)
import pandas as data_lib 
import matplotlib.pyplot as plot_lib
import seaborn as sns
import os
from Bio import SeqIO,Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.SeqUtils import ProtParam
from Bio import pairwise2
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
data_lib.plotting.register_matplotlib_converters()

def print_dna_sequence():
	for dna_sequence in SeqIO.parse('MN908947.fna', "fasta"):
		print(dna_sequence.id)
		print(dna_sequence.seq)
		print(len(dna_sequence),'nucliotides')
		
def convert_dna_rna(COVID_DNA):
	COVID_mRNA=COVID_DNA.transcribe()
	print("\n COVID RNA is :\n",COVID_mRNA)
	print("\n Length of COVID RNA is :\n",len(COVID_mRNA))
	return COVID_mRNA
	
def get_animo_acid_sequence(COVID_mRNA):
	COVID_amino_Acid = COVID_mRNA.translate(table=1, cds=False)
	print('COVID Amino Acid : ', COVID_amino_Acid)
	print("Length of COVID Protein : ",len(COVID_amino_Acid))
	print("Length of COVID Original mRNA : ",len(COVID_amino_Acid))
	return COVID_amino_Acid
	
def protein_analysis(COVID_proteins):
	protein_of_interest_list = []
	Molecular_Weight_list = []
	for protein_record in COVID_proteins[:]:
		print("\n")
		X = ProtParam.ProteinAnalysis(str(protein_record))
		POI = X.count_amino_acids()
		protein_of_interest_list.append(POI)
		MW = X.molecular_weight()
		Molecular_Weight_list.append(MW)
		print("Protein of Interest = ", POI) 
		print("Amino acids percent = ", str(X.get_amino_acids_percent())) 
		print("Molecular weight = ", MW)
		print("Aromaticity = ", X.aromaticity()) 
		print("Flexibility = ", X.flexibility()) 
		print("Isoelectric point = ", X.isoelectric_point()) 
		print("Secondary structure fraction = ", X.secondary_structure_fraction())
	MoW = data_lib.DataFrame(data = Molecular_Weight_list,columns = ["Molecular Weights"] )
	print("\n Molecular weight DF is :\n", MoW.head())
	return protein_of_interest_list
	
def plot_poi_list(protein_of_interest_list):
	poi_list = protein_of_interest_list[48]
	plot_lib.figure(figsize=(10,6));
	plot_lib.bar(poi_list.keys(), list(poi_list.values()), align='center')
	plot_lib.show()
	
	# Plot lengths
	#plot_lib.figure(figsize=(20,5))
	#plot_lib.subplot(111)
	#plot_lib.hist(functional_proteins['length'])
	#plot_lib.title('Length of proteins -- histogram')
	# Remove the extremes
	#plot_lib.figure(figsize=(20,5))
	#wo = functional_proteins.loc[functional_proteins['length'] < 60]
	#plot_lib.subplot(121)
	#plot_lib.hist(wo['length'])
	#plot_lib.title('Lenght of proteins (where < 60)')
	
	#wo = functional_proteins.loc[functional_proteins['length'] > 1000]
	#plot_lib.subplot(122)
	#plot_lib.hist(wo['length'])
	#plot_lib.title('Length of proteins (where > 1000)')
	
	# See what's about that huge protein
	#large_prot = functional_proteins.loc[functional_proteins['length'] > 2700]
	#l = large_prot['sequence'].tolist()[0]
	#print('Sequence sample:', '...',l[1000:1150],'...')
	
def draw_COVID_genome_diagram():
	Entrez.email="test@test.com"
	handle=Entrez.efetch(db='nuccore',id='NC_045512.2',rettype='gb',retmode='text')
	COVID_DNA_sequence=SeqIO.read(handle,"genbank")
	handle.close()
	print("\n DNA Sequence is :\n",COVID_DNA_sequence)
	gd_covid_diagram=GenomeDiagram.Diagram("Severe Acute resoiratory Syndrome coronavirus 2 isolate Wuhan-Hu-1")
	gd_track_covid_features=gd_covid_diagram.new_track(1,name="Annoted Features")
	gd_covid_feature_set=gd_track_covid_features.new_set()
	for covid_dna_sequence in COVID_DNA_sequence.features:
		if covid_dna_sequence.type != "gene":
			continue
		if len(gd_covid_feature_set) % 2 == 0:
			color=colors.blue
		else:
			color=colors.green
		gd_covid_feature_set.add_feature(covid_dna_sequence,sigil="ARROW",color=color,label=True,label_size=14,label_angle=0)
	gd_covid_diagram.draw(		
	format="linear",
	orientation="landscape",
	pagesize=(30*cm,30*cm),
	fragments=4,
	start=0,
	end=len(COVID_DNA_sequence),
	)
	gd_covid_diagram.write("COVID_linear.pdf","PDF")
	gd_covid_diagram.write("COVID_linear.png","PNG")
	
	gd_covid_diagram.draw(
	format="circular",
	circular=True,
	pagesize=(30*cm,30*cm),
	start=0,
	end=len(COVID_DNA_sequence),
	circle_core=0.7,
	)
	
	gd_covid_diagram.write("COVID_circular.pdf","PDF")
	
print_dna_sequence()	

COVID_DNA_sequence=SeqIO.read('MN908947.fna', "fasta")
COVID_DNA = COVID_DNA_sequence.seq
print("\n COVID DNA Seqneuce is:\n", COVID_DNA)



#The difference between the complementary DNA and the mRNA is just that the 
#bases T (for Thymine) is replaced with U (for Uracil).
COVID_RNA=convert_dna_rna(COVID_DNA)
COVID_amino_acid=get_animo_acid_sequence(COVID_RNA)

#Codons Cells decode mRNAs by reading their nucleotides in groups of three, called codons. 
#Here are some features of codons:
#Most codons specify an amino acid.Three "stop" codons mark the end of a protein
#One "start" codon, AUG, marks the beginning of a protein and also encodes the amino acid methionine
print(CodonTable.unambiguous_rna_by_name['Standard'])

#Let's now identify all the polypeptides so basically separating at the stop codon, marked by * . 
#Then let's remove any sequence less than 20 amino acids long, as this is the smallest known functional protein (if curious). 
#Note: In humans the smallest known functional protien is 44 amino acids long.

#Identify all the Proteins (chains of amino acids)
COVID_proteins = COVID_amino_acid.split('*') # * is translated stop codon
for index in COVID_proteins[:]:
	if len(index) < 20:	
		COVID_proteins.remove(index)
		
protein_of_interest_list=protein_analysis(COVID_proteins)
plot_poi_list(protein_of_interest_list)

SARS = SeqIO.read("sars.fasta", "fasta")
MERS = SeqIO.read("mers.fasta", "fasta")
COV2 = SeqIO.read("cov2.fasta", "fasta")

print('Sequence Lengths:')
print('SARS:', len(SARS.seq))
print('COV2:', len(COV2.seq))
print('MERS:', len(MERS.seq))

#Visualize DNA sequence using Squiggle 
os.system("Squiggle cov2.fasta sars.fasta mers.fasta --method=gates --separate")

SARS_COV = pairwise2.align.globalxx(SARS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('SARS/COV Similarity (%):', SARS_COV / len(SARS.seq) * 100)
MERS_COV = pairwise2.align.globalxx(MERS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('MERS/COV Similarity (%):', MERS_COV / len(MERS.seq) * 100)
MERS_SARS = pairwise2.align.globalxx(MERS.seq, SARS.seq, one_alignment_only=True, score_only=True)
print('MERS/SARS Similarity (%):', MERS_SARS / len(SARS.seq) * 100)

# Plot the data
X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
Y = [SARS_COV/ len(SARS.seq) * 100, MERS_COV/ len(MERS.seq)*100, MERS_SARS/len(SARS.seq)*100]
plot_lib.title('Sequence identity (%)')
plot_lib.bar(X,Y)
plot_lib.show()

draw_COVID_genome_diagram()