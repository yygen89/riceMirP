The training and the testing datasets are included in this folder.

training set:

1. mirPlantPre19.txt: All non-redundant plant pre-miRNAs published in miRbase 19 were downloaded from miRBase 19.0

2. mirPlantMat19.txt: All non-redundant plant miRNAs published in miRbase 19 were downloaded from miRBase 19.0

3. mirPlantPre19_single.txt: All non-redundant plant miRNAs and their pre-miRNAs with stem-loop structure in the mirPlantPre19.txt. Each of miRNA is at the front of its pre-miRNA.

testing set:

1. mirPlantPre20.txt: All non-redundant plant pre-miRNAs published in miRbase 20 were downloaded from miRBase 20.0

2. mirPlantMat20.txt: All non-redundant plant miRNAs published in miRbase 20 were downloaded from miRBase 20.0

3. mirPlantPre20_single.txt: All non-redundant plant miRNAs and their pre-miRNAs with stem-loop structure in the mirPlantPre20.txt. Each of miRNA is at the front of its pre-miRNA.

negative set:

1. negData.txt: Pseudo pre-miRNAs we selected from cDNAs. It is found that most of known plants pre-miRNAs are 120nt in length. Thus, a sliding window of width ranging from 60nt to 150nt randomly was used to scan the CDSs to produce sequence segments. The sequence segments should be folded into single stem-loop structures. Further, they should satisfy five criteria on the number of base pairs in hairpins, %G+C, MFEI, complementary base pairing on mature miRNAs and the stability of the precursor related to MFE rate. The criteria are determined by observing real intercepted plant pre-miRNAs. The criteria for selecting the pseudo miRNAs in length are: minimum of 19 base pairings in hairpin structure, %G+C>0.242 and <0.825, MFEI >0.522 and <1.39, no me miRNAs.