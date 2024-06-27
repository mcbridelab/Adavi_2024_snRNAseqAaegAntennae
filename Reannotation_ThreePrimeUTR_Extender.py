#!/usr/bin/env python 
#usage: python gffparse_vitor.py your.gff bedtools_intersect.txt
#example: python ./nhr_sandbox/gffparse_vitor.py mismatchCorrected_mergedManualVectorbaseAnnotation.gff 250
# python ./nhr_sandbox/gffparse_vitor.py mismatchCorrected_mergedManualVectorbaseAnnotation.gff 250 &
# python ./nhr_sandbox/gffparse_vitor.py mismatchCorrected_mergedManualVectorbaseAnnotation.gff 500 &
# python ./nhr_sandbox/gffparse_vitor.py mismatchCorrected_mergedManualVectorbaseAnnotation.gff 750


# Read in libraries
import numpy
import gffutils
import sys
import sqlite3
import subprocess
from collections import OrderedDict


# Read in terminal arguments into local variables
currgff=sys.argv[1] # 'mismatchCorrected_mergedManualVectorbaseAnnotation.gff' #sys.argv[1]
extension=sys.argv[2] # 250 #sys.argv[2]
# intersect=sys.argv[3] # 'utrSubsetUtrs250.gff' #sys.argv[3]

#### Make or read in database's ####
# Write database if first time
annotation=gffutils.create_db(currgff,dbfn='annotation'+str(extension)+'.db',force=True,keep_order=True, sort_attribute_values=True)
# stringtie_trans=gffutils.create_db('utr250_transcripts.gff',dbfn='stringtie_trans'+str(extension)+'.db',force=True,keep_order=False, sort_attribute_values=True)

# # Read database if not first time
# annotation=gffutils.FeatureDB('annotation.db')
# stringtie_trans=gffutils.FeatureDB('stringtie_trans.db')
####.####


# #### NOTE: use subprocess.run to run a bash command to assemble stringtie models ####
# # Assemble the UTR with stringtie
# bashCommand = "stringtie -p 40 -m 30 -j 1000000000 --fr ./utrSubsetReads" + str(extension) + ".bam > utrSubsetUtrs" + str(extension) + ".gff"
# # Run command with subprocess.run
# subprocess.run(bashCommand, shell=True, check=True)
# ####.####




#### Wait for stringtie model to finish before creating database####
# stringtie=gffutils.create_db("utrSubsetUtrs"+str(extension)+".gff", dbfn='stringtie'+str(extension)+'.db',\
#     force=True, keep_order=True, sort_attribute_values=True, disable_infer_transcripts=True)
# # Read database if not first time
stringtie=gffutils.FeatureDB('stringtie'+str(extension)+'.db')
####.####




#### Defines function to test for strandedness mismatch between a feature and its children ####
def countMismatch(db):
    count=0
    mismatched_features=[]
    # Loop through all features
    for feature in db.all_features():
        # Get associated children for all features
        children=list(db.children(id=feature.id))
        # Append strand information into temp variable of which to check via for loop for consistency across relevent features
        strand_temp=list(feature.strand)
        for strand_child in children:
            strand_temp.append(strand_child.strand)
        # Check for strand mismatch then print ID of parent feature 
        #       and append mismatched features to a list
        if len(numpy.unique(strand_temp)) != 1:
            print(feature.id)
            mismatched_features.append(feature)
    return(mismatched_features)


##### Module to change the gene strand if it mismatches both mRNA and exon strands####
# Set a boolean flag to keep track of whether correction happened as a result of mismatch
corrected=False
# Use user defined function to check if there are any strand mismatches between genes and their children
miss=countMismatch(annotation)
count=0
for m in miss:
    count+=1
    m.id

print("Number of mismatched strand genes before correction="+str(count))

if (count > 0):
    corrected=True
    temp="mismatchCorrected_"+currgff
    f = open(temp, "w")
    # Module to change the gene strand if it mismatches both mRNA and exon strands
    for feature in annotation.all_features():
        # Check if feature is of type gene
        if feature.featuretype in ['gene','protein_coding_Gene_ID']:
            # Get associated children of parent
            mRNA=list(annotation.children(id=feature, featuretype='mRNA'))
            rest=list(annotation.children(id=feature))
            # Append strand information into temp variable of which to check via for loop for consistency across relevent features
            strand_temp=list(feature.strand)
            for strand_mRNA in mRNA:
                strand_temp.append(strand_mRNA.strand)
            for strand_rest in rest:
                strand_temp.append(strand_rest.strand)
            if len(numpy.unique(strand_temp)) <= 2:
                feature.strand=strand_mRNA.strand
            f.write(str(feature)+'\n')
        else:
            f.write(str(feature)+'\n')
    f.close()

if (count > 0):
    # Use user defined function to check if there are any strand mismatches between genes and their children
    ann=gffutils.create_db(temp,dbfn='correctedAnnotation.db',force=True,keep_order=True, sort_attribute_values=True)
    miss=countMismatch(ann)
    count=0
    for m in miss:
        count+=1
        m.id
    print("Number of mismatched strand genes after correction="+str(count))
##### Module to change the gene strand if it mismatches both mRNA and exon strands####









#### MERGING LOGIC ####
# Defines the maximum extension allowed
if (corrected == True):
    temp="mismatchCorrected_"+currgff
    ann=gffutils.create_db(temp,dbfn='correctedAnnotation.db',force=True,keep_order=True, sort_attribute_values=True)
else:
    temp=currgff
    ann=annotation

## Sub module to get a lis of the features needed to handel every case correctly by moving from gene 
#               level to transcript level to features that need to be updated if there is an extension
# Get a list of all of the features
allFeature=[]
for i in annotation.featuretypes():
    print(i)
    allFeature.append(i)

# Gene level features
gene_level=['gene','pseudogene','protein_coding_Gene_ID','ncRNA_Gene_ID','pseudoGene_ID']
print(gene_level)
# Transcript level features
transcripts_level=['transcript','pseudogenic_transcript','rRNA','snRNA','snoRNA','tRNA','mRNA','RNase_MRP_RNA','RNase_P_RNA','SRP_RNA','lnc_RNA','ncRNA','pre_miRNA']
print(transcripts_level)
# Features that should be updated if gene model changes
to_update=['exon','three_prime_UTR']
print(to_update)
# 'Missing' features
missing=[]

# Prints out a list of features not being changed
al=[]
for i in gene_level:
    # print(i)
    al.append(i)

for i in transcripts_level:
    # print(i)
    al.append(i)

for i in to_update:
    # print(i)
    al.append(i)

for i in allFeature:
    if (i not in al):
        # print(i)
        missing.append(i)

print(missing)
##.##


# #### UNIT TEST LOADING IN ####
# temp='unit_test_merge.gff'
# ann=gffutils.create_db('unit_tests_gene.gff', dbfn='unit_annotation.db', force=True, keep_order=True, sort_attribute_values=True)
# stringtie=gffutils.create_db('unit_tests_stringtie.gff', dbfn='unit_stringtie.db', force=True, keep_order=True, sort_attribute_values=True)
# #### UNIT TEST LOADING IN ####


# #### LINE FOR TESTING A FEW GENE MODELS####
# ann=gffutils.create_db('oneGeneTest.gff',dbfn='test_gene.db',force=True,keep_order=True, sort_attribute_values=True)
# #### LINE FOR TESTING A FEW GENE MODELS####


#### RECORD OF NOTES I MADE DURING CODE GENERATION
# NOTE: figure out what to do about the double exon issue
# THE PROBLEM: in the original file there is fewer exons than CDS features because some exons have two transcript parent IDs associated with it
# For instance:
# # AaegL5_1	VEuPathDB	exon	291228668	291228799	.	+	.	ID=exon_AAEL003267_CG16791-RB-E12;Parent=AAEL003267_CG16791-RB,AAEL003267_CG16791-RC;Gene_ID_id=AAEL003267_CG16791
 
# NOTE: The extra weird models is actually only because of a gff vs gtf formatting issue. In the gtf version of the merge file those 'extra' models go away
#### RECORD OF NOTES I MADE DURING CODE GENERATION




#### Read in all features into an ordereddict python object which I will search by key to change start and stop columns####
gffOrderedDict=OrderedDict()
for fe in ann.all_features():
    gffOrderedDict[fe.id]=fe
####.####




## Start of merging logic
exception=False
for gene in ann.all_features():
    # TODO: REMOVE THIS CHECK MAYBE?
    if (exception == True):
        break
    if gene.featuretype in gene_level:
        # Make empty list that everything trancript and below will be appended to and printed after transcript for loop
        tran_print=[]
        # Flag to check if gene changes
        gene_changed=False 
        # Define an impossible location as temp variable
        new_length=-1 
        # Loop through the associated transcripts from gene level
        for trans in ann.children(id=gene.id, featuretype=transcripts_level):
            # Flag to check if transcript changes
            transcript_changed=False
            # Flag to check if current transcript has a 3'UTR or not
            hasUTR=False
            # Check if a stringtie model exists at 3' end of gene + N bases
            match_number=0
            if (trans.strand == '+'):
                # Save original stopping point to recover later if exon overlaps other gene's exons 
                original_stop=trans.stop
                for i in stringtie.region(seqid=trans.seqid, start=(trans.stop-1), end=(trans.stop+int(extension)), strand=trans.strand, completely_within=False, featuretype='transcript'):
                    match_number+=1
                if (match_number > 0):
                    matches=stringtie.region(seqid=trans.seqid, start=(trans.stop-1),\
                        end=(trans.stop+int(extension)), strand=trans.strand, completely_within=False, featuretype='transcript')
            elif (trans.strand == '-'):
                # Save original starting point to recover later if exon overlaps other gene's exons 
                original_start=trans.start 
                for i in stringtie.region(seqid=trans.seqid, start=(trans.start-int(extension)), end=(trans.start+1), strand=trans.strand, completely_within=False, featuretype='transcript'):
                    match_number+=1
                if (match_number > 0):
                    matches=stringtie.region(seqid=trans.seqid, start=(trans.start-int(extension)),\
                        end=(trans.start+1), strand=trans.strand, completely_within=False, featuretype='transcript')
            # Logic to get the new (longest) possible stopping point
            #### Positive strand code
            if ((trans.strand == '+') and (match_number > 0)):
                # Save original stopping point to recover later if exon overlaps other gene's exons 
                original_stop=trans.stop
                for match in matches:
                    if (new_length < match.stop or new_length == -1):
                        new_length=match.stop
                # Check if new stop is greater than current and save if true
                if (trans.stop < new_length and new_length != -1):
                    transcript_changed=True
                    # Checks if the proposed extension is greater than the maximum allowed
                    if ((new_length - trans.stop) > int(extension)):
                        new_length=(trans.stop+int(extension))
                    # Find the correct exon to extend 
                    toUpdate=ann.children(id=trans.id, featuretype=to_update,order_by='end')
                    for fe in toUpdate:
                        if(fe.featuretype == 'exon'):
                            last_exon=fe
                        if(fe.featuretype == 'three_prime_UTR'):
                            hasUTR=True
                            last_3utr=fe
                    # Check if the exon chosen is correct
                    if (trans.stop != last_exon.stop):
                        print("last_exon DOES NOT HAVE CURRENT TRANSCRIPTS FINAL BP")
                        print(trans)
                        print(last_exon)
                        print(trans.stop, last_exon.stop)
                        print(trans.stop == last_exon.stop)
                        exception=True
                    # Check if the three_prime_UTR chosen is correct
                    if (hasUTR == True):
                        if (trans.stop != last_3utr.stop):
                            print("three_prime_UTR DOES NOT HAVE CURRENT TRANSCRIPTS FINAL BP")
                            print(trans)
                            print(last_3utr)
                            print(trans.stop, last_3utr.stop)
                            print(trans.stop == last_3utr.stop)
                            exception=True
                    ## Implemented logic for checking if another gene's exon overlaps new boundary
                    # Make empty list of overlapping exons that do not match gene.id
                    transcript_conflict=[]
                    exon_conflict=[]
                    # Reset match number and check that there is any relevent matchs
                    match_number=0
                    # Loop through transcripts overlapping last_exon
                    for t in ann.region(seqid=last_exon.seqid, start=trans.stop, end=new_length, strand=last_exon.strand, featuretype=transcripts_level, completely_within=False):
                        # Check if the parent of the transcripts do not match the current gene's transcripts
                        parent=ann.parents(id=t.id)
                        for p in parent:
                            if (p.id != gene.id):
                                match_number+=1
                                transcript_conflict.append(t.id)
                    # Check if overlapping exons are from transcripts not associated with the given gene
                    if (match_number > 0):
                        for ex in ann.region(seqid=last_exon.seqid, start=trans.stop, end=new_length, strand=last_exon.strand, featuretype='exon', completely_within=False):
                            parent=ann.parents(id=ex.id)
                            for p in parent:
                                if (p.id in transcript_conflict):
                                    exon_conflict.append(ex)
                    # Take the overlapping exons from said transcripts
                    for ex in exon_conflict:
                        if (ex.start < new_length):
                            new_length=(ex.start-1)
                            if (new_length < original_stop):
                                new_length=original_stop
            #### Negative strand code
            if ((trans.strand == '-') and (match_number > 0)):
                # Save original stopping point to recover later if exon overlaps other gene's exons 
                original_start=trans.start
                for match in matches:
                    if (new_length > match.start or new_length == -1):
                        new_length=match.start
                # Check if new start is greater than current and save if true
                if (trans.start > new_length and new_length != -1):
                    transcript_changed=True
                    # Checks if the proposed extension is greater than the maximum allowed
                    if ((trans.start - new_length) > int(extension)):
                        new_length=(trans.start-int(extension))
                    # Find the correct exon to extend 
                    toUpdate=ann.children(id=trans.id, featuretype=to_update,order_by='start')
                    for fe in toUpdate:
                        if(fe.featuretype == 'exon'):
                            last_exon=fe
                            break
                        if(fe.featuretype == 'three_prime_UTR'):
                            hasUTR=True
                            last_3utr=fe
                            break
                    # Check if the exon chosen is correct
                    if (trans.start != last_exon.start):
                        print("last_exon DOES NOT HAVE CURRENT TRANSCRIPTS FINAL BP")
                        print(trans)
                        print(last_exon)
                        print(trans.start, last_exon.start)
                        print(trans.start == last_exon.start)
                        exception=True
                    # Check if the three_prime_UTR chosen is correct
                    if (hasUTR == True):
                        if (trans.start != last_3utr.start):
                            print("three_prime_UTR DOES NOT HAVE CURRENT TRANSCRIPTS FINAL BP")
                            print(trans)
                            print(last_3utr)
                            print(trans.start, last_3utr.start)
                            print(trans.start == last_3utr.start)
                            exception=True
                    ## Implemented logic for checking if another gene's exon overlaps new boundary
                    # Make empty list of overlapping exons that do not match gene.id
                    transcript_conflict=[]
                    exon_conflict=[]
                    # Reset match number and check that there is any relevent matchs
                    match_number=0
                    # Loop through transcripts overlapping last_exon
                    for t in ann.region(seqid=last_exon.seqid, end=trans.start, start=new_length, strand=last_exon.strand, featuretype=transcripts_level, completely_within=False):
                        # Check if the parent of the transcripts do not match the current gene's transcripts
                        parent=ann.parents(id=t.id)
                        for p in parent:
                            if (p.id != gene.id):
                                match_number+=1
                                transcript_conflict.append(t.id)
                    # Check if overlapping exons are from transcripts not associated with the given gene
                    if (match_number > 0):
                        for ex in ann.region(seqid=last_exon.seqid, end=trans.start, start=new_length, strand=last_exon.strand, featuretype='exon', completely_within=False):
                            parent=ann.parents(id=ex.id)
                            for p in parent:
                                if (p.id in transcript_conflict):
                                    exon_conflict.append(ex)
                    # Take the overlapping exons from said transcripts
                    for ex in exon_conflict:
                        if (ex.stop > new_length):
                            new_length=(ex.stop+1)
                            if (new_length > original_start):
                                new_length=original_start
            # TODO: IS THIS IF CHECK DOING ANYTHING? VARIABLE SET TO TRUE REGARDLESS IN BEGINNING OF THIS LOOP
            if (transcript_changed == True and trans.strand == '+'):
                # Replace transcript length with new longer length
                if (trans.stop < new_length):
                    gffOrderedDict[trans.id].stop=new_length##**
                # Check if new transcript length is greater than current gene model stop and update if true
                if (gene.stop < new_length):
                    gffOrderedDict[gene.id].stop=new_length##**
            if (transcript_changed == True and trans.strand == '-'):
                # Replace transcript length with new longer length
                if (trans.start > new_length):
                    gffOrderedDict[trans.id].start=new_length##**
                # Check if new transcript length is greater than current gene model stop and update if true
                if (gene.start > new_length):
                    gffOrderedDict[gene.id].start=new_length##**        
            # Extract all features from transcript to update and print them in for loop
            for update in ann.children(id=trans.id, featuretype=to_update, order_by='end'):
                # Check if transcript changed and therefore need to update associated features 
                if (transcript_changed == True):
                    # Update last exon
                    if (last_exon.id == update.id and last_exon.strand == '+'):
                        # Check if need to update all parents of last_exon as well:
                        for parent in ann.parents(id=last_exon.id, featuretype=transcripts_level):
                            if (parent.stop < new_length):
                                gffOrderedDict[parent.id].stop=new_length##**
                        gffOrderedDict[last_exon.id].stop=new_length##**
                        gffOrderedDict[last_exon.id].featuretype='3_'+last_exon.featuretype#'3_exon'##**
                    elif (last_exon.id == update.id and last_exon.strand == '-'):
                        # Check if need to update all parents of last_exon as well:
                        for parent in ann.parents(id=last_exon.id, featuretype=transcripts_level):
                            if (parent.start > new_length):
                                gffOrderedDict[parent.id].start=new_length##**
                        gffOrderedDict[last_exon.id].start=new_length##**
                        gffOrderedDict[last_exon.id].featuretype='3_'+last_exon.featuretype#'3_exon'##**
                    # Find UTR feature and update it with new stop point
                    elif (hasUTR == True):
                        if (last_3utr.id == update.id and last_3utr.strand == '+'):
                            gffOrderedDict[last_3utr.id].stop=new_length##**
                            gffOrderedDict[last_3utr.id].featuretype='3_'+last_3utr.featuretype##**
                        elif (last_3utr.id == update.id and last_3utr.strand == '-'):
                            gffOrderedDict[last_3utr.id].start=new_length##**
                            gffOrderedDict[last_3utr.id].featuretype='3_'+last_3utr.featuretype##**
                    else:
                        gffOrderedDict[update.id].featuretype='3_'+update.featuretype##**
                elif (transcript_changed == False):
                    gffOrderedDict[update.id].featuretype='3_'+update.featuretype##**
        # TODO: REMOVE THIS CHECK MAYBE?
        if (exception == True):
            break
        # Print first the genes then the transcript then the exons
        gffOrderedDict[gene.id].featuretype='1_'+gene.featuretype##**
    # Print the features not addressed in above code
    if gene.featuretype in missing:
        gffOrderedDict[gene.id].featuretype='4_'+gene.featuretype##**
    # TODO: REMOVE THIS CHECK MAYBE?
    if (exception == True):
        break




#### Print altered ordereddict into file that will be sorted and/or converterd ####
temp="mergeExtended"+str(extension)+"_"+currgff
f = open(temp, "w")
for key, value in gffOrderedDict.items():
    # print(str(value))
    f.write(str(value)+'\n')

f.close()

####.####




#### NOTE: use subprocess.run to sort gffutils output then feed into gffread to convort to gtf ####
# Assemble the UTR with stringtie
#** below
# bashCommand1 = "sort -k1,1V -k4,4n -k5,5rn -k3,3r mergeExtended"+str(extension)+"_"+currgff+" > sorted_mergeExtended"+str(extension)+"_"+currgff
# subprocess.run("sort mergeExtended"+str(extension)+"_"+currgff+" | uniq > small_sorted_mergeExtended"+str(extension)+"_"+currgff, shell=True, check=True)
# bashCommand1 = "sort -k1,1V -k5,5r -k4,4n -k3,3rn small_sorted_mergeExtended"+str(extension)+"_"+currgff+" | uniq > sorted_mergeExtended"+str(extension)+"_"+currgff
#** ^

# bashCommand1 = "sort -k1,1V -k5,5r -k4,4n -k3,3rn mergeExtended"+str(extension)+"_"+currgff+" > sorted_mergeExtended"+str(extension)+"_"+currgff
bashCommand1 = "gff3sort.pl mergeExtended"+str(extension)+"_"+currgff+" > sorted_mergeExtended"+str(extension)+"_"+currgff
bashCommand2 = "gffread -E sorted_mergeExtended"+str(extension)+"_"+currgff+" -T -o sorted_mergeExtended"+str(extension)+"_mergedManualVectorbaseAnnotation.gtf"
# Run command with subprocess.run
subprocess.run(bashCommand1, shell=True, check=True)
subprocess.run(bashCommand2, shell=True, check=True)
####.####




subprocess.run("gffread -E VectorBase-55_AaegyptiLVP_AGWG.gff -T -o VectorBase-55_AaegyptiLVP_AGWG.gtf", shell=True, check=True)
subprocess.run("gffread -E mismatchCorrected_mergedManualVectorbaseAnnotation.gff -T -o mismatchCorrected_mergedManualVectorbaseAnnotation.gtf", shell=True, check=True)
subprocess.run("gffread -E dmelJoinedvectorbaseAaeg.gff -T -o dmelJoinedvectorbaseAaeg.gtf", shell=True, check=True)


