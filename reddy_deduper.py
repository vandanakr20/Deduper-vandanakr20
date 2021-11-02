#!/usr/bin/env python

import argparse
import gzip
import sys
import re

# use argparse to get user input of the file names, if the data is paired, and the umi file
def get_args():
    parse = argparse.ArgumentParser(description="A program to remove the PCR duplicates from a chromosome sorted SAM file")
    parse.add_argument("-f", "--file", help="absolute file path of sorted SAM file", required=True, type=str)
    parse.add_argument("-p", "--paired", help="designate if file is paired end reads, pass -p if paired", required=False, action='store_true')
    parse.add_argument("-u", "--umi", help="absolute file path of list of UMI file", required=True, type=str)
    return parse.parse_args()

args = get_args()

#quit the program if the paired option is set
if args.paired:
    print('This program does not yet have the ability to remove PCR duplicates from a file of paired end reads. Sorry :(')
    sys.exit(1)

def fix_start(cigar: str, start: str, strand: str, length: str) -> int:
        
    '''
    This function will take a string of the cigar string, a string of the start position, the strand and the length of the read and return the updated 5' start position as an int after accounting for soft clipping on both strands. 

    Input: 3S5M, 13, F, 100
    Expected output: 10

    parameter: string cigar, string start, string strand, string length
    return: int corr_start 

    return int
    '''
    corr_start = start
    clip_amt = 0
    insert_amt = 0
    splice_amt = 0
    delete_amt = 0
    s_match = 0
    split_match = {}

    #break cigar string into list of tuples (#, letter)
    pattern = "([0-9]+)([ISNDMHP])" #do i really need the MHP?
    list_match = re.findall(pattern, cigar[:-1])

    #turn it into a dict, with a letters as keys and values as a list of digits
    for item in list_match:
        if item[1] in split_match.keys():
            split_match[item[1]].append(int(item[0]))
        else:
            split_match[item[1]] = [int(item[0])]
    
    #find all the true values 
    for value in split_match:
        if value == 'S':
            clip_amt = split_match['S'][0]
        if value == 'I':
            insert_amt = sum(split_match['I'])
        if value == 'N':
            splice_amt = sum(split_match['N'])
        if value == 'D':
            delete_amt = sum(split_match['D'])

    #correct the strand based on direction
    if strand == 'F':
            #strand is forward
            corr_start = int(start) - clip_amt
    elif strand == 'R':
            #strand is reverse and need the 5' end 
            correct = int(length) - clip_amt - insert_amt 
            corr_start = int(start) + correct + delete_amt + splice_amt
    
        
    return corr_start

    
def extract(line: str) -> list:
    
    '''
    This function will take a string and return of list of 5 items from the line, UMI, Start Position, Chromosome Number, Cigar String, and Strand (bitwise flag)

    Input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	5M	*	0	0	TCCAG	6AEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
    Expected output: ['CTGTTCAC', '2', '5M', 'R', '76814284']

    parameter: string line
    return: list of strings ret_list

    return list
    '''

    ret_list = []

    # extract UMI
    qname = line.split()[0]
    umi = qname[-8:]
    ret_list.append(umi)
    
    #extract chromosome 
    chrom = line.split()[2]
    ret_list.append(chrom)

    #extract cigar 
    cigar = line.split()[5]
    ret_list.append(cigar)

    #extract strand and set string to match 
    #forward is F and reverse is R
    strand = line.split()[1]
    if ((int(strand) & 16) == 16):
        strand = 'R'
    else:
        strand = 'F'
    ret_list.append(strand)
   
    #extract start position and pass it to fix start to correct start position 
    start = line.split()[3]
    fixed = fix_start(cigar, start, strand, len(line.split()[9]))
    ret_list.append(str(fixed))

    #return list [UMI, chrom, cigar, strand, position]
    return ret_list

    
def check(rec: list, all_recs: dict) -> bool:
    
    '''
    This function will take a list of strings of the sUMI, chrom, cigar, strand, and corrected position and a dictionary, and return a bool (True/False) if record is in dict

    Input: ['CTGTTCAC', '76814284', '2', '5M', '36'], {7681: 2CTGTTCACR}
    Expected output: False 

    parameter: list of strings, dictionary 
    return: bool 

    return bool 
    '''
    #create string for value
    value = rec[1]+rec[0]+rec[3]
 
    #check if both key and value match a record in the dict
    if rec[4] in all_recs.keys():
        for sting in all_recs[rec[4]]:
            if sting == value:
                return True
        return False
    else:
        return False


#create the set of Umis
umi_set = set()
with open(args.umi, 'r') as fh:
    for line in fh:
        umi_set.add(line.strip())

#open reading and writing file
fhr = open(args.file, "r")
fhw = open('C1_SE_uniqAlign_deduped.sam', 'w')


prev_chrom = '0'
all_recs = {}
in_dict = True
value = ''
bad_umi = 0

for line in fhr:
    #write out header lines to output file
    if line.startswith('@'):
        fhw.write(line)
    else:
        #extract useful info from line
        rec = extract(line)

        #check if umi is in the set of known umis
        if rec[0] in umi_set:

            #check if record is in dict
            in_dict = check(rec, all_recs)
            if in_dict == False:
                # if record is not in dict
                fhw.write(line) #write out line
                value = rec[1]+rec[0]+rec[3]

                #if new chromosome, clear dict and start over
                if prev_chrom != rec[1]:
                    prev_chrom = rec[1]
                    all_recs = {}
                    all_recs[rec[4]] = [value]
                
                #else, add to existing dict
                else:
                    if rec[4] in all_recs.keys():
                        all_recs[rec[4]].append(value)
                    else:
                        all_recs[rec[4]] = [value]
            else:
                #record is in dict --> ignore 
                continue
        else:
            bad_umi += 1

#print(bad_umi)
#close the files
fhr.close()
fhw.close()

