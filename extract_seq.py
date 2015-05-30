#coding:utf-8
__author__ = 'shawn_au'
"""This script extracts the query miRNA from any database to extract_output.txt
*Input:extract_input.txt(windows)
*Output:extract_output.txt(linux)

Beta:2015-4-22
Release:2015-5-15
This line was left blank."""
import re


def extraction():
    print('Start searching...')
    with open("input/query.txt", "r+") as query:
        query_list = query.readlines()
    with open("database/mature.txt", "r+") as seq:
        seq_list = seq.readlines()
    with open("output/extract_output.txt", "w+") as tar:
        match = 0
        failed = 0
        for i in range(len(query_list)):  # iterate from query list
            for j in range(0, len(seq_list), 2):  # iterate from database
                processed = (query_list[i])[:-2]  # -2 is used in txt under windows, -1 used under linux
                reg = re.compile('>'+processed+'\D'+'.*', re.I)  # use '\D' to prevent 123 matches 12345
                searchobj = reg.search(seq_list[j])
                if searchobj:
                    full_name = re.search(processed+'\S*', searchobj.group(), re.I)  # cut the name
                    tar.write(full_name.group()+'\n')   # output the query in *(not)* FASTA format
                    tar.write(seq_list[j+1])  # output the sequence
                    match += 1
                else:
                    failed += 1  # count failed matches
    print('Done. Find '+str( match)+' matches.')
    print('Results saved into /output/extract_output.txt')
    return 0
