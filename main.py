#-*- coding:utf-8 –*-
__author__ = 'chobits'
"""This script:
*Input format:dict.txt(linux)
*Output format:free_energy.csv (csv format)

*calculates the RNA cofold energy
*reorder results in energy order
*filtering results out of the miRNA in the same hairpin
*select the pair with the lower free energy in the same miRNA pairs
*plots it out in line form

beta:2015-4-22
release:2015-5-15
This line was left blank intentionally."""
import operator
import re
import RNA
import extract_seq
import cluster_distance
import plot_list
#  all of the functions accepts a list with the structure of [str, str, float]

extract_seq.extraction()

print('Calculating folding energy...')


def print_list_csv(lst, tar):
    for i in range(len(lst)):
        tar.write(str(lst[i][0])+','+str(lst[i][1])+','+str(lst[i][2])+'\n')


def screen(lst):  # screen out the pairs in the same hairpin
    screened_list_xy = []
    for i in range(len(lst)):
        tmp = re.compile('(?<=(hsa-\w\w\w-))\d*')
        fig1 = tmp.search(lst[i][0])  # to solve the bug: 17 will match 170:solved
        fig2 = tmp.search(lst[i][1])
        if not fig1.group() == fig2.group():
            screened_list_xy.append(lst[i])
    return screened_list_xy


def pick_lower(lst):  # pick the pair with the lower del G, can't use for with range() since del() exists
    m = 0
    while m < len(lst):
        n = m+1
        while n < len(lst):
            am = re.search('.*(?=-\dp$)', lst[m][0])
            bm = re.search('.*(?=-\dp$)', lst[m][1])
            an = re.search('.*(?=-\dp$)', lst[n][0])
            bn = re.search('.*(?=-\dp$)', lst[n][1])
            if (am and bm) and (an and bn):
                if (am.group() == bn.group()) and (bm.group() == an.group()):
                    if float(lst[m][2]) < float(lst[n][2]):  # attention: negative numbers, |a|>|b| means a<b
                        del lst[n]
                        n -= 1
                    else:
                        del lst[m]
                        m -= 1
                        break
            elif (am and an) and (bm == None) and (bn == None):
                if (am.group() == an.group()) and (lst[m][1] == lst[n][1]):
                    if float(lst[m][2]) < float(lst[n][2]):
                        del lst[n]
                        n -= 1
                    else:
                        del lst[m]
                        m -= 1
                        break
            elif (bm and bn) and (am == None) and (an == None):
                if (bm.group() == bn.group()) and (lst[m][0] == lst[n][0]):
                    if float(lst[m][2]) < float(lst[n][2]):
                        del lst[n]
                        n -= 1
                    else:
                        del lst[m]
                        m -= 1
                        break
            n += 1
        m += 1


def del_tail(lst):  # delete the -5p & -3p
    for i in range(len(lst)):
        for j in range(2):
            match = re.search('.*(?=-\dp$)', lst[i][j])
            if match:
                lst[i][j] = match.group()  # result_xy would be edited

with open("output/extract_output.txt", "r+") as database:
    db_list = database.readlines()

mirna = []  # 2D list store the mirna name & seq
result_xy = []  # 3D list store the mirna pairs with their del G


for i in range(0, len(db_list)-1, 2):  # notice the last line is blank
    mirna.append([db_list[i][:-1],  # -1 in linux, -2 in windows
                  db_list[i+1][:-1]
                  ])  # mirna:[name,seq]

for i in range(len(mirna)):
    for j in range(i, len(mirna)):
        result_xy.append([mirna[i][0],
                         mirna[j][0],
                         (RNA.cofold(mirna[i][1]+'&'+mirna[j][1]))[1]
                         ])

screened_result_xy = screen(result_xy)


pick_lower(screened_result_xy)
del_tail(screened_result_xy)  # result_xy would be edited
sorted_result_xy = sorted(screened_result_xy, key=operator.itemgetter(2))  # sort by element 3

with open("output/free_energy.csv", "w+") as output_sorted_xy:
    output_sorted_xy.write('mirna_1'+'\t'+'mirna_2'+'\t'+'free energy/kcal'+'\n')
    print_list_csv(sorted_result_xy, output_sorted_xy)

print('Done.Results saved in output/free_energy.csv')

cluster_distance.cluster_distance_calculation()
plot_list.plot_list(sorted_result_xy)


