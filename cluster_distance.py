__author__ = 'chobits'
import csv
import re
import plot_list
"""This script:
*Input format:cluster.csv,free_energy.csv
*Output format:results.txt(xls format)

*from a csv file with the structure [mirna, mirna, free energy]
calculate the cluster distance with an output format as [mirna, mirna, free energy, distance]

beta:2015-5-19
This line was left blank intentionally."""


def cluster_distance_calculation():
    print('Calculating cluster distance...')
    with open("database/cluster.csv", "r+") as database:
        data_csv = csv.reader(database)
        data_list = list(data_csv)
    with open("output/free_energy.csv", "r+") as query:
        query_csv = csv.reader(query)
        query_list = list(query_csv)

    result = []
    for i in range(1, len(query_list)):
        for a in range(1, len(data_list)):
            flag_a = re.search(query_list[i][0], data_list[a][0])
            if flag_a:
                break
        for b in range(1, len(data_list)):
            flag_b = re.search(query_list[i][1], data_list[a][0])
            if flag_b:
                break

        if flag_a and flag_b:
            cluster_a_search = re.search('^\d*', data_list[a][3])
            cluster_b_search = re.search('^\d*', data_list[b][3])
            cluster_a = cluster_a_search.group()
            cluster_b = cluster_b_search.group()
            if data_list[a][1] == data_list[b][1]:
                result.append([query_list[i][0],
                               query_list[i][1],
                               query_list[i][2],
                               abs(float(cluster_a)-float(cluster_b))]
                              )
            else:
                result.append([query_list[i][0],
                               query_list[i][1],
                               query_list[i][2],
                               2.0e8]
                              )
    plot_list.plot_bar(result)
    plot_list.plot_scatter(result)

    with open("output/results.txt", "w+") as output:
        output.write('mirna_1'+'\t'+'mirna_2'+'\t'+'free energy/kcal'+'\t'+'cluster distance'+'\n')
        for i in range(len(result)):
            output.write(str(result[i][0])+'\t'+result[i][1]+'\t'+str(result[i][2])+'\t'+str(result[i][3])+'\n')
    print('Done.Results saved in output/results.txt')
    return 0