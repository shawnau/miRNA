__author__ = 'chobits'
import numpy as np
import matplotlib.pyplot as plt
import math


def plot_list(list_line):
    n = len(list_line)
    x = np.arange(n)
    y = []
    for i in range(n):
        y.append(list_line[i][2])
    plt.plot(x, y, color="blue", linewidth=1.5, linestyle="-")

    plt.xlabel('miRNA pair')
    plt.ylabel('del G/kcal')
    plt.savefig("output/free_energy.png", dpi=72)
    print('*Figure saved to output/free_energy.png')
    plt.show()


def plot_scatter(lst):
    x = []
    y = []
    for i in range(len(lst)):
        x.append(-float(lst[i][2]))
        y.append(float(lst[i][3]))
    plt.plot(x, y, 'o', color="blue")
    plt.xlabel('-G/kcal')
    plt.ylabel('distance')
    plt.savefig("output/distance_energy.png", dpi=72)
    print('*Figure saved to output/distance_energy')
    plt.show()


def plot_bar(lst):
    max_number = -int(float(lst[0][2]))+1
    count_chr= [0 for x in range(max_number)]
    count_nonchr = [0 for x in range(max_number)]
    average = []
    section = [x for x in range(max_number)]
    results = [[] for x in range(max_number)]

    for i in range(len(lst)):
        results[-int(float(lst[i][2]))].append(lst[i][3])

    for i in range(max_number):
        for j in range(len(results[i])):
            if results[i][j] == 2.0e8:
                count_nonchr[i] += 1
            else:
                count_chr[i] += 1

    for i in range(len(results)):
        if len(results[i]) != 0:
            average.append(sum(results[i])/len(results[i]))
        else:
            average.append(0)

    p1 = plt.bar(section, count_chr, color='r')
    p2 = plt.bar(section, count_nonchr, color='y', bottom=count_chr)

    plt.xlabel('-G/kcal')
    plt.ylabel('amount')
    plt.legend((p1[0], p2[0]), ('same chr', 'diff chr'))
    plt.savefig("output/bar.png", dpi=72)
    #plt.plot(section, average)
    #plt.savefig("output/line.png", dpi=72)
    plt.show()