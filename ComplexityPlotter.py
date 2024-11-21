# keys : Number of boids are 5, 20, 50, 100, 200, 500, 2000, 5000
# values : Time taken to compute the flocking behavior
# values that are arrays stretch from 100, 1000, 10000, 100000, 1000000, 10000000
# all times are in seconds


import matplotlib.pyplot as plt

cellListSerial = {
    5 : [0.000431061, 0.0012629, 0.0145791, 0.131556, 1.39329, 13.4156]
    , 20 : [0.000876904, 0.011167, 0.0767398, 0.964271, 9.76051, 107.803]
    , 50 : [0.00351906, 0.0331268, 0.399619, 3.89028, 41.4308]
    , 100 : [0.010392, 0.109823, 1.36479, 11.6787, 122.553]
    , 200 : [0.0320132, 0.41608, 4.3159, 42.3148, 479.551]
    , 500 : [0.153826, 2.89095, 32.9393, 341.661, 3161.09]
    , 2000 : [4.10836, 47.6298, 557.841, 5506.09]
    , 5000 : [19.0703, 382.212, 3759.29]
}

nonCellListSerial = {
    5 : [0.000264168, 0.000865936, 0.0085299, 0.0870869, 0.866721, 8.69318]
    , 20 : [0.000760078, 0.00455403, 0.044435, 0.439804, 4.40948, 44.3134]
    , 50 : [0.00300407, 0.0277591, 0.276082, 2.75873, 27.6924, 278.176]
    , 100 : [0.010901, 0.104441, 1.04284, 10.4412, 105.001, 1052.53]
    , 200 : [0.0416071, 0.410179, 4.10404, 41.114, 413.672]
    , 500 : [0.253247, 2.52217, 25.2543, 253.568, 2548.02]
    , 2000 : [2.65505, 51.4714, 544.083]
    , 5000 : [24.3355, 258.982, 2440.89]
}

initialParallel = {
    5 : [0.0134799, 0.0322652, 0.178481, 0.14241],
    20 : [0.028856, 0.154878, 0.115671, 0.369385],
    50 : [0.0284851, 0.0637732, 0.169308, 0.203225],
    100 : [0.047617, 0.129852, 0.231031, 0.399773],
    200 : [0.38805, 0.119245, 0.388084,  0.944709],
    500 : [0.73525, 0.528363, 0.900905, 0.651199],
    2000 : [9.86751, 5.22129, 3.52873, 1.99147],
    5000 : [63.0375, 38.2399, 16.1931, 11.4019]
}

def compareSets(set1, set2):

    x_values = [100, 1000, 10000, 100000, 1000000, 10000000]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange']

    for idx, key in enumerate(set1):
        y_values = set1[key]
        plt.plot(x_values[:len(y_values)], y_values, label=f'cellListSerial {key}', linestyle='-', marker='o', color=colors[idx])

    for idx, key in enumerate(set2):
        y_values = set2[key]
        plt.plot(x_values[:len(y_values)], y_values, label=f'nonCellListSerial {key}', linestyle='--', marker='x', color=colors[idx])

    plt.xlabel('Number of iterations')
    plt.ylabel('Time (seconds)')
    plt.title('Time taken to compute the flocking behavior')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    # plt.legend()
    plt.show()

compareSets(cellListSerial, nonCellListSerial)


def plotInitialParallel(data):
    x_values = [4, 8, 16, 28]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange']

    for idx, key in enumerate(data):
        y_values = data[key]
        plt.plot(x_values[:len(y_values)], y_values, label=f'initialParallel {key}', linestyle='-', marker='o', color=colors[idx])

    plt.xlabel('Number of threads')
    plt.ylabel('Time (seconds)')
    plt.title('Time taken to compute the flocking behavior with initial parallelization')
    # plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

plotInitialParallel(initialParallel)
