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

initialSerial = {
    5 : [0.000865936, 0.00135824, 0.000953965],
    20 : [0.00455403, 0.00297878, 0.00375438],
    50 : [0.0277591, 0.0126213, ],
    100 : [0.104441, 0.0460696, 0.0461556],
    200 : [0.410179, 0.170698, 0.170372],
    500 : [2.52217, 0.985575, 0.0125752],
    2000 : [51.4714, 14.6849, 14.6938],
    5000 : [258.982, 87.366, 87.3052]
}

initialParallel4 = {
    5: [0.0134799, 0.0205304, 0.0155911],
    20: [0.028856, 0.0187193, 0.024578],
    50: [0.0284851, 0.036857, 0.0269392],
    100: [0.047617, 0.0509305, 0.0594744],
    200: [0.38805, 0.11171, 0.160339],
    500: [0.73525,  0.509405, 0.500131],
    2000: [9.86751, 4.44451, 4.46593],
    5000: [63.0375, 28.1508, 28.1095]
}

initialParallel8 = {
    5: [0.0322652, 0.374933, 0.202163],
    20: [0.154878, 0.477128, 0.143093] ,
    50: [0.0637732, 0.0490083, 0.0483636],
    100: [0.129852, 0.254739, 0.259544],
    200: [0.119245, 0.312193, 0.333327],
    500: [0.528363, 0.530821, 0.868392],
    2000: [5.22129, 2.51511, 2.51267],
    5000: [38.2399, 17.0975, 14.1518]
}

initialParallel16 = {
    5: [0.178481, 0.0707547, 0.464457, ],
    20: [0.115671, 0.728587, 0.279131],
    50: [0.169308, 0.160313, 0.190245, ],
    100: [0.231031, 0.876578, 0.919611],
    200: [0.388084, 0.805043, 0.984688],
    500: [0.900905, 1.33905, 1.09713, ],
    2000: [3.52873, 2.02753, 2.02092],
    5000: [16.1931, 8.6953, 8.63836]
}

initialParallel28 = {
    5: [0.14241, 0.425157, 0.303253],
    20: [0.369385, 0.218089, 0.449517],
    50: [0.203225, 0.293378, 0.39907],
    100: [0.399773, 0.59982, 0.835831],
    200: [0.944709, 1.48291, 1.36603],
    500: [0.651199, 1.61528, 1.65186],
    2000: [1.99147, 1.96014, 2.00924],
    5000: [11.4019, 10.0178, 9.90301]
}


mpiInitial8 = { 
    5 : [0.301409, 0.294897, 0.195516],
    20 : [0.30144, 0.518421, 0.459304],
    50 : [0.376975, 0.445665, 0.440459],
    100 : [0.757991, 0.319718, 0.50314],
    200 : [0.48591, 0.387078, 0.443098],
    500 : [0.548027, 0.969682, 0.596153],
    2000 : [1.28303, 1.33233, 1.20081],
    5000 : [2.84257, 2.78449, 2.75487]
}

mpiInitial4 = {
    5 : [0.134424, 0.362292, 0.456371],
    20 : [0.206021, 0.383926, 0.556805],
    50 : [0.616278, 0.511937, 0.430955],
    100 : [0.426967, 0.45179, 0.460359],
    200 : [0.518171, 0.715657, 0.831657],
    500 : [0.908089, 0.767825, 0.767372],
    2000 : [1.08876, 1.89128, 1.85191],
    5000 : [4.008, 4.2853, 4.22755]
}

mpiAsync4 = {
    5 : [0.23032, 0.363753, ],
    20 : [0.277509, 0.528361, ],
    50 : [0.340986, 0.278243, ],
    100 : [0.737634, 0.679606, ],
    200 : [1.09002, 1.06094, ],
    500 : [1.4798, 1.70892, ],
    2000 : [2.05099, 1.86242, ],
    5000 : [6.76111, 6.78084, ]
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
