from matplotlib import pyplot as plt
import matplotlib.animation as animation
import polars as pl
from mpl_toolkits.mplot3d import Axes3D

file_path = '/user/home/oz21652/Comp401/BoidSimulation.csv'
dfs = []

columnDict = {
    0 : 'posX',
    1 : 'posY',
    2 : 'posZ',
    3 : 'dirX',
    4 : 'dirY',
    5 : 'dirZ',
    6 : 'speed',
    7 : 'mass'
}

with open(file_path, 'r') as file:
    for line in file:
        dataDict = {
            'posX': [],
            'posY': [],
            'posZ': [],
            'dirX': [],
            'dirY': [],
            'dirZ': [],
            'speed': [],
            'mass': []
            }
            
        line = line.strip()

        for boid in line.split('],'):
            boid = boid.replace('[', '')

            for index, value in enumerate(boid.split(',')):
                if value == '':
                    continue

                dataDict[columnDict[index]].append(float(value))

        dfs.append(pl.DataFrame(dataDict))

figure = plt.figure()

ax = figure.add_subplot(111, projection='3d')


allXVals = [val for df in dfs for val in df['posX'].to_list()]
allYVals = [val for df in dfs for val in df['posY'].to_list()]
allZVals = [val for df in dfs for val in df['posZ'].to_list()]

def update(frame):
    ax.clear()
    scat = ax.scatter(dfs[frame]['posX'], dfs[frame]['posY'], dfs[frame]['posZ'], c='b', marker='o')
    ax.set_xlim(min(allXVals), max(allXVals))
    ax.set_ylim(min(allYVals), max(allYVals))
    ax.set_zlim(min(allZVals), max(allZVals))
    ax.set_title('Boid Simulation')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.grid(True)

    return [scat]


ani = animation.FuncAnimation(figure, update, frames=len(dfs), repeat=False, blit=True)

ani.save('simulation.gif', writer='pillow', fps=100)
plt.show()
