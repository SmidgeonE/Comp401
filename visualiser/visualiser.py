from matplotlib import pyplot as plt
import matplotlib.animation as animation
import polars as pl

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

figure, ax = plt.subplots()

allXVals = [val for df in dfs for val in df['posX'].to_list()]
allYVals = [val for df in dfs for val in df['posY'].to_list()]


def update(frame):
    ax.clear()
    asd = ax.plot(dfs[frame]['posX'], dfs[frame]['posY'], 'bo')
    ax.set_xlim(min(allXVals), max(allXVals))
    ax.set_ylim(min(allYVals), max(allYVals))
    ax.set_title('Boid Simulation')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.grid(True)

    return asd


ani = animation.FuncAnimation(figure, update, frames=len(dfs), repeat=False, blit=True)

ani.save('simulation.gif', writer='pillow', fps=10)
plt.show()
