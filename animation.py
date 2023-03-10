import numpy as np
import pandas as pd
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt


print("Quante iterazioni sono state fatte? ")
numeroIterazioni = int(input())

print("Quante particelle ha il sistema? ")
numeroParticelle = int(input())


#matrice tridimensionale contenente tutte le coordinate di tutte le particelle in ogni istante di tempo
iterations = np.empty(shape=[numeroIterazioni, numeroParticelle,  3]) 

for i in range(numeroIterazioni):
    
    istante = pd.read_table(f"iterazioni/{i+1}.txt")

    for j in range(numeroParticelle):

        iterations[i][j][0] = istante["x"][j]
        iterations[i][j][1] = istante["y"][j]
        iterations[i][j][2] = istante["z"][j]


fig = plt.figure()
ax = fig.add_subplot(projection="3d")


def animate(i):

    plt.cla()

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.set_xlim(0, 0.01)
    ax.set_ylim(0, 0.01)
    ax.set_zlim(0, 0.01)

    for j in range(numeroParticelle):
        ax.scatter(iterations[i][j][0], iterations[i][j][1], iterations[i][j][2], color='blue')


anim = FuncAnimation(fig, animate, interval=1000, repeat=False)

writergif = animation.PillowWriter(fps=5)

anim.save("animation.gif", writer=writergif)
