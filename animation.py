import numpy as np
import pandas as pd
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import PIL.Image as image
import array

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

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")


def animate(i):

    for j in range(numeroParticelle):
        ax.scatter(iterations[i][j][0], iterations[i][j][1], iterations[i][j][2], color='blue')


anim = FuncAnimation(fig, animate, interval=100, repeat=False)

writergif = animation.PillowWriter(fps=30)

anim.save("animation.gif", writer=writergif)




#for j in range(numeroIterazioni):
#    for i in range(numeroParticelle):
#
#        ax.scatter(iterations[j][i][0], iterations[j][i][1], iterations[j][i][2], color='blue')
#
#    fig.savefig(f"frames/{j}.png", dpi='figure', format='png')
#
#images = []
#
#for i in range(numeroIterazioni):
#    images.append(image.open(f"frames/{i}.png"))
#
#images[0].save("prova.gif", save_all=True, append_images=images[1:], optimize=False, duration=numeroIterazioni, loop=0)
