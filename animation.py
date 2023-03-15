import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print("Inserire il numero di particelle: ")
numeroParticelle = int(input())

print("Inserire il numero di iterazioni: ")
numeroIterazioni = int(input())

outputInterazione = pd.read_table("outputInterazione.txt")
outputNoInterazione = pd.read_table("outputNoInterazione.txt")

iterazione = outputInterazione["Iterazione"]

energiaInterazione = outputInterazione["Energia"]
pvInterazione = outputInterazione["pV"]
nkbtInterazione = outputInterazione["NkBT"]

energiaNoInterazione = outputNoInterazione["Energia"]
pvNoInterazione = outputNoInterazione["pV"]
nkbtNoInterazione = outputNoInterazione["NkBT"]

massimoEnergia = 0

for e in energiaInterazione:
    if (e>massimoEnergia):
        massimoEnergia = e

minimoEnergia = massimoEnergia

for e in energiaInterazione:
    if (e < minimoEnergia):
        minimoEnergia = e

print(energiaInterazione)


fig, axs = plt.subplots(2)
fig.set_size_inches(12, 10)
fig.suptitle("Andamento dell'energia totale del sistema")

axs[0].set_xlabel("Iterazioni")
axs[0].set_ylabel(r"Energia (J)")

axs[0].set_xlim(0, numeroIterazioni)
axs[1].set_xlim(0, numeroIterazioni)
#axs[1].set_ylim(minimoEnergia, massimoEnergia)


axs[1].set_xlabel("Iterazioni")
axs[1].set_ylabel(r"Energia (J)")


axs[0].set_title("Senza interazione tra particelle")
axs[1].set_title("Con interazione tra particelle")

axs[0].plot(iterazione, energiaNoInterazione)
axs[1].plot(iterazione, energiaInterazione)

fig.savefig("comparazioneEnergia.png")

plt.cla()

fig, axs = plt.subplots(2)
fig.set_size_inches(12, 10)
fig.suptitle("Andamento legge dei gas ideali")

axs[0].set_xlabel("Iterazioni")
axs[0].set_ylabel(r"$pV \longleftrightarrow N k_B T$")

axs[1].set_xlabel("Iterazioni")
axs[1].set_ylabel(r"$pV \longleftrightarrow N k_B T$")

axs[0].set_title("Senza interazione tra particelle")
axs[1].set_title("Con interazione tra particelle")

axs[0].set_xlim(0, numeroIterazioni)
axs[1].set_xlim(0, numeroIterazioni)

axs[0].plot(iterazione, pvNoInterazione, label=r"$pV$")
axs[0].plot(iterazione, nkbtNoInterazione, label=r"$N k_b T$")

axs[1].plot(iterazione, pvInterazione, label=r"$pV$")
axs[1].plot(iterazione, nkbtInterazione, label=r"$N k_B T$")

axs[0].legend()

axs[1].legend()

fig.savefig("comparazioneLegge.png")

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,9), subplot_kw={"projection": "3d"})

lato = 0.01

def animate(i):

    axs[0].clear()
    axs[1].clear()

    fig.suptitle("Visualizzazione del sistema")

    axs[0].set_title("Senza interazione")
    axs[1].set_title("Con interazione")

    axs[0].set_xlabel("X")
    axs[0].set_ylabel("Y")
    axs[0].set_zlabel("Z")

    axs[1].set_xlabel("X")
    axs[1].set_ylabel("Y")
    axs[1].set_zlabel("Z")

    axs[0].set_xlim(0, lato)
    axs[0].set_ylim(0, lato)
    axs[0].set_zlim(0, lato)

    axs[1].set_xlim(0, lato)
    axs[1].set_ylim(0, lato)
    axs[1].set_zlim(0, lato)



    cI = pd.read_table(f".coordinateInterazione/{i}.txt")
    cNI = pd.read_table(f".coordinateNoInterazione/{i}.txt")
    vI = pd.read_table(f".velocitàInterazione/{i}.txt")
    vNI = pd.read_table(f".velocitàNoInterazione/{i}.txt")

    axs[0].scatter(cNI["x"], cNI["y"], cNI["z"], color="blue")
    axs[1].scatter(cI["x"], cI["y"], cI["z"], color="blue")   

anim = animation.FuncAnimation(fig, animate, numeroIterazioni, blit=False)
anim.save("comparazioneParticelle.gif", writer=animation.PillowWriter(fps=5))