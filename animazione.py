import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print("Inserire il numero di particelle: ")
numeroParticelle = int(input())
print("Inserire il numero di iterazioni: ")
numeroIterazioni = int(input())

#import dei file con i risultati per entrambe le casistiche
outputInterazione = pd.read_table("risultatiInterazione.txt")
outputNoInterazione = pd.read_table("risultatiNoInterazione.txt")

#estrazione dei dati dai dataframes
iterazioni = outputInterazione["Iterazione"]

energiaI = outputInterazione["Energia"]
sinistraI = outputInterazione["p * V"]
destraI = outputInterazione["N * kB * T"]

energiaN = outputNoInterazione["Energia"]
sinistraN = outputNoInterazione["p * V"]
destraN = outputNoInterazione["N * kB * T"]


#generazione primo plot con andamento energia totale del sistema
fig, axs = plt.subplots(2)
fig.set_size_inches(12, 9)
fig.suptitle("Andamento dell'energia totale del sistema")

axs[0].set_xlabel("Iterazioni")
axs[0].set_ylabel(r"Energia (J)")

axs[1].set_xlabel("Iterazioni")
axs[1].set_ylabel(r"Energia (J)")

axs[0].set_title("Senza interazione tra particelle")
axs[1].set_title("Con interazione tra particelle")

axs[0].set_xlim(0, numeroIterazioni)
axs[1].set_xlim(0, numeroIterazioni)

axs[0].plot(iterazioni, energiaN)
axs[1].plot(iterazioni, energiaI)

fig.savefig("comparazioneEnergia.png")

plt.cla()

#generazione secondo plot con visualizzazione dell'aderimento alla legge dei gas ideali
fig, axs = plt.subplots(2)
fig.set_size_inches(12, 9)
fig.suptitle("Andamento legge dei gas ideali")

axs[0].set_xlabel("Iterazioni")
axs[0].set_ylabel(r"$pV \longleftrightarrow N k_B T$")

axs[1].set_xlabel("Iterazioni")
axs[1].set_ylabel(r"$pV \longleftrightarrow N k_B T$")

axs[0].set_title("Senza interazione tra particelle")
axs[1].set_title("Con interazione tra particelle")

axs[0].set_xlim(0, numeroIterazioni)
axs[1].set_xlim(0, numeroIterazioni)

axs[0].plot(iterazioni, sinistraN, label=r"$pV$")
axs[0].plot(iterazioni, destraN, label=r"$N k_b T$")

axs[1].plot(iterazioni, sinistraI, label=r"$pV$")
axs[1].plot(iterazioni, destraI, label=r"$N k_B T$")

axs[0].legend()
axs[1].legend()

fig.savefig("comparazioneLegge.png")

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,9), subplot_kw={"projection": "3d"})

lato = 1e-6

def animate(i):

    axs[0].clear()
    axs[1].clear()

    fig.suptitle("Visualizzazione del sistema")

    axs[0].set_title("Senza interazione")
    axs[1].set_title("Con interazione")

    axs[0].set_xlabel("X (m)")
    axs[0].set_ylabel("Y (m)")
    axs[0].set_zlabel("Z (m)")

    axs[1].set_xlabel("X (m)")
    axs[1].set_ylabel("Y (m)")
    axs[1].set_zlabel("Z (m)")

    axs[0].set_xlim(0, lato)
    axs[0].set_ylim(0, lato)
    axs[0].set_zlim(0, lato)

    axs[1].set_xlim(0, lato)
    axs[1].set_ylim(0, lato)
    axs[1].set_zlim(0, lato)



    cI = pd.read_table(f"posizioniInterazione/{iterazioni[i]}.txt")
    cNI = pd.read_table(f"posizioniNoInterazione/{iterazioni[i]}.txt")
    #vI = pd.read_table(f"velocitàInterazione/{iterazioni[i]}.txt")
    #vNI = pd.read_table(f"velocitàNoInterazione/{iterazioni[i]}.txt")

    axs[0].scatter(cNI["X"], cNI["Y"], cNI["Z"], color="blue")
    axs[1].scatter(cI["X"], cI["Y"], cI["Z"], color="blue")   

anim = animation.FuncAnimation(fig, animate, len(iterazioni), blit=False)
anim.save("comparazioneParticelle.gif", writer=animation.PillowWriter(fps=5))