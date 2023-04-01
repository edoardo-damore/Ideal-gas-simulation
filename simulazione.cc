#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;


#define kB 1.380649e-23
#define Na 6.02214129e23
#define UMA 1.66054e-27 //kg
#define epsilonArgon 111.84 //K
#define sigmaArgon 362.3e-12 //m
#define numeroMassaArgon 39.948
//sono stati scelti i parametri dell'argon per simulare il gas ideale

long modulo = 2147483647, a = 16807;
double massa = UMA * numeroMassaArgon;
bool interazioneParticelle = false;

void debug (double x)
{
    cout << x << endl;
    return;
}

double distribuzioneCoordinate (double x);

double distribuzioneMaxwellBoltzmann (double temperature, double moduloVelocity);

void metropolis (double &valoreControllo, long &numeroCasualeGrezzo, string oggetto);

void generazionePosizioni (long seed, int N, double lato, double** valoriGenerati);

void generatoreLCG (long &valoreGrezzo, double &valoreNormalizzato);

void generazioneVelocity (long seed, int N, double temperatura, double** velocity);

void forzaLennardJones (int N, int i, double** posizioni, double* forza);

void motoVelocityVerlet(int N, double dt, double lato, double** posizioni, double** velocity);

void differenceCheck (long seed, int N, double lato, double** posizioni);

double funzioneDistanza (int i, int j, double** posizioni);

double moduloQuadro (int N, double** velocity);

double funzionePressione (int N, double volume, double** velocity);

double funzioneTemperatura (int N, double** velocity);

double funzioneEnergia (int N, double** posizioni, double** velocity);

void print (int N, int iterazioneAttuale, string path, double** elementi);

void printRisultati (int N, int iterazioneAttuale, ofstream &fileRisultati, double volume, double** posizioni, double** velocity);


int main()
{
    double lato = 1.e-2; 
    double volume = pow(lato, 3);

    double dt = 1e-6;

    char configurazione;
    cout << "Vuoi che le particelle interagiscano tra loro (y/n)? ";
    cin >> configurazione;

    string cartellaPosizioni, cartellaVelocity, nomeFile;

    if (configurazione == 'y') 
    {
        interazioneParticelle = true;
        cartellaPosizioni = "posizioniInterazione";
        cartellaVelocity = "velocitàInterazione";
        nomeFile = "risultatiInterazione.txt";

    }
    else if (configurazione == 'n') 
    {
        interazioneParticelle = false;
        cartellaPosizioni = "posizioniNoInterazione";
        cartellaVelocity = "velocitàNoInterazione";
        nomeFile = "risultatiNoInterazione.txt";
    }
    else 
    {
        cout << "Invalid input. Program aborted." << endl;
        return 1;
    }

    int N;
    cout << "Inserire il numero di atomi di Argon da generare: ";
    cin >> N;

    double temperatura;
    cout << "Inserire la temperatura (K): ";
    cin >> temperatura;

    int numeroIterazioni;
    cout << "Inserire il numero di iterazioni: ";
    cin >> numeroIterazioni;

    int stepStampa;
    cout << "Inserire ogni quante iterazioni va stampato lo stato del sistema: ";
    cin >> stepStampa;

    int seedPosizioni;
    cout << "Inserire il seed per le posizioni: ";
    cin >> seedPosizioni;

    int seedVelocity;
    cout << "Inserire il seed per le velocità: ";
    cin >> seedVelocity;


    //creazione delle matrici contenenti posizioni e velocità di ogni particella in un dato istante

    double** posizioni = new double*[N];
    for (int i = 0; i < N; i++)
    {
        posizioni[i] = new double[3];
    }

    double** velocity = new double*[N];
    for (int i = 0; i < N; i++)
    {
        velocity[i] = new double[3];
    }

    //generazione casuale delle posizioni e delle velocità delle particelle
    generazionePosizioni(seedPosizioni, N, lato, posizioni);
    generazioneVelocity(seedVelocity, N, temperatura, velocity);

    // controllo non vi siano particelle sovrapposte
    differenceCheck(seedPosizioni, N, lato, posizioni);


    // creazione del file su cui andranno scritti i parametri che ci interessano del sistema
    ofstream fileRisultati (nomeFile);

    fileRisultati << "Iterazione\tEnergia\tPressione (p)\tTemperatura (T)\tp * V\tN * kB * T" << endl;


    int iterazioneAttuale = 0;
    while (iterazioneAttuale < numeroIterazioni)
    {
        // stampa dei risulati ogni tot iterazioni
        if (iterazioneAttuale % stepStampa == 0)
        {
            printRisultati(N, iterazioneAttuale, fileRisultati, volume, posizioni, velocity);
            print(N, iterazioneAttuale, cartellaPosizioni, posizioni);
            print(N, iterazioneAttuale, cartellaVelocity, velocity);
        }

        // utilizzo dell'algoritmo velocity-verlet per l'evoluzione del sistema
        motoVelocityVerlet(N, dt, lato, posizioni, velocity);
 
        iterazioneAttuale++;
    }

    fileRisultati.close();

    // distruzione dei vari array dinamici creati
    for (int i = 0; i < N; i++)
    {
        delete [] posizioni[i];
        delete [] velocity[i];
    }
    delete [] posizioni;
    delete [] velocity;

    return 0;
}

double distribuzioneCoordinate (double x)
{
    return 1;
}

double distribuzioneMetropolis (double x, string oggetto)
{
    if (oggetto=="posizione") return 1;
    else if (oggetto=="theta") return sin(x);
    else if (oggetto=="phi") return 1;
    else return 1;
}

//algoritmo di metropolis che genera un numero casuale secondo una distribuzione scelta in base all'oggetto
void metropolis (double &valoreControllo, long &numeroCasualeGrezzo, string oggetto)
{
    double valoreTrial;

    generatoreLCG(numeroCasualeGrezzo, valoreTrial);

    if (distribuzioneMetropolis(valoreTrial, oggetto) > distribuzioneMetropolis(valoreControllo, oggetto))
    {
        valoreControllo = valoreTrial;
    }
    else 
    {
        double rapporto;
        generatoreLCG(numeroCasualeGrezzo, rapporto);

        if (distribuzioneMetropolis(valoreTrial, oggetto)/distribuzioneMetropolis(valoreControllo, oggetto) > rapporto)
        {
            valoreControllo = valoreTrial;
        }        
    }

    return;
}

//utilizzo di metropolis per generare le posizioni iniziali
void generazionePosizioni (long seed, int N, double lato, double** valoriGenerati)
{
    long numeroCasualeGrezzo = seed;
    double numeroCasualeNormalizzato;

    generatoreLCG(numeroCasualeGrezzo, numeroCasualeNormalizzato);

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            metropolis(numeroCasualeNormalizzato, numeroCasualeGrezzo, "posizione");
            valoriGenerati[i][k] = numeroCasualeNormalizzato * lato;
        }
    }

    return;
}

double distribuzioneMaxwellBoltzmann (double temperatura, double moduloVelocity)
{
    return sqrt(2/M_PI) * pow(massa/(kB * temperatura), 1.5) * pow(moduloVelocity, 2) * exp(- massa * pow(moduloVelocity, 2)/(2 * kB * temperatura));
}

// algoritmo LCG in configurazione di park-miller
void generatoreLCG (long &valoreGrezzo, double &valoreNormalizzato)
{
    valoreGrezzo = (a * valoreGrezzo) % modulo;
    valoreNormalizzato = double(valoreGrezzo)/double(modulo);

    return;
}

void generazioneVelocity (long seed, int N, double temperatura, double** velocity)
{
    double moduloVelocityPeak = sqrt(2 * kB * temperatura / massa); //modulo della velocità con probabilità più alta
    double moduloVelocityMedia = 2/sqrt(M_PI) * moduloVelocityPeak; //modulo della velocità medio
    double moduloVelocitySigma = sqrt((3 * M_PI - 8)/(2 * M_PI)) * moduloVelocityPeak; //deviazione standard

    double massimoDistribuzione = distribuzioneMaxwellBoltzmann(temperatura, moduloVelocityPeak); //valore massimo della distribuzione di M-B

    double rangeModuloVelocity = moduloVelocityMedia + 3 * moduloVelocitySigma; //range in cui generare i valori del modulo della velocità

    long numeroCasualeGrezzo = seed;
    double numeroCasualeNormalizzato;

    double moduloVelocityGenerato, valoreControllo;

    for (int i = 0; i < N; i++)
    {
        //generazione del valore di controllo e del modulo della velocità
        generatoreLCG(numeroCasualeGrezzo, numeroCasualeNormalizzato);
        valoreControllo = numeroCasualeNormalizzato * massimoDistribuzione;

        generatoreLCG(numeroCasualeGrezzo, numeroCasualeNormalizzato);
        moduloVelocityGenerato = numeroCasualeNormalizzato * rangeModuloVelocity;

        //se c > f(v) allora si deve generare un altro modulo
        if (valoreControllo > distribuzioneMaxwellBoltzmann(temperatura, moduloVelocityGenerato))
        {
            i--;
            continue;
        }

        //generare casualmente gli angoli theta e phi secondo la distribuzione sferica
        generatoreLCG(numeroCasualeGrezzo, numeroCasualeNormalizzato);
        metropolis(numeroCasualeNormalizzato, numeroCasualeGrezzo, "theta");
        double theta = numeroCasualeNormalizzato * M_PI;

        generatoreLCG(numeroCasualeGrezzo, numeroCasualeNormalizzato);
        metropolis(numeroCasualeNormalizzato, numeroCasualeGrezzo, "phi");
        double phi = numeroCasualeNormalizzato * 2 * M_PI;

        //trovare le componenti della velocità tramite coordinate sferiche
        velocity[i][0] = moduloVelocityGenerato * sin(theta) * cos(phi);
        velocity[i][1] = moduloVelocityGenerato * sin(theta) * sin(phi);
        velocity[i][2] = moduloVelocityGenerato * cos(theta);
    }

    return;
}

// data la particella i-esima, calcola la forza totale di lennard-jones agente su essa sui tre assi
void forzaLennardJones (int N, int i, double** posizioni, double* forza)
{
    for (int j = 0; j < N; j++)
    {
        if (j==i) continue;

        double distanza = funzioneDistanza(i, j, posizioni);

        for (int k = 0; k < 3; k++)
        {
            forza[k] += 24. * epsilonArgon/pow(distanza, 2) * (2 * pow(sigmaArgon/distanza, 12) - pow(sigmaArgon/distanza, 6)) * (posizioni[i][k] - posizioni[j][k]);
        }
    }

    return;
}

// propagatore velocity-verlet utilizzato per far evolvere il sistema
void motoVelocityVerlet(int N, double dt, double lato, double** posizioni, double** velocity)
{
    //posizioni e velocità all'istante n+1
    double** posizioniFuture = new double*[N];
    for (int i = 0; i < N; i++)
    {
        posizioniFuture[i] = new double[3];
    }

    double** velocityFuture = new double*[N];
    for (int i = 0; i < N; i++)
    {
        velocityFuture[i] = new double[3];
    }


    double* forzaLJ = new double[3];

    for (int i = 0; i < N; i++) //generazione delle posizioni all'instante di tempo n+1
    {
        for (int k = 0; k < 3; k++)
        {
            forzaLJ[k] = 0;
        }

        if (interazioneParticelle) forzaLennardJones(N, i, posizioni, forzaLJ);

        for (int k = 0; k < 3; k++)
        {
            posizioniFuture[i][k] = posizioni[i][k] + velocity[i][k] * dt + 1/(2 * massa) * forzaLJ[k] * pow(dt, 2); // + O(dt^3)
        }

    }

    double* forzaFuturaLJ = new double[3]; 
    
    for (int i = 0; i < N; i++) //generazione delle velocità all'istante di tempo n+1
    {
        for (int k = 0; k < 3; k++)
        {
            forzaLJ[k] = 0;
            forzaFuturaLJ[k] = 0;
        }

        if (interazioneParticelle)
        {
            forzaLennardJones(N, i, posizioni, forzaLJ);
            forzaLennardJones(N, i, posizioniFuture, forzaFuturaLJ);
        }

        for (int k = 0; k < 3; k++)
        {
            velocityFuture[i][k] = velocity[i][k] + 1/(2 * massa) * (forzaLJ[k] + forzaFuturaLJ[k]) * dt; // + O(dt^2)

            //se la particella ha superato i bordi del cubo, la velocità è invertita
            if ((posizioniFuture[i][k] <= 0) or (posizioniFuture[i][k] >= lato)) velocityFuture[i][k] = - velocityFuture[i][k];
        }

    }

    //poniamo come presente l'istante n+1
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            posizioni[i][k] = posizioniFuture[i][k];
            velocity[i][k] = velocityFuture[i][k];
        }
    }

    for (int i = 0; i < N; i++)
    {
        delete [] posizioniFuture[i];
        delete [] velocityFuture[i];
    }
    delete [] posizioniFuture;
    delete [] velocityFuture;
    delete [] forzaLJ;
    delete [] forzaFuturaLJ;

   return;
}

// controllo se due particelle sono nello stesso punto e in quel caso sposto la seconda
void differenceCheck (long seed, int N, double lato, double** posizioni)
{
    double** posizioniSostitutive = new double*[1];
    posizioniSostitutive[0] = new double[3];

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j==i) continue;

            if (funzioneDistanza(i, j, posizioni) == 0)
            {
                seed++;

                generazionePosizioni(seed, 1, lato, posizioniSostitutive);

                for (int k = 0; k < 3; k++)
                {
                    posizioni[j][k] = posizioniSostitutive[0][k];
                }

                j--;
            }
        }
    }

    delete [] posizioniSostitutive[0];
    delete [] posizioniSostitutive;

    return;
}

double funzioneDistanza (int i, int j, double** posizioni)
{
    double distanza = 0;
    for (int k = 0; k < 3; k++)
    {
        distanza += pow(posizioni[i][k] - posizioni[j][k], 2);
    }

    return sqrt(distanza);
}

double moduloQuadro (int N, double** velocity)
{
    double sommaQuadra = 0;

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            sommaQuadra += pow(velocity[i][k], 2);
        }
    }

    return sommaQuadra;
}

double funzionePressione (int N, double volume, double** velocity)
{
    double moduloVelocityQuadro = moduloQuadro(N, velocity);

    return 1/(3 * volume) * massa * moduloVelocityQuadro;
}

double funzioneTemperatura (int N, double** velocity)
{
    double moduloVelocityQuadro = moduloQuadro(N, velocity);

    return 1/((3 * N - 3) * kB) * massa * moduloVelocityQuadro;
}

double funzioneEnergia (int N, double** posizioni, double** velocity)
{
    double moduloVelocityQuadro = moduloQuadro(N, velocity);

    double energiaCinetica = 0.5 * massa * moduloVelocityQuadro;


    double potenzialeLennardJones = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j==i) continue;

            double distanza = funzioneDistanza(i, j, posizioni);

            potenzialeLennardJones += 0.5 * 4. * epsilonArgon * (pow(sigmaArgon/distanza, 12) - pow(sigmaArgon/distanza, 6));
        }
    }
    if (!interazioneParticelle) potenzialeLennardJones = 0;

    double energiaTotale = energiaCinetica + potenzialeLennardJones;
    
    return energiaTotale;
}

// scrive un file con posizioni o velocità in un dato momento nella cartella scelta
void print (int N, int iterazioneAttuale, string path, double** elementi)
{
    ofstream file (path + "/" + to_string(iterazioneAttuale) + ".txt");

    file << "X\tY\tZ" << endl;
    
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            file << elementi[i][k];

            if (k==2) continue;
            file << "\t";
        }

        file << endl;
    }

    file.close();   
}

// stampa su un file i parametri principali del sistema in un determinato istante
void printRisultati (int N, int iterazioneAttuale, ofstream &fileRisultati, double volume, double** posizioni, double** velocity)
{
    double energiaTotaleIstantanea = funzioneEnergia(N, posizioni, velocity);
    double pressioneIstantanea  = funzionePressione(N, volume, velocity);
    double temperaturaIstantanea = funzioneTemperatura(N, velocity);

    fileRisultati << iterazioneAttuale << "\t" << energiaTotaleIstantanea << "\t" << pressioneIstantanea << "\t" << temperaturaIstantanea << "\t" << pressioneIstantanea * volume << "\t" << N * kB * temperaturaIstantanea << endl;

}