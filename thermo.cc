#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;


#define G 6.67e-11
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

void metropolisPosizioni (long seed, int N, double lato, double** valoriGenerati);

void generatoreLCG (long &valoreGrezzo, double &valoreNormalizzato);

void generazioneVelocity (long seed, int N, double temperatura, double** velocity);

void forzaLennardJones (int N, int i, double** posizioni, double* forza);

void motoVelocityVerlet(int N, double dt, double lato, double** posizioni, double** velocity);

void differenceCheck (long seed, int N, double lato, double** posizioni);

double funzioneDistanza (int N, int i, int j, double** posizioni);

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

    if (configurazione == 'y') interazioneParticelle = true;
    else if (configurazione == 'n') interazioneParticelle = false;
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
    metropolisPosizioni(seedPosizioni, N, lato, posizioni);
    generazioneVelocity(seedVelocity, N, temperatura, velocity);

    // controllo non vi siano particelle sovrapposte
    differenceCheck(seedPosizioni, N, lato, posizioni);


    // creazione del file su cui andranno scritti i parametri che ci interessano del sistema
    ofstream fileRisultati ("risultati.txt");

    fileRisultati << "Iterazione\tEnergia\tPressione (p)\tTemperatura (T)\tp * V\tN * kB * T" << endl;


    int iterazioneAttuale = 0;
    while (iterazioneAttuale < numeroIterazioni)
    {
        // stampa dei risulati ogni tot iterazioni
        if (iterazioneAttuale % stepStampa == 0)
        {
            printRisultati(N, iterazioneAttuale, fileRisultati, volume, posizioni, velocity);
            print(N, iterazioneAttuale, "posizioni", posizioni);
            print(N, iterazioneAttuale, "velocità", velocity);
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

// utilizzo dell'algoritmo di metropolis per la creazione delle posizioni delle particelle
// i numeri casuali vengono generati tramite LCG con parametri di park-miller partendo da un seed scelto a caso dall'utente
void metropolisPosizioni (long seed, int N, double lato, double** valoriGenerati)
{

    double valoreControllo;
    double valoreTrial;

    long numeroCasualeGrezzo = (a * seed) % modulo;
    double numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.) / (modulo * 1.);

    valoreControllo = numeroCasualeNormalizzato;

    for (int k = 0; k < 3; k++)
    {
        for (int i = 0; i < N; i++)
        {
            numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
            numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

            valoreTrial = numeroCasualeNormalizzato;

            if (distribuzioneCoordinate(valoreTrial) > distribuzioneCoordinate(valoreControllo))
            {
                valoreControllo = valoreTrial;
                valoriGenerati[i][k] = valoreTrial * lato;
            }
            else 
            {
                numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
                numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

                double rapporto = numeroCasualeNormalizzato;

                if (distribuzioneCoordinate(valoreTrial)/distribuzioneCoordinate(valoreControllo) > rapporto)
                {
                    valoreControllo = valoreTrial;
                    valoriGenerati[i][k] = valoreTrial * lato;
                }
                else
                {
                    valoriGenerati[i][k] = valoreControllo * lato;
                }
            }
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

// generazione velocità di partenza secondo la distribuzione di maxwell-boltzmann
// i numeri casuali sono generati tramite LCG
void generazioneVelocity (long seed, int N, double temperatura, double** velocity)
{

    double velocityPeak = sqrt(2 * kB * temperatura / massa); //velocità con probabilità più alta
    double massimoDistribuzione = distribuzioneMaxwellBoltzmann(temperatura, velocityPeak); 

    double velocityMedia = 2/sqrt(M_PI) * velocityPeak; //velocità media
    //double velocityQuadroMedia = 1.5 * pow(velocityPeak, 2);

    double velocitySigma = sqrt((3 * M_PI - 8)/(2 * M_PI)) * velocityPeak; //deviazione standard

    double rangeVelocity = velocityMedia + 3 * velocitySigma;

    long valoreGeneratoGrezzo = seed;
    double valoreGeneratoNormalizzato;
    double valoreGeneratoFinale;
    double valoreGeneratoControllo;

    double* velocityParticellaCorrente = new double[3];

    for (int i = 0; i < N; i++) 
    {
        for (int k = 0; k < 3; k++) //generazione delle velocità lungo i tre assi per una particella
        {
            generatoreLCG(valoreGeneratoGrezzo, valoreGeneratoNormalizzato);
            valoreGeneratoFinale = (valoreGeneratoNormalizzato * 2 - 1) * rangeVelocity / sqrt(3);
            velocityParticellaCorrente[k] = valoreGeneratoFinale;
        }

        double moduloVelocity = sqrt(pow(velocityParticellaCorrente[0],2) + 
                                     pow(velocityParticellaCorrente[1],2) +
                                     pow(velocityParticellaCorrente[2],2));
    
        generatoreLCG(valoreGeneratoGrezzo, valoreGeneratoNormalizzato);
        valoreGeneratoControllo = valoreGeneratoNormalizzato * massimoDistribuzione;

        if (valoreGeneratoControllo <= distribuzioneMaxwellBoltzmann(temperatura, moduloVelocity)) //check 
        {
            for (int k = 0; k < 3; k++)
            {
                velocity[i][k] = velocityParticellaCorrente[k];
            }
        }        
        else i--;
    }

    delete [] velocityParticellaCorrente;

    return;
}

// data la particella i-esima, calcola la forza totale di lennard-jones agente su essa sui tre assi
void forzaLennardJones (int N, int i, double** posizioni, double* forza)
{

    for (int j = 0; j < N; j++)
    {
        if (j==i) continue;

        double distanza = funzioneDistanza(N, i, j, posizioni);

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

            if (funzioneDistanza(N, i, j, posizioni) == 0)
            {
                seed++;

                metropolisPosizioni(seed, 1, lato, posizioniSostitutive);

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

double funzioneDistanza (int N, int i, int j, double** posizioni)
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

            double distanza = funzioneDistanza(N, i, j, posizioni);

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