#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


#define G 6.67e-11
#define kB 1.380649e-23
#define Na 6.02214129e23
#define UMA 1.66054e-27 //kg
#define epsilonArgon 111.84 //K
#define sigmaArgon 362.3e-12 //m
#define numeroMassaArgon 39.948
//sono stati scelti i parametri dell'argon per simulare il gas ideale

bool interazioneParticelle = false;

void debug(double prova)
{
    cout << prova << endl;
    return;
}

double distribuzioneCoordinate (double x)
{
    return 1;
}

double distribuzioneVelocity (double v)
{
    return 1;
}

//algoritmo di metropolis generante N punti casuali nel range [offset, offset + range] 
void metropolisCoordinate (long seed, long modulo, long a, int N, 
                           double offset, double range, double* valoriGenerati)
{

    double valoreControllo;
    double valoreTrial;

    long numeroCasualeGrezzo = (a * seed) % modulo;
    double numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.) / (modulo * 1.);

    valoreControllo = numeroCasualeNormalizzato;


    for (int i = 0; i < N; i++)
    {
        numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
        numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

        valoreTrial = numeroCasualeNormalizzato;

        if (distribuzioneCoordinate(valoreTrial) > distribuzioneCoordinate(valoreControllo))
        {
            valoreControllo = valoreTrial;
            valoriGenerati[i] = valoreTrial * range + offset;
        }
        else 
        {
            numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
            numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

            double rapporto = numeroCasualeNormalizzato;

            if (distribuzioneCoordinate(valoreTrial)/distribuzioneCoordinate(valoreControllo) > rapporto)
            {
                valoreControllo = valoreTrial;
                valoriGenerati[i] = valoreTrial * range + offset;
            }
            else
            {
                valoriGenerati[i] = valoreControllo * range + offset;
            }
        }

    }

    return;
}

void metropolisVelocity (long seed, long modulo, long a, int N,
                         double offset, double range, double* valoriGenerati)
{

    double valoreControllo;
    double valoreTrial;

    long numeroCasualeGrezzo = (a * seed) % modulo;
    double numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.) / (modulo * 1.);

    valoreControllo = numeroCasualeNormalizzato;


    for (int i = 0; i < N; i++)
    {
        numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
        numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

        valoreTrial = numeroCasualeNormalizzato;

        if (distribuzioneVelocity(valoreTrial) > distribuzioneVelocity(valoreControllo))
        {
            valoreControllo = valoreTrial;
            valoriGenerati[i] = valoreTrial * range + offset;
        }
        else 
        {
            numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
            numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

            double rapporto = numeroCasualeNormalizzato;

            if (distribuzioneVelocity(valoreTrial)/distribuzioneVelocity(valoreControllo) > rapporto)
            {
                valoreControllo = valoreTrial;
                valoriGenerati[i] = valoreTrial * range + offset;
            }
            else
            {
                valoriGenerati[i] = valoreControllo * range + offset;
            }
        }

    }

    return;
}

int counter = 1;

void metropolisSingolo (long seed, long modulo, long a, 
                        double offset, double range, double &valoreSingolo)
{
    seed += counter;
    counter++;

    double valoreControllo;
    double valoreTrial;

    long numeroCasualeGrezzo = (a * seed) % modulo;
    double numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.) / (modulo * 1.);

    valoreControllo = numeroCasualeNormalizzato;

    numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
    numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);

    valoreTrial = numeroCasualeNormalizzato;
    
    if (distribuzioneCoordinate(valoreTrial) > distribuzioneCoordinate(valoreControllo))
    {
        valoreControllo = valoreTrial;
        valoreSingolo = valoreTrial * range + offset;
    }
    else 
    {
        numeroCasualeGrezzo = (a * numeroCasualeGrezzo) % modulo;
        numeroCasualeNormalizzato = (numeroCasualeGrezzo * 1.)/(modulo * 1.);
    
        double rapporto = numeroCasualeNormalizzato;
    
        if (distribuzioneCoordinate(valoreTrial)/distribuzioneCoordinate(valoreControllo) > rapporto)
        {
            valoreControllo = valoreTrial;
            valoreSingolo = valoreTrial * range + offset;
        }
        else
        {
            valoreSingolo = valoreControllo * range + offset;
        }
    }

    return;
}

void differenceCheck (long seedX, long seedY, long seedZ, long modulo, long a, int N,
                      double offsetX, double offsetY, double offsetZ, 
                      double rangeX, double rangeY, double rangeZ,
                      double* coordinateX, double* coordinateY, double* coordinateZ)
{
    for (int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if (i==j) continue;

            if ((coordinateX[i]==coordinateX[j]) and (coordinateY[i]==coordinateY[j]) and (coordinateZ[i]==coordinateZ[j]))
            {
                metropolisSingolo(seedX, modulo, a, offsetX, rangeX, coordinateX[i]);
                metropolisSingolo(seedY, modulo, a, offsetY, rangeY, coordinateY[i]);
                metropolisSingolo(seedZ, modulo, a, offsetZ, rangeZ, coordinateZ[i]);
                j--;
            }
        }
    }

    return;
}

double pressione(int N, double volume, double massa, 
                 double* velocityX, double* velocityY, double* velocityZ)
{
    double somma = 0;

    for (int i = 0; i < N; i++)
    {
        somma += pow(velocityX[i], 2) + pow(velocityY[i], 2) + pow(velocityZ[i], 2);
    }

    return 1/(3 * volume) * massa * somma;
}

double temperatura (int N, double massa,
                    double* velocityX, double* velocityY, double* velocityZ)
{
    double somma = 0;

    for (int i = 0; i < N; i++)
    {
        somma += pow(velocityX[i], 2) + pow(velocityY[i], 2) + pow(velocityZ[i], 2);
    }

    return 1/((3 * N - 3) * kB) * massa * somma;
    
}

double energiaTotale (int N, double massa,
                      double* coordinateX, double* coordinateY, double* coordinateZ,
                      double* velocityX, double* velocityY, double* velocityZ)
{
    double sommaVelocity = 0;

    for (int i = 0; i < N; i++)
    {
        sommaVelocity += pow(velocityX[i], 2) + pow(velocityY[i], 2) + pow(velocityZ[i], 2);
    }

    double energiaCinetica = 0.5 * massa * sommaVelocity;

    double sommaDistanze = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j==i) continue;

            double distanza = sqrt(pow(coordinateX[i] - coordinateX[j],2) + pow(coordinateY[i] - coordinateY[j],2) + pow(coordinateZ[i] - coordinateZ[j],2) );
            sommaDistanze += 1./distanza;
        }
    }
    
    double energiaPotenzialeGravitazionale = G * massa * massa * sommaDistanze;
    energiaPotenzialeGravitazionale = 0;

    double potenzialeLennardJones = 0;
 
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j==i) continue;

            double distanza = sqrt(pow(coordinateX[i] - coordinateX[j],2) + pow(coordinateY[i] - coordinateY[j],2) + pow(coordinateZ[i] - coordinateZ[j],2) );
        
            potenzialeLennardJones += 4 * epsilonArgon * (pow(sigmaArgon/distanza, 12) - pow(sigmaArgon/distanza, 6));
        }
    }

    if (!interazioneParticelle) potenzialeLennardJones = 0;

    double energia = energiaCinetica + energiaPotenzialeGravitazionale + potenzialeLennardJones;

    return energia;    
}

void print (int N, double volume, double massa, int iterazioneAttuale,
            double* coordinateX, double* coordinateY, double* coordinateZ,
            double* velocityX, double* velocityY, double* velocityZ, 
            ofstream &fileOutput)
{
    double p = pressione(N, volume, massa, velocityX, velocityY, velocityZ);
    double T = temperatura(N, massa, velocityX, velocityY, velocityZ);
    double energia = energiaTotale(N, massa, coordinateX, coordinateY, coordinateZ, velocityX, velocityY, velocityZ);

    fileOutput << iterazioneAttuale << '\t' << energia << '\t' << p*volume << '\t' << N * kB * T << endl;  
    
    ofstream filePosizioni;
    if (interazioneParticelle) filePosizioni.open(".coordinateInterazione/" + to_string(iterazioneAttuale) + ".txt");
    else filePosizioni.open(".coordinateNoInterazione/" + to_string(iterazioneAttuale) +".txt");

    filePosizioni << "x\ty\tz" << endl;
    for (int i = 0; i < N; i++)
    {
        filePosizioni << coordinateX[i] << '\t' << coordinateY[i] << '\t' << coordinateZ[i] << endl;   
    }
    
    filePosizioni.close();

    ofstream fileVelocity;
    if (interazioneParticelle) fileVelocity.open(".velocitàInterazione/" + to_string(iterazioneAttuale) + ".txt");
    else fileVelocity.open(".velocitàNoInterazione/" + to_string(iterazioneAttuale) + ".txt");

    fileVelocity << "Vx\tVy\tVz" << endl;
    for (int i = 0; i < N; i++)
    {
        fileVelocity << velocityX[i] << '\t' << velocityY[i] << '\t' << velocityZ[i] << endl;
    }
    
    fileVelocity.close();

    return;
}

void forzaLennardJones (int N, int i, int j, 
                        double* coordinateX, double* coordinateY, double* coordinateZ,
                        double* forzaLJ)
{
    if(!interazioneParticelle) return;

    double distanza = sqrt(pow(coordinateX[i] - coordinateX[j], 2) + 
                           pow(coordinateY[i] - coordinateY[j], 2) + 
                           pow(coordinateZ[i] - coordinateZ[j], 2));

    forzaLJ[0] = 24. * epsilonArgon/pow(distanza, 2) * (2 * pow(sigmaArgon/distanza, 12) - pow(sigmaArgon/distanza, 6)) * (coordinateX[i] - coordinateX[j]); 
    forzaLJ[1] = 24. * epsilonArgon/pow(distanza, 2) * (2 * pow(sigmaArgon/distanza, 12) - pow(sigmaArgon/distanza, 6)) * (coordinateY[i] - coordinateY[j]); 
    forzaLJ[2] = 24. * epsilonArgon/pow(distanza, 2) * (2 * pow(sigmaArgon/distanza, 12) - pow(sigmaArgon/distanza, 6)) * (coordinateZ[i] - coordinateZ[j]); 

    return;
}

void motoVelocityVerlet(int N, double dt, double massa, double volume,
                        double* coordinateX, double* coordinateY, double* coordinateZ,
                        double* velocityX, double* velocityY, double* velocityZ,
                        double offsetX, double offsetY, double offsetZ,
                        double rangeX, double rangeY, double rangeZ)
{
    //array contenenti posizioni e velocità all'istante n+1
    double* futureX = new double[N]; 
    double* futureY = new double[N]; 
    double* futureZ = new double[N];

    double* futureVX = new double[N]; 
    double* futureVY = new double[N]; 
    double* futureVZ = new double[N];

    double* forzaLJ = new double[3];

    forzaLJ[0] = 0;
    forzaLJ[1] = 0;
    forzaLJ[2] = 0;

    for (int i = 0; i < N; i++) //generazione delle posizioni all'instante di tempo n+1
    {
        double forzaX = 0, forzaY = 0, forzaZ = 0;

        for (int j = 0; j < N; j++) //calcolo forza totale agente sulla particella
        {
            if (j==i) continue;

            forzaLennardJones(N, i, j, coordinateX, coordinateY, coordinateZ, forzaLJ);
            forzaX += forzaLJ[0];
            forzaY += forzaLJ[1];
            forzaZ += forzaLJ[2];

        }

        futureX[i] = coordinateX[i] +  velocityX[i] * dt + 1/(2 * massa) * forzaX * pow(dt, 2); //O(dt^3)
        futureY[i] = coordinateY[i] +  velocityY[i] * dt + 1/(2 * massa) * forzaY * pow(dt, 2);
        futureZ[i] = coordinateZ[i] +  velocityZ[i] * dt + 1/(2 * massa) * forzaZ * pow(dt, 2);   

    }
    
    for (int i = 0; i < N; i++) //generazione delle velocità all'istante di tempo n+1
    {
        double forzaX = 0, forzaY = 0, forzaZ = 0;
        double forzaFutureX = 0, forzaFutureY = 0, forzaFutureZ = 0;

        for (int j = 0; j < N; j++) //calcolo forza totale agente sulla particella agli istanti n e n+1
        {
            if (j==i) continue;

            forzaLennardJones(N, i, j, coordinateX, coordinateY, coordinateZ, forzaLJ);
            forzaX += forzaLJ[0];
            forzaY += forzaLJ[1];
            forzaZ += forzaLJ[2];
            
            forzaLennardJones(N, i, j, futureX, futureY, futureZ, forzaLJ);
            forzaFutureX += forzaLJ[0];
            forzaFutureY += forzaLJ[1];
            forzaFutureZ += forzaLJ[2];
        }

        futureVX[i] = velocityX[i] + 1/(2 * massa) * (forzaX + forzaFutureX) * dt; //O(dt^2)
        futureVY[i] = velocityY[i] + 1/(2 * massa) * (forzaY + forzaFutureY) * dt;
        futureVZ[i] = velocityZ[i] + 1/(2 * massa) * (forzaZ + forzaFutureZ) * dt;   

        //se la particella ha superato i lati del cubo, la velocità è invertita

        if ((futureX[i] <= offsetX) or (futureX[i] >= offsetX + rangeX))
        {
            futureVX[i] = - futureVX[i];
        }
        
        if ((futureY[i] <= offsetY) or (futureY[i] >= offsetY + rangeY))
        {
            futureVY[i] = - futureVY[i];
        }
        
        if ((futureZ[i] <= offsetZ) or (futureZ[i] >= offsetZ + rangeZ))
        {
            futureVZ[i] = - futureVZ[i];
        }

    }

    for (int i = 0; i < N; i++) //porre l'istante n+1 = n
    {
        coordinateX[i] = futureX[i];
        coordinateY[i] = futureY[i];
        coordinateZ[i] = futureZ[i];

        velocityX[i] = futureVX[i];
        velocityY[i] = futureVY[i];
        velocityZ[i] = futureVZ[i];
    }

   delete [] futureX; 
   delete [] futureY; 
   delete [] futureZ; 
   delete [] futureVX; 
   delete [] futureVY; 
   delete [] futureVZ;
   delete [] forzaLJ;

   return;
}

int main()
{
    //si è usato il long perchè spesso il prodotto a * seed andava in overflow e dava un numero negativo
    long modulo = 2147483647, a = 16807;

    double offset = 0., range = 1.e-2; 

    double offsetX = offset, offsetY = offset, offsetZ = offset;
    double rangeX = range, rangeY = range, rangeZ = range;

    //il volume è un cubo di lato range

    double volume = rangeX * rangeY * rangeZ;

    double massa = UMA * numeroMassaArgon;
    double dt = 1e-4;

    double offsetV = -1e0, rangeV = 2e0;

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
    cout << "Inserire il numero di punti da generare: ";
    cin >> N;

    int numeroIterazioni;
    cout << "Inserire il numero di iterazioni: ";
    cin >> numeroIterazioni;

    long seedX, seedY, seedZ;
    cout << "Inserire i seed per le tre coordinate: "; 
    cin >> seedX >> seedY >> seedZ;

    long seedVX, seedVY, seedVZ;
    cout << "Inserire i seed per le velocità: ";
    cin >> seedVX >> seedVY >> seedVZ;


    //generazione casuale delle posizioni e delle velocità iniziali

    double* coordinateX = new double[N];
    double* coordinateY = new double[N];
    double* coordinateZ = new double[N];

    metropolisCoordinate(seedX, modulo, a, N, offsetX, rangeX, coordinateX);
    metropolisCoordinate(seedY, modulo, a, N, offsetY, rangeY, coordinateY);
    metropolisCoordinate(seedZ, modulo, a, N, offsetZ, rangeZ, coordinateZ);

    double* velocityX = new double[N];
    double* velocityY = new double[N];
    double* velocityZ = new double[N];

    metropolisVelocity(seedVX, modulo, a, N, offsetV, rangeV, velocityX);
    metropolisVelocity(seedVY, modulo, a, N, offsetV, rangeV, velocityY);
    metropolisVelocity(seedVZ, modulo, a, N, offsetV, rangeV, velocityZ);

    //controllo che non vi siano particelle nello stesso punto e nel caso ne cambio la posizione
    differenceCheck(seedX, seedY, seedZ, modulo, a, N, offsetX, offsetY, offsetZ, rangeX, rangeY, rangeZ, coordinateX, coordinateY, coordinateZ);

    int iterazioneAttuale = 0;

    ofstream fileOutput;

    if (interazioneParticelle) fileOutput.open("outputInterazione.txt");
    else fileOutput.open("outputNoInterazione.txt");

    fileOutput << "Iterazione\tEnergia\tpV\tNkBT" << endl;

    while (iterazioneAttuale < numeroIterazioni)
    {
        print(N, volume, massa, iterazioneAttuale, coordinateX, coordinateY, coordinateZ, velocityX, velocityY, velocityZ, fileOutput);

        motoVelocityVerlet(N, dt, massa, volume, 
                           coordinateX, coordinateY, coordinateZ, 
                           velocityX, velocityY, velocityZ, 
                           offsetX, offsetY, offsetZ, 
                           rangeX, rangeY, rangeZ);
        
        iterazioneAttuale++;      
    }

    fileOutput.close();

    delete [] coordinateX;
    delete [] coordinateY;
    delete [] coordinateZ;
    delete [] velocityX;
    delete [] velocityY;
    delete [] velocityZ;

    return 0;
}