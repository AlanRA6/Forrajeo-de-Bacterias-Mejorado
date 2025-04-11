from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
import numpy
import copy
import pandas as pd
import os
from fastaReader import fastaReader

def main():
    numeroDeBacterias = 4
    numRandomBacteria = 1
    iteraciones = 3
    tumbo = 5  # Número de gaps a insertar
    nado = 3
    secuencias = list()

    secuencias = fastaReader().seqs
    names = fastaReader().names

    # Convierte las secuencias en listas de caracteres
    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])

    globalNFE = 0  # Número de evaluaciones de la función objetivo

    dAttr = 0.5
    wAttr = 0.01
    hRep = 0.3
    wRep = 0.005

    manager = Manager()
    numSec = len(secuencias)
    print("numSec: ", numSec)

    poblacion = manager.list(range(numeroDeBacterias))
    names = manager.list(names)
    NFE = manager.list(range(numeroDeBacterias))

    # Lista para almacenar los resultados de cada iteración
    resultados = []

    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)

    operadorBacterial = bacteria(numeroDeBacterias)
    veryBest = [None, None, None]  # índice, fitness, secuencias

    # Registra el tiempo de inicio
    start_time = time.time()

    print("Población inicial ...")
    poblacionInicial()

    for it in range(iteraciones):
        print(f"Iteración {it + 1} de {iteraciones} - Tumbo ...")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)
        print("Tumbo realizado - Cuadrando ...")
        operadorBacterial.cuadra(numSec, poblacion)
        print("Población cuadrada - Creando gran lista de pares...")
        operadorBacterial.creaGranListaPares(poblacion)
        print("Gran lista creada - Evaluando BLOSUM en paralelo...")
        operadorBacterial.evaluaBlosum()  # Paralelo
        print("BLOSUM evaluado - Creando tablas de atracción y repulsión...")
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        operadorBacterial.creaTablaInteraction()
        print("Tabla de interacción creada - Creando tabla de fitness...")
        operadorBacterial.creaTablaFitness()
        print("Tabla de fitness creada")
        globalNFE += operadorBacterial.getNFE()
        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
        if (veryBest[0] is None) or (bestFitness > veryBest[1]):
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = copy.deepcopy(poblacion[bestIdx])
        operadorBacterial.replaceWorst(poblacion, veryBest[0])

        # Guardar datos de la iteración
        resultados.append({
            "Iteración": it + 1,
            "Fitness": veryBest[1],
            "Tiempo": time.time() - start_time,
            "Interacción": operadorBacterial.tablaInteraction[veryBest[0]],
            "BLOSUM": operadorBacterial.blosumScore[veryBest[0]]
        })

        operadorBacterial.resetListas(numeroDeBacterias)

    print("Very Best: ", veryBest)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    # Crear la carpeta si no existe
    if not os.path.exists("resultados_corridas"):
        os.makedirs("resultados_corridas")
        
    # Guardar resultados indiviudales en un archivo CSV dentro de la carpeta
    df = pd.DataFrame(resultados)
    csv_path = os.path.join("resultados_corridas", f"resultados_corrida_{time.strftime('%Y%m%d_%H%M%S')}.csv")
    df.to_csv(csv_path, index=False)

    return resultados

if __name__ == "__main__":
    main()