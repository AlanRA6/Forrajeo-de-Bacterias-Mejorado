import math
from multiprocessing import Manager, Pool, managers
from pickle import FALSE, TRUE
from evaluadorBlosum import evaluadorBlosum
import numpy
from fastaReader import fastaReader
import random

#from copy import copy
import copy
import concurrent.futures
import numpy as np

class bacteria():
    
    def __init__(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
        self.evaluador = evaluadorBlosum()  # Inicializar evaluador BLOSUM
        self.max_fitness = 1e6  # Límite superior para valores de fitness
            
    def resetListas(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
  
    def cuadra(self, numSec, poblacion):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])
            bacterTmp = bacterTmp[:numSec]
            maxLen = max(len(seq) for seq in bacterTmp)
            for t in range(numSec):
                gap_count = maxLen - len(bacterTmp[t])
                if gap_count > 0:
                    bacterTmp[t].extend(["-"] * gap_count)
            poblacion[i] = tuple(bacterTmp)
        return poblacion

    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1

    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]

    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True
    
    def tumbo(self, numSec, poblacion, numGaps):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])
            
            for _ in range(numGaps):
                seqnum = random.randint(0, len(bacterTmp)-1)
                pos = random.randint(0, len(bacterTmp[seqnum]))
                bacterTmp[seqnum].insert(pos, "-")
            
            # Asegurar igual longitud
            max_len = max(len(seq) for seq in bacterTmp)
            for seq in bacterTmp:
                while len(seq) < max_len:
                    seq.append('-')
            
            poblacion[i] = tuple(bacterTmp)
        return poblacion
       
    def creaGranListaPares(self, poblacion):   
        for i in range(len(poblacion)):
            pares = list()
            bacterTmp = list(poblacion[i])
            
            for j in range(len(bacterTmp[0])):  # Usar longitud de la primera secuencia
                column = self.getColumn(bacterTmp, j)
                pares += self.obtener_pares_unicos(column)
            
            self.granListaPares[i] = pares
        return self.granListaPares

    def evaluaFila(self, args):
        fila, num = args
        score = 0
        evaluador = evaluadorBlosum()
        for par in fila:
            score += evaluador.getScore(par[0], par[1])
        # Normalizar dividiendo por el número de pares
        normalized_score = score / (len(fila) + 1e-8)  # +1e-8 evita división por cero
        return (num, normalized_score)
        

    def evaluaBlosum(self):
        with Pool() as pool:
            # Preparar argumentos serializables
            args = [(copy.deepcopy(self.granListaPares[i]), i) 
                for i in range(len(self.granListaPares))]
            
            # Usar map normal en lugar de starmap
            results = pool.map(self.evaluaFila, args)
            
            # Actualizar los scores
            for num, score in results:
                self.blosumScore[num] = score

    def getColumn(self, bacterTmp, colNum):
        return [seq[colNum] for seq in bacterTmp]
            
    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                # Saltar pares con valores None o inválidos
                if columna[i] is not None and columna[j] is not None:
                    par = tuple(sorted([columna[i], columna[j]]))
                    pares_unicos.add(par)
        return list(pares_unicos)

    def compute_diff(self, args):
        indexBacteria, otherBlosumScore, self.blosumScore, d, w = args
        
        # 1. Calcular diferencia con protección contra divisiones/overflow
        diff = (self.blosumScore[indexBacteria] - otherBlosumScore)
        
        # 2. Normalizar la diferencia usando el rango de scores
        max_possible_diff = 1000  # Valor máximo esperado entre scores
        normalized_diff = diff / max_possible_diff
        
        # 3. Aplicar función sigmoide para mantener valores acotados
        safe_exp_arg = min(w * (normalized_diff ** 2), 100)  # Limitar argumento exponencial
        exp_result = numpy.exp(safe_exp_arg)
        
        # 4. Calcular resultado final con protección numérica
        result = d * exp_result
        
        self.NFE[indexBacteria] += 1
        
        # 5. Verificar y corregir posibles overflows residuales
        if numpy.isinf(result) or result > 1e100:
            return 1e100  # Valor grande pero manejable
        return result
    
    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        with Pool() as pool:
            args = [(indexBacteria, otherBlosumScore, self.blosumScore, d, w) 
                   for otherBlosumScore in self.blosumScore]
            results = pool.map(self.compute_diff, args)
            pool.close()
            pool.join()
    
        total = sum(results)
        if atracTrue:
            self.tablaAtract[indexBacteria] = total
        else:
            self.tablaRepel[indexBacteria] = total
  
    def creaTablaAtract(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, TRUE)

    def creaTablaRepel(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, FALSE)
    
    def creaTablasAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(self.creaTablaAtract, poblacion, dAttr, wAttr)
            executor.submit(self.creaTablaRepel, poblacion, dRepel, wRepel)
            
    def creaTablaInteraction(self):
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def creaTablaFitness(self):
        for i in range(len(self.tablaInteraction)):
            # 1. Obtener valores base con protección
            valorBlsm = self.blosumScore[i] if not numpy.isinf(self.blosumScore[i]) else 0
            valorInteract = self.tablaInteraction[i] if not numpy.isinf(self.tablaInteraction[i]) else 0
            
            # 2. Escalar valores para mantener equilibrio
            scale_factor = 0.1  # Ajustar según necesidad
            scaled_blosum = self.blosumScore[i] / 100.0
            scaled_interaction = self.tablaInteraction[i] * 10.0
            
            # 3. Calcular fitness con protección
            try:
                valorFitness = scaled_blosum + scaled_interaction
            except:
                valorFitness = scaled_blosum  # Fallback seguro
                
            # 4. Verificar límites numéricos
            if numpy.isinf(valorFitness) or numpy.isnan(valorFitness):
                valorFitness = 0  # Valor por defecto seguro
                
            # 5. Asignar valor final
            self.tablaFitness[i] = valorFitness
            
            # 6. Opcional: Loggear valores extremos para diagnóstico
            if abs(valorFitness) > 1e6:
                print(f"Valor extremo en bacteria {i}: BLOSUM={valorBlsm}, Interacción={valorInteract}")
    
    def getNFE(self):
        return sum(self.NFE)
        
    def obtieneBest(self, globalNFE):
        bestIdx = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] > self.tablaFitness[bestIdx]:
                bestIdx = i
        print(f"Best: {bestIdx}, Fitness: {self.tablaFitness[bestIdx]}, "
              f"BLOSUM: {self.blosumScore[bestIdx]}, "
              f"Interaction: {self.tablaInteraction[bestIdx]}, "
              f"NFE: {globalNFE}")
        return bestIdx, self.tablaFitness[bestIdx]

    def replaceWorst(self, poblacion, best):
        worst = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] < self.tablaFitness[worst]:
                worst = i
        poblacion[worst] = copy.deepcopy(poblacion[best])
        
    
    def busquedaLocal(self, poblacion, numSec, max_iter=5):
            poblacion = list(poblacion)
            
            for i in range(len(poblacion)):
                if self.blosumScore[i] < -5000:
                    continue
                    
                bacterTmp = list(poblacion[i])
                mejorado = False
                
                max_len = max(len(seq) for seq in bacterTmp)
                for seq in bacterTmp:
                    seq.extend(['-'] * (max_len - len(seq)))
                
                for _ in range(max_iter):
                    if not bacterTmp or not bacterTmp[0]:
                        break
                        
                    puntajes_columnas = []
                    for col in range(len(bacterTmp[0])):
                        try:
                            columna = [bacterTmp[seq][col] for seq in range(numSec)]
                            pares = self.obtener_pares_unicos(columna)
                            puntaje = sum(self.evaluador.getScore(p[0], p[1]) for p in pares)
                            puntajes_columnas.append(puntaje)
                        except:
                            continue
                    
                    if not puntajes_columnas:
                        break
                        
                    try:
                         # 1. Identificar peor columna
                        peor_col = numpy.nanargmin(puntajes_columnas)
                    except:
                        break
                        
                    for seq_idx in range(numSec):
                        if bacterTmp[seq_idx][peor_col] == "-":
                            vecino = copy.deepcopy(bacterTmp)
                            del vecino[seq_idx][peor_col]
                            pos = random.randint(0, len(vecino[seq_idx]))
                            vecino[seq_idx].insert(pos, "-")
                            
                            try:
                                pares_vecino = self.creaGranListaPares([vecino])[0]
                                nuevo_blosum = sum(self.evaluador.getScore(p[0], p[1]) for p in pares_vecino)
                                
                                if nuevo_blosum > self.blosumScore[i] + 50:
                                    bacterTmp = vecino
                                    self.blosumScore[i] = nuevo_blosum
                                    mejorado = True
                                    break
                            except:
                                continue
                    
                    if not mejorado:
                        break
                
                poblacion[i] = tuple(bacterTmp)
            
            return poblacion