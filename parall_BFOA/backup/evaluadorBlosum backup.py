import blosum as bl #Permite usar la matriz BLOSUM (Asigna un puntaje a cada par de aminoacidos segun su simulitud)

#La clase que manejará la matriz BLOSUM y calculará puntutuaciones
class evaluadorBlosum():
    
    #Constructor, se ejecuta cuando se crea un objeto de evaluadorBlosum
    def __init__(self):

        matrix = bl.BLOSUM(62) #Para calcular puntuaciones entre pares de minoacidos 
        
        self.matrix = matrix #Carga la matriz BLOSUM62 en self.matrix
    
    #Muestra la matriz completa en la pantalla    
    def showMatrix(self):
        print(self.matrix)
    
    #Recibe dos aminoacidos A y B, devuelve la puntuacion de sustitucioin entre ellos segun BLOSUM62
    def getScore(self, A, B):
        #si alguno de los dos es un gap
        #En bioinformática, los gaps representan inserciones o eliminaciones en secuencias de proteínas
        if A == "-" or B == "-":
            return -8
        score = self.matrix[A][B] #Busca en self.matrix la puntuación entre A y B
        return score #lo devuelve como resultado
    
    
    pass




