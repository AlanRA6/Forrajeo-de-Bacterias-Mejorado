import blosum as bl #Permite usar la matriz BLOSUM (Asigna un puntaje a cada par de aminoacidos segun su simulitud)

#La clase que manejará la matriz BLOSUM y calculará puntutuaciones
class evaluadorBlosum():
    
    #Constructor, se ejecuta cuando se crea un objeto de evaluadorBlosum
    def __init__(self):
        # Convertir la matriz BLOSUM a un diccionario serializable
        matrix = bl.BLOSUM(62)
        self.scores = {}
        for A in matrix.keys():
            for B in matrix.keys():
                self.scores[(A, B)] = matrix[A][B]
    
    def getScore(self, A, B):
        if A == "-" or B == "-":
            return -8
        return self.scores.get((A, B), -4)  # Default penalty for unknown pairs
    
    #Muestra la matriz completa en la pantalla    
    def showMatrix(self):
        print(self.matrix)
    
    pass




