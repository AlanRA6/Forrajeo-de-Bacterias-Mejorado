from evaluadorBlosum import evaluadorBlosum


#Crear unsa instancia de la clase 
evaluator = evaluadorBlosum()
#MOstrar la matriz BLOSUM62
print("Matriz BLOSUM62:")
evaluator.showMatrix()
#Calcular puntutaciones
print("\nPuntutaciones de ejemplos:")
print(f"Score entre A y G: {evaluator.getScore('A', 'G')}") # A (Alanina) y G (Glicina)
print(f"Score entre A y gap: {evaluator.getScore('A', '-')}") # Evalúa un gap (-), devolviendo -8
print(f"Score entre W y L: {evaluator.getScore('W', 'L')}") #Compara Triptófano (W) con Leucina (L).


