import pandas as pd
from parallel_BFOA import main  # Importa la función principal del algoritmo

if __name__ == "__main__":
    resultados_totales = []

    for corrida in range(30):
        print(f"Ejecutando corrida {corrida + 1}...")
        resultados_corrida = main()  # Ejecuta el algoritmo y obtén los resultados
        resultados_totales.extend(resultados_corrida)

    # Guardar todos los resultados en un archivo CSV
    df = pd.DataFrame(resultados_totales)
    df.to_csv("resultados_30_corridas.csv", index=False)