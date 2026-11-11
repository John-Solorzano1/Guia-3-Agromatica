
from Bio import SeqIO
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

# Ruta del archivo FASTA
carpeta = Path(__file__).parent
ruta_archivo = carpeta / "ls_orchid.fasta"

# Carpeta donde se guardaran los graficos
carpeta_graficos = carpeta / "graficos"
carpeta_graficos.mkdir(exist_ok=True)

# Leer todas las secuencias
registros = list(SeqIO.parse(ruta_archivo, "fasta"))

if not registros:
    print("No se encontraron secuencias en el archivo FASTA.")
    raise SystemExit

# GRAFICO 1: LONGITUD DE LAS PRIMERAS 10

ids = [record.id for record in registros[:10]]
longitudes = [len(record.seq) for record in registros[:10]]

plt.figure(figsize=(12, 6))

plt.bar(ids, longitudes)

plt.xlabel("Identificador de la secuencia")
plt.ylabel("Longitud (bases)")
plt.title("Longitud de las primeras 10 secuencias")

plt.xticks(rotation=90)
plt.tight_layout()

plt.savefig(
    carpeta_graficos / "grafico_barras_longitudes.png",
    dpi=300
)

plt.show()

# GRAFICO 2: PROPORCION DE BASES A, T, G, C
# Utilizar la segunda secuencia para mantener
# la continuidad con el analisis de la sesion 1
registro = registros[1] if len(registros) > 1 else registros[0]

secuencia = str(registro.seq).upper()

conteo = Counter(secuencia)

bases = ["A", "T", "G", "C"]
cantidades = [conteo.get(base, 0) for base in bases]

plt.figure(figsize=(8, 8))

plt.pie(
    cantidades,
    labels=bases,
    autopct="%1.1f%%",
    startangle=90
)

plt.title("Proporcion de bases nitrogenadas")

plt.tight_layout()

plt.savefig(
    carpeta_graficos / "grafico_pastel_bases.png",
    dpi=300
)

plt.show()

print(" VISUALIZACIONES GENERADAS")
print("Grafico de barras guardado correctamente.")
print("Grafico de pastel guardado correctamente.")
print(f"Secuencia utilizada para el pastel: {registro.id}")
print(f"Longitud analizada: {len(secuencia)} bases")
print(f"Resultados guardados en: {carpeta_graficos}")
