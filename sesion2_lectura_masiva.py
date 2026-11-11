
from Bio import SeqIO
from pathlib import Path

# Ruta del archivo FASTA
ruta_archivo = Path(__file__).parent / "ls_orchid.fasta"

# Listas para almacenar los datos
longitudes = []
secuencias = []

# Leer automáticamente todas las secuencias
for record in SeqIO.parse(ruta_archivo, "fasta"):
    secuencias.append(record)
    longitudes.append(len(record.seq))

# Calcular estadísticas generales
total_secuencias = len(secuencias)

if total_secuencias > 0:

    longitud_promedio = sum(longitudes) / total_secuencias

    secuencia_mas_larga = max(
        secuencias,
        key=lambda record: len(record.seq)
    )

    # Mostrar resultados
    print(" ANALISIS AUTOMATIZADO DE SECUENCIAS")
    print(f"Total de secuencias: {total_secuencias}")

    print(f"Longitud promedio: {longitud_promedio:.2f} bases")

    print("\nSECUENCIA MAS LARGA")
    print(f"Identificador: {secuencia_mas_larga.id}")
    print(f"Longitud: {len(secuencia_mas_larga.seq)} bases")

else:
    print("No se encontraron secuencias en el archivo FASTA.")
