# sesion4_traduccion_patrones.py


from Bio import SeqIO  
from Bio.Seq import Seq  


def encontrar_motivo(secuencia, motivo):
    """Devuelve todas las posiciones (base 1) donde aparece el motivo."""
    secuencia = str(secuencia).upper()
    posiciones = []
    inicio = secuencia.find(motivo)
    while inicio != -1:
        posiciones.append(inicio + 1)
        inicio = secuencia.find(motivo, inicio + 1)
    return posiciones


# Traducción de una secuencia de ADN a proteína

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
proteina = dna.translate()

print("PARTE 1: TRADUCCIÓN DE ADN A PROTEÍNA")
print("ADN      :", dna)
print("Longitud :", len(dna), "nucleótidos")
print("Proteína :", proteina)
print("Hasta el primer codón de parada:", dna.translate(to_stop=True))

# Búsqueda de motivos genéticos

motivos = ["ATG", "GCGG", "TGA"]
resultados = []

print("-" * 60)
print("PARTE 2: BÚSQUEDA DE MOTIVOS EN LA SECUENCIA")
print(f"{'Motivo':<10}{'Veces':<10}{'Posiciones (base 1)'}")

for motivo in motivos:
    posiciones = encontrar_motivo(dna, motivo)
    texto_pos = ", ".join(map(str, posiciones)) if posiciones else "No encontrado"
    resultados.append((motivo, len(posiciones), texto_pos))
    print(f"{motivo:<10}{len(posiciones):<10}{texto_pos}")

# Aplicación a una secuencia real del archivo FASTA

print("-" * 60)
print("PARTE 3: ANÁLISIS DE UNA SECUENCIA DE ls_orchid.fasta")


registro = next(SeqIO.parse("ls_orchid.fasta", "fasta"))
secuencia_real = registro.seq
recorte = len(secuencia_real) - (len(secuencia_real) % 3)  # múltiplo de 3
proteina_real = secuencia_real[:recorte].translate()

print("ID          :", registro.id)
print("Longitud    :", len(secuencia_real), "nucleótidos")
print("Proteína (primeros 50 aa):", str(proteina_real)[:50])

resultados_real = []
print(f"\n{'Motivo':<10}{'Veces':<10}{'Primeras posiciones (base 1)'}")
print("-" * 60)
for motivo in motivos:
    posiciones = encontrar_motivo(secuencia_real, motivo)
    texto_pos = ", ".join(map(str, posiciones[:10])) if posiciones else "No encontrado"
    resultados_real.append((motivo, len(posiciones), texto_pos))
    print(f"{motivo:<10}{len(posiciones):<10}{texto_pos}")

# Guardar la tabla de resultados en un archivo .txt

with open("resultados_sesion4.txt", "w", encoding="utf-8") as archivo:
    archivo.write("RESULTADOS - SESIÓN 4: TRADUCCIÓN Y BÚSQUEDA DE PATRONES\n")
    archivo.write("=" * 60 + "\n\n")
    archivo.write(f"ADN de ejemplo : {dna}\n")
    archivo.write(f"Proteína       : {proteina}\n\n")
    archivo.write(f"{'Motivo':<10}{'Veces':<10}{'Posiciones (base 1)'}\n")
    archivo.write("-" * 60 + "\n")
    for motivo, veces, pos in resultados:
        archivo.write(f"{motivo:<10}{veces:<10}{pos}\n")

    archivo.write(f"\nSecuencia real: {registro.id} ({len(secuencia_real)} nt)\n")
    archivo.write(f"{'Motivo':<10}{'Veces':<10}{'Primeras posiciones (base 1)'}\n")
    archivo.write("-" * 60 + "\n")
    for motivo, veces, pos in resultados_real:
        archivo.write(f"{motivo:<10}{veces:<10}{pos}\n")

print("\nResultados guardados en 'resultados_sesion4.txt'")