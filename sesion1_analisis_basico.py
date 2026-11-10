# 1. Leer y extraer la segunda secuencia del archivo FASTA
numero_secuencia = 0
encabezado = ""
secuencia = ""

with open("ls_orchid.fasta", "r", encoding="utf-8") as archivo:
    for linea in archivo:
        linea = linea.strip()

        # Identificar el inicio de cada secuencia
        if linea.startswith(">"):
            numero_secuencia += 1

            # Guardar el encabezado de la segunda secuencia
            if numero_secuencia == 2:
                encabezado = linea

            # Detener la lectura cuando comienza la tercera
            elif numero_secuencia > 2:
                break

        # Unir las líneas de ADN de la segunda secuencia
        elif numero_secuencia == 2:
            secuencia += linea.upper()


# 2. Mostrar la secuencia extraída
if numero_secuencia >= 2 and secuencia:
    print("SEGUNDA SECUENCIA EXTRAÍDA")
    print("Encabezado:", encabezado)
    print("Secuencia de ADN:")
    print(secuencia)
    print("Longitud:", len(secuencia))
else:
    print("No se encontró la segunda secuencia.")
    
    

def contar_bases(secuencia):
    bases = {
        "A": secuencia.count("A"),
        "T": secuencia.count("T"),
        "C": secuencia.count("C"),
        "G": secuencia.count("G")
    }

    return bases

def calcular_porcentaje_gc(secuencia):
    total = len(secuencia)
    cantidad_gc = secuencia.count("G") + secuencia.count("C")

    if total == 0:
        return 0

    porcentaje = (cantidad_gc / total) * 100

    return porcentaje

# Calcular resultados
longitud = len(secuencia)
conteo = contar_bases(secuencia)
porcentaje_gc = calcular_porcentaje_gc(secuencia)

# Mostrar resultados
print("\nRESULTADOS DEL ANÁLISIS")
print("----------------------------")
print("Longitud:", longitud, "bases")

print("\nConteo de bases:")
for base, cantidad in conteo.items():
    print(f"{base}: {cantidad}")

print(f"\nPorcentaje GC: {porcentaje_gc:.2f}%")