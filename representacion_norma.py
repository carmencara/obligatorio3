import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del archivo
datos = np.loadtxt("norma.txt")

# Extraer la columna "j" y la columna "abs(phi_j_n)"
j = datos[:, 0]
phi_j_n = datos[:, 1]

# Crear el gráfico
fig, ax = plt.subplots()
ax.plot(j, phi_j_n)

# Establecer las etiquetas de los ejes
ax.set_xlabel("j")
ax.set_ylabel(r"|$\phi_{j,n}$|")

# Mostrar y guardad el gráfico
plt.savefig('norma.png', dpi=300, bbox_inches='tight')
plt.show()
