import numpy as np
import matplotlib.pyplot as plt

def solve_poisson(V, rho, dx, dy, epsilon_0, max_iter=1000, tol=1e-6):
    N, M = V.shape
    V_new = np.zeros_like(V)
    iter_count = 0
    diff = tol + 1
    
    while iter_count < max_iter and diff > tol:
        V_old = V.copy()
        for i in range(1, N-1):
            for j in range(1, M-1):
                V_vc = 0.25 * (V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1])
                V_vsc = 0.25 * (V[i+1,j+1] + V[i-1,j+1] + V[i-1,j-1] + V[i+1,j-1])
                V_new[i,j] = (2/3) * V_vc + (1/3) * V_vsc + (dx**2 / (4 * epsilon_0)) * rho[i,j]
                
        diff = np.max(np.abs(V_new - V_old))
        V, V_new = V_new, V
        iter_count += 1
        
    return V
# Ejemplo de uso
N = 50  # Número de puntos en x
M = 50  # Número de puntos en y
a = 1.0  # Longitud del lado del cuadrado
v0 = 1.0  # Potencial en los bordes
rho = np.zeros((N, M))  # Densidad de carga
rho[N//4:3*N//4, M//4:3*M//4] = 1.0  # Agregando una distribución de carga en el centro

# Paso de la discretización
dx = a / (N - 1)
dy = a / (M - 1)

# Inicialización del potencial
V = np.zeros((N, M))
V[:,0] = -v0  # Condiciones de frontera
V[:,-1] = v0   # Condiciones de frontera

# Resolviendo la ecuación de Poisson
epsilon_0 = 1.0  # Permitividad del vacío
V = solve_poisson(V, rho, dx, dy, epsilon_0)

# Visualización del potencial
plt.figure(figsize=(12, 5))
plt.imshow(V, cmap='viridis', origin='lower', extent=[0, a, 0, a])
plt.colorbar(label='Potencial (V)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potencial eléctrico: Solucion Analítica')
plt.show()

def V0(x, y, a, b, v0):
    tolerance = 1e-10  # Tolerancia para la comparación de punto flotante
    if abs(x - 0) < tolerance or abs(x - a) < tolerance:
        return 0
    elif abs(y - 0) < tolerance:
        return -v0
    elif abs(y - b) < tolerance:
        return v0
    else:
        return 0  # Valor por defecto si no se cumplen ninguna de las condiciones anteriores



def solve_laplace(V, dx, dy, V0, a, b, v0, tol=1e-4, max_iter=10000):
    iter_count = 0
    diff = tol + 1
    while iter_count < max_iter and diff > tol:
        V_old = V.copy()
        for i in range(1, N-1):
            for j in range(1, M-1):
                V[i,j] = 0.25 * (V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1]) + V0(i*dx, j*dy, a, b, v0)
        diff = np.max(np.abs(V - V_old))
        iter_count += 1
    return V

# Dimensiones de la región rectangular
a = 1.0
b = 1.0
N = 50 # Número de puntos en x
M = 50 # Número de puntos en y

# Paso de la discretización
dx = a / (N - 1)
dy = b / (M - 1)

# Inicialización de V
V = np.zeros((N, M))

# Condiciones de frontera
for i in range(N):
    for j in range(M):
        V[i,j] = V0(i*dx, j*dy, a, b, 1.0)

# Resolviendo la ecuación de Laplace
V = solve_laplace(V, dx, dy, V0, a, b, 1.0)

# Definiendo la malla
x = np.linspace(0, a, N)
y = np.linspace(0, b, M)
X, Y = np.meshgrid(x, y)


# Gráfico
plt.figure(figsize=(12, 5))
plt.title("Potencial Eléctrico: Solución numérica")
plt.contourf(X, Y, V, cmap='viridis')
plt.colorbar(label='Potencial (V)')
plt.xlabel('x')
plt.ylabel('y')


