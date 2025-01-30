import numpy as np
from scipy.linalg import solve

def compute_zeta(s):
    sigma, t = s.real, s.imag
    
    beta = 0.45
    K = int(np.ceil(beta * t))
    N = 15
    
    # Вычисление первой суммы
    n = np.arange(1, K)
    sum1 = np.sum((-1)**n / np.power(n, s))
    
    # Создание коэффициентов
    idx = np.arange(0, 2*N + 1)
    coeff = 1 / np.power(K + idx, s)
    
    # Построение матрицы Q и вектора U
    Q = np.zeros((N, N), dtype=np.complex128)
    U = np.zeros(N, dtype=np.complex128)
    
    for k in range(N):
        U[k] = -coeff[N + k + 1]
        for m in range(N):
            Q[k, m] = coeff[N + k - m]
    
    # Решение системы уравнений
    q = solve(Q, U)
    
    # Вычисление коэффициентов A
    A = np.zeros(N + 1, dtype=np.complex128)
    A[0] = coeff[0]
    
    for n in range(1, N + 1):
        k_slice = np.arange(1, min(n, N) + 1)
        valid_idx = n - k_slice
        sum_val = np.sum(coeff[valid_idx] * q[k_slice - 1])
        A[n] = coeff[n] + sum_val
    
    # Вычисление сумм
    i = np.arange(1, N + 2)
    sumA = np.sum((-1)**(i + 1) * A[:N + 1])
    
    j = np.arange(1, N + 1)
    sumB = 1 + np.sum((-1)**j * q[:N])
    
    r = sumA / sumB
    LI = sum1 + (-1)**K * r
    zeta = -LI / (1 - 2**(1 - s))
    
    return np.round(zeta, 10)

# Запрос ввода комплексного числа
s_input = input("Введите комплексное число s в формате a+bj (например, 0.5+156j): ").replace(' ', '')
try:
    s = complex(s_input)
except ValueError:
    print("Ошибка: некорректный формат числа. Пожалуйста, используйте формат a+bj.")
    exit()

# Вызов функции и вывод результата
result = compute_zeta(s)
print(f"Result: {result:.10f}")