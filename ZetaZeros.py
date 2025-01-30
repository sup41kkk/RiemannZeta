import numpy as np

def eta(s, max_terms=100000):
    """Вычисление η-функции Дирихле с использованием векторизации."""
    n = np.arange(1, max_terms + 1)
    terms = (-1) ** (n - 1) / np.power(n, s)
    return np.sum(terms)

def zeta(s, max_terms=100000):
    """Вычисление ζ-функции через η-функцию."""
    return eta(s, max_terms) / (1 - 2 ** (1 - s))

def find_zero_in_interval(a, b, tol=1e-5, max_iter=100):
    """Поиск нуля в интервале [a, b] методом золотого сечения."""
    gr = (np.sqrt(5) + 1) / 2
    c = b - (b - a)/gr
    d = a + (b - a)/gr
    
    for _ in range(max_iter):
        if abs(c - d) < tol:
            break
            
        fc = np.abs(zeta(0.5 + 1j*c))
        fd = np.abs(zeta(0.5 + 1j*d))
        
        if fc < fd:
            b = d
        else:
            a = c
        
        c = b - (b - a)/gr
        d = a + (b - a)/gr
    
    return (a + b)/2

def find_zeros(a, b, step=0.5, threshold=1e-3):
    """
    Поиск всех нулей в диапазоне [a, b]
    Параметры:
        step - шаг первичного сканирования
        threshold - порог значения |ζ| для идентификации нуля
    """
    zeros = []
    t_values = np.arange(a, b, step)
    moduli = [np.abs(zeta(0.5 + 1j*t)) for t in t_values]
    
    # Ищем локальные минимумы модуля
    for i in range(1, len(moduli)-1):
        if moduli[i] < moduli[i-1] and moduli[i] < moduli[i+1]:
            # Уточняем ноль в окрестности минимума
            t_zero = find_zero_in_interval(t_values[i-1], t_values[i+1])
            z_value = zeta(0.5 + 1j*t_zero)
            
            if np.abs(z_value) < threshold:
                zeros.append(t_zero)
    
    # Удаление дубликатов и сортировка
    zeros = np.unique(np.round(zeros, decimals=3))
    return zeros

zeros = find_zeros(14, 100, step=0.2)
print("Найденные нули на критической линии:")
for i, t in enumerate(zeros, 1):
    print(f"Номер {i}: t ≈ {t:.4f}")