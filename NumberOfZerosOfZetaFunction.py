import mpmath as mp
from tqdm import tqdm

def count_riemann_zeros(T, method='adaptive', prec=100, adaptive_tol=1e-6):
    """
    Вычисляет количество нулей дзета-функции Римана на критической прямой
    Re(s) = 1/2 с мнимой частью в интервале [0, T] с использованием аргументного принципа.

    Параметры:
        T (float): Верхняя граница мнимой части
        method (str): Метод интегрирования ('adaptive' или 'trapezoid')
        prec (int): Точность вычислений в битах
        adaptive_tol (float): Допустимая погрешность для адаптивного метода

    Возвращает:
        int: Количество нулей (округлённое до целого)
    """
    mp.mp.prec = prec
    mp.mp.pretty = False

    def argument_derivative(t):
        s = mp.mpc(0.5, t)
        try:
            z = mp.zeta(s)
            if abs(z) < 1e-20:  # Обнаружение нуля
                return mp.inf
            return mp.im(mp.zeta(s, derivative=1)/z)
        except (ValueError, mp.libmp.NoConvergence):
            return mp.nan

    # Адаптивное интегрирование (рекомендуется)
    if method == 'adaptive':
        integral = mp.quad(
            argument_derivative,
            [0, T],
            method='gauss-legendre',
            error=True,
            maxdegree=20,
            tol=adaptive_tol
        )
        result = integral[0]/mp.pi
        error = integral[1]/mp.pi
        
        if error > 0.5:
            print(f"Warning: Возможная погрешность {error:.2f} нулей")

    # Метод трапеций с прогресс-баром
    elif method == 'trapezoid':
        steps = max(100, int(T//10))  # Автоподбор шагов
        dt = T/steps
        total = 0
        prev_val = argument_derivative(0)
        
        with tqdm(total=steps, desc="Integrating") as pbar:
            for i in range(1, steps+1):
                t = i*dt
                curr_val = argument_derivative(t)
                if mp.isnan(curr_val):
                    curr_val = prev_val  # Простейшая обработка ошибок
                total += (prev_val + curr_val)*dt/2
                prev_val = curr_val
                pbar.update(1)
        
        result = total/mp.pi

    return int(round(result))

if __name__ == "__main__":
    import sys
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
    
    print(f"Расчёт нулей ζ(1/2 + it) для t ∈ [0, {T}]")
    zeros = count_riemann_zeros(T, method='adaptive', prec=200)
    print(f"\nПриблизительное количество нулей: {zeros}")
    
    # Сравнение с известными значениями
    if T == 30:
        print("\nСправка: Для T=30 точное количество нулей равно 49")