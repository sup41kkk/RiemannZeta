# Advanced Riemann Zeta Function Computational Toolkit

This repository provides three Python scripts for advanced computation and analysis of the Riemann zeta function, ζ(s), and its critical properties, including the identification and counting of its zeros along the critical line. These programs employ sophisticated numerical methods, such as Padé approximation, argument principle-based zero counting, and Dirichlet η-function transformations. Below is a detailed breakdown of their functionality and usage.

## Program 1: **High-Precision Computation of ζ(s) Using Padé Approximation**

### Overview
This program utilizes Padé approximations combined with series expansions to evaluate the Riemann zeta function for a given complex input `s`. The method incorporates the solution of a linear system to accurately calculate series coefficients, making it effective for inputs with high imaginary components.

### Features
- High precision computation for inputs with large imaginary parts.
- Accepts complex numbers in the standard format `a+bj`.
- Relies on numpy and scipy for robust matrix algebra and numerical operations.

### Instructions
1. Execute the script.
2. Provide the desired complex number `s` in the format `a+bj` (e.g., `0.5+156j`).
3. The computed value of ζ(s) is returned with a precision of 10 decimal places.

Example:
```
Введите комплексное число s в формате a+bj (например, 0.5+156j): 0.5+30j
Result: -0.0936996268
```

---

## Program 2: **Counting Zeros of ζ(s) via Argument Principle**

### Overview
This program determines the number of zeros of the Riemann zeta function on the critical line Re(s) = 1/2 over the interval [0, T] by leveraging the argument principle. It supports both adaptive and fixed-step integration methods to ensure precision and efficiency.

### Features
- Adjustable computational precision using arbitrary precision arithmetic.
- Warnings to indicate potential numerical inaccuracies.
- Outputs the total number of zeros in the specified range.

### Instructions
1. Run the script with an optional argument `T` specifying the upper bound of the interval. Defaults to `T = 30` if no argument is provided.
2. The program calculates and displays the number of zeros.

Example:
```
python NumberOfZerosOfZetaFunction.py 30
Расчёт нулей ζ(1/2 + it) для t ∈ [0, 30]
Приблизительное количество нулей: 49
```

---

## Program 3: **Locating Zeros of ζ(s) via the Dirichlet η-Function**

### Overview
This script identifies the zeros of the Riemann zeta function on the critical line Re(s) = 1/2 by leveraging the Dirichlet η-function and advanced numerical root-finding algorithms, including the golden-section search for refining approximations.

### Features
- Detects zeros within a user-defined range.
- Configurable step size for initial scanning and accuracy threshold for zero identification.
- Outputs a list of approximate zero locations with high precision.

### Instructions
1. Specify the range `[a, b]` for the imaginary part of `s` within the script.
2. Execute the program.
3. The zeros identified within the specified range are displayed.

Example:
```
Найденные нули на критической линии:
Номер 1: t ≈ 14.1347
Номер 2: t ≈ 21.0220
...
```

---

## Technical Dependencies
The scripts require the following Python libraries:
- `numpy`: Efficient numerical operations.
- `scipy`: Linear algebra tools for solving equations.
- `mpmath`: High-precision arithmetic for critical calculations.
- `tqdm`: Visual progress bars for long-running computations (used in Program 2).

Install the dependencies using pip:
```
pip install numpy scipy mpmath tqdm
```

---

## Important Considerations
- These programs are intended for advanced academic research and may need optimization for computations involving extensive data ranges.
- Numerical precision is limited by the chosen settings and inherent rounding errors, particularly for large imaginary values of `s`.

---

## Licensing
The repository is open-source and available under the MIT License for unrestricted academic and research use.

