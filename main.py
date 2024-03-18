from sympy import symbols, exp, diff
from sympy import solve
from sympy import Matrix

x, y = symbols('x y')

f1 = x * (y**2) - 1
f2 = y - exp(x)


u = x * (y**2) - 1
v = y - exp(x)

def f_1(x, y):
    return x * (y**2) - 1

def f_2(x, y):
    return y - exp(x)

def simple_iterations():
    df1_dx = diff(f1, x)
    df1_dy = diff(f1, y)

    df2_dx = diff(f2, x)
    df2_dy = diff(f2, y)

    print(f"f1 = {f1}")
    print(f"f2 = {f2}")

    print(f"Часткові похідні від f1: ")
    print("df1_dx = ", df1_dx)
    print("df1_dy = ", df1_dy)

    print(f"Часткові похідні від f2: ")
    print("df2_dx = ", df2_dx)
    print("df2_dy = ", df2_dy)

    x_val = 0.5
    y_val = 1.5
    df1_dx_val = df1_dx.subs({x:x_val, y:y_val})
    df1_dy_val = df1_dy.subs({x:x_val, y:y_val})
    df2_dx_val = df2_dx.subs({x:x_val, y:y_val})
    df2_dy_val = df2_dy.subs({x:x_val, y:y_val})

    alpha, betta, gamma, delta = symbols('alpha betta gamma delta')

    equations = [
        1 + alpha * df1_dx_val + betta * df2_dx_val,
        alpha * df1_dy_val + betta * df2_dy_val,
        1 + gamma * df1_dy_val + delta * df2_dy_val,
        gamma * df1_dx_val + delta * df2_dx_val
    ]
    solutions = solve(equations, (alpha, betta, gamma, delta))
    print("Корені системи рівнянь:")
    for var, val in solutions.items():
        rounded_val = round(val.evalf(), 4)
        print(f"{var} = {rounded_val}")

    x_prev = 0.5
    y_prev = 1.5

    alpha_val = solutions[alpha].evalf()
    betta_val = solutions[betta].evalf()
    gamma_val = solutions[gamma].evalf()
    delta_val = solutions[delta].evalf()

    for i in range(1, 11):
        print("________________")
        x_i = (x_prev + alpha_val * f_1(x_prev, y_prev) + betta_val * f_2(x_prev, y_prev)).evalf()
        y_i = (y_prev + gamma_val * f_1(x_prev, y_prev) + delta_val * f_2(x_prev, y_prev)).evalf()

        print(f"x_{i} = {x_i}")
        print(f"y_{i} = {y_i}")

        x_prev = x_i
        y_prev = y_i

    print(f"\nПохибка обчислення x = {f_1(x_prev, y_prev)}"
          f"\nПохибка обчислення y = {f_2(x_prev, y_prev)}")

def newton_method():
    print(f"u(x;y) = {u} \n v(x;y) = {v}")

    du_dx = diff(u, x)
    du_dy = diff(u, y)

    dv_dx = diff(v, x)
    dv_dy = diff(v, y)

    print("Часткові похідні від u(x;y): ")
    print(f"du_dx = {du_dx}")
    print(f"du_dy = {du_dy}")

    print("Часткові похідні від v(x;y): ")
    print(f"dv_dx = {dv_dx}")
    print(f"dv_dy = {dv_dy}")

    delta_k = Matrix([
        [du_dx, du_dy],
        [dv_dx, dv_dy]
    ])
    delta_k_x = Matrix([
        [u, du_dy],
        [v, dv_dy]
    ])
    delta_k_y = Matrix([
        [du_dx, u],
        [dv_dx, v]
    ])

    determ_delta_k = delta_k.det()
    determ_delta_k_x = delta_k_x.det()
    determ_delta_k_y = delta_k_y.det()

    print("Детермінант матриці delta_k:")
    print(determ_delta_k)
    print("Детермінант матриці delta_k_x:")
    print(determ_delta_k_x)
    print("Детермінант матриці delta_k_y:")
    print(determ_delta_k_y)

    x_prev = 0.5
    y_prev = 1.5
    e = 0.0001

    s_n = float('inf')
    n = 0
    while s_n > e:
        n += 1
        print("_____________")
        determ_delta_k_value = determ_delta_k.subs({x: x_prev, y: y_prev}).evalf()
        determ_delta_k_x_value = determ_delta_k_x.subs({x: x_prev, y: y_prev}).evalf()
        determ_delta_k_y_value = determ_delta_k_y.subs({x: x_prev, y: y_prev}).evalf()
        s_n = 1 / abs(determ_delta_k_value) * ((determ_delta_k_x_value ** 2 + determ_delta_k_y_value ** 2) ** 0.5)

        print("Детермінант матриці delta_k:")
        print(round(determ_delta_k_value, 4))

        print("Детермінант матриці delta_k_x:")
        print(round(determ_delta_k_x_value, 4))

        print("Детермінант матриці delta_k_y:")
        print(round(determ_delta_k_y_value, 4))

        x_k = x_prev - determ_delta_k_x_value / determ_delta_k_value
        y_k = y_prev - determ_delta_k_y_value / determ_delta_k_value

        print(f"\nx_{n} = {x_k} \ny_{n} = {y_k}")
        x_prev = x_k
        y_prev = y_k

    print(f"\nПохибка обчислення x = {f_1(x_prev, y_prev)}"
        f"\nПохибка обчислення y = {f_2(x_prev, y_prev)}")
def main():
    print("Комп'ютерний практикум №3 \nВаріант №11 \nРозумняк Руслан ")
    print("\n________Метод простих ітерацій________")
    simple_iterations()
    print("\n________Метод Ньютона________")
    newton_method()
if __name__ == "__main__":
    main()