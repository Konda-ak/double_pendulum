# Вывод уравнений движения двойного маятника

От лагранжиана до системы дифференциальных уравнений первого порядка, с вязким
трением (функция Рэлея) и внешним моментом в точке опоры.

## 1. Модель

Два точечных груза $m_1$, $m_2$ на невесомых жёстких стержнях длиной $L_1$, $L_2$.
Первый стержень закреплён в точке опоры $O$ (шарнир), второй — на первом грузе.
Обобщённые координаты — углы от вертикали вниз:

$$
q = (\theta_1,\ \theta_2)^T, \qquad \omega_i = \dot\theta_i .
$$

Декартовы координаты грузов (ось $y$ вверх, начало в точке опоры):

$$
\begin{aligned}
x_1 &= L_1 \sin\theta_1, & y_1 &= -L_1 \cos\theta_1,\\
x_2 &= L_1 \sin\theta_1 + L_2 \sin\theta_2, & y_2 &= -L_1 \cos\theta_1 - L_2 \cos\theta_2 .
\end{aligned}
$$

## 2. Кинетическая и потенциальная энергия

Скорости грузов:

$$
\begin{aligned}
\dot x_1 &= L_1 \cos\theta_1\,\omega_1, & \dot y_1 &= L_1 \sin\theta_1\,\omega_1,\\
\dot x_2 &= \dot x_1 + L_2 \cos\theta_2\,\omega_2, & \dot y_2 &= \dot y_1 + L_2 \sin\theta_2\,\omega_2 .
\end{aligned}
$$

Квадрат скорости второго груза:

$$
v_2^2 = L_1^2\omega_1^2 + L_2^2\omega_2^2 + 2L_1L_2\,\omega_1\omega_2(\cos\theta_1\cos\theta_2 + \sin\theta_1\sin\theta_2)
= L_1^2\omega_1^2 + L_2^2\omega_2^2 + 2L_1L_2\,\omega_1\omega_2\cos(\theta_1-\theta_2).
$$

Кинетическая энергия $T = \tfrac12 m_1 v_1^2 + \tfrac12 m_2 v_2^2$:

$$
T = \tfrac12 (m_1+m_2) L_1^2 \omega_1^2 + \tfrac12 m_2 L_2^2 \omega_2^2
  + m_2 L_1 L_2\, \omega_1 \omega_2 \cos(\theta_1-\theta_2).
$$

Потенциальная энергия $V = m_1 g y_1 + m_2 g y_2$:

$$
V = -(m_1+m_2)\, g L_1 \cos\theta_1 - m_2\, g L_2 \cos\theta_2 .
$$

Нуль потенциала выбран на уровне точки опоры. Он ни на что не влияет в уравнениях
движения, поэтому для оценки ошибок энергии используется не $E(0)$, а характерная
энергия системы

$$
E_{scale} = g\,\bigl(m_1 L_1 + m_2 (L_1 + L_2)\bigr),
$$

то есть работа подъёма обоих грузов из нижнего положения на уровень опоры.

## 3. Лагранжиан

$$
\mathcal{L} = T - V .
$$

## 4. Трение: функция Рэлея

Вязкое трение в точке опоры (коэффициент $c_1$) и в среднем шарнире (коэффициент $c_2$,
сила пропорциональна **относительной** угловой скорости стержней):

$$
R = \tfrac12 c_1\, \omega_1^2 + \tfrac12 c_2\, (\omega_2 - \omega_1)^2 .
$$

Обобщённая сила трения для координаты $q_i$ равна $-\partial R/\partial\dot q_i$.
Внешний момент $\tau$ приложен в точке опоры, то есть работает на координате $\theta_1$:

$$
Q_1 = \tau, \qquad Q_2 = 0 .
$$

Сейчас $\tau = 0$. Параметр оставлен в `derivatives(...)` для будущего регулятора.

## 5. Уравнения Эйлера–Лагранжа

$$
\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot q_i}
- \frac{\partial \mathcal{L}}{\partial q_i}
+ \frac{\partial R}{\partial \dot q_i} = Q_i, \qquad i = 1, 2 .
$$

Обозначим $\delta = \theta_1 - \theta_2$.

**Уравнение для $\theta_1$.** Импульс и производные:

$$
\frac{\partial \mathcal{L}}{\partial \omega_1} = (m_1+m_2)L_1^2\omega_1 + m_2L_1L_2\,\omega_2\cos\delta,
$$

$$
\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \omega_1} = (m_1+m_2)L_1^2\dot\omega_1 + m_2L_1L_2\cos\delta\,\dot\omega_2
- m_2L_1L_2\,\omega_2(\omega_1-\omega_2)\sin\delta,
$$

$$
\frac{\partial \mathcal{L}}{\partial \theta_1} = -m_2L_1L_2\,\omega_1\omega_2\sin\delta - (m_1+m_2)\,gL_1\sin\theta_1,
\qquad
\frac{\partial R}{\partial \omega_1} = c_1\omega_1 - c_2(\omega_2-\omega_1).
$$

Подставляем и сокращаем: $-\omega_2(\omega_1-\omega_2) + \omega_1\omega_2 = \omega_2^2$.

$$
(m_1+m_2)L_1^2\,\dot\omega_1 + m_2L_1L_2\cos\delta\,\dot\omega_2 + m_2L_1L_2\,\omega_2^2\sin\delta
+ (m_1+m_2)\,gL_1\sin\theta_1 + c_1\omega_1 - c_2(\omega_2-\omega_1) = \tau .
$$

**Уравнение для $\theta_2$.**

$$
\frac{\partial \mathcal{L}}{\partial \omega_2} = m_2L_2^2\omega_2 + m_2L_1L_2\,\omega_1\cos\delta,
$$

$$
\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \omega_2} = m_2L_2^2\dot\omega_2 + m_2L_1L_2\cos\delta\,\dot\omega_1
- m_2L_1L_2\,\omega_1(\omega_1-\omega_2)\sin\delta,
$$

$$
\frac{\partial \mathcal{L}}{\partial \theta_2} = m_2L_1L_2\,\omega_1\omega_2\sin\delta - m_2\,gL_2\sin\theta_2,
\qquad
\frac{\partial R}{\partial \omega_2} = c_2(\omega_2-\omega_1).
$$

Сокращаем: $-\omega_1(\omega_1-\omega_2) - \omega_1\omega_2 = -\omega_1^2$.

$$
m_2L_2^2\,\dot\omega_2 + m_2L_1L_2\cos\delta\,\dot\omega_1 - m_2L_1L_2\,\omega_1^2\sin\delta
+ m_2\,gL_2\sin\theta_2 + c_2(\omega_2-\omega_1) = 0 .
$$

## 6. Матричная форма

Оба уравнения линейны по ускорениям:

$$
M(\theta)\begin{pmatrix}\dot\omega_1\\ \dot\omega_2\end{pmatrix} = \begin{pmatrix} f_1 \\ f_2 \end{pmatrix},
\qquad
M(\theta) = \begin{pmatrix}
(m_1+m_2)L_1^2 & m_2L_1L_2\cos\delta\\
m_2L_1L_2\cos\delta & m_2L_2^2
\end{pmatrix},
$$

$$
\begin{aligned}
f_1 &= -m_2L_1L_2\,\omega_2^2\sin\delta - (m_1+m_2)\,gL_1\sin\theta_1 - c_1\omega_1 + c_2(\omega_2-\omega_1) + \tau,\\
f_2 &= \ \ m_2L_1L_2\,\omega_1^2\sin\delta - m_2\,gL_2\sin\theta_2 - c_2(\omega_2-\omega_1).
\end{aligned}
$$

Определитель матрицы масс не обращается в ноль:

$$
\det M = m_2 L_1^2 L_2^2\,(m_1 + m_2\sin^2\delta) > 0,
$$

поэтому $M$ обратима при любых углах:

$$
M^{-1} = \frac{1}{\det M}\begin{pmatrix}
m_2L_2^2 & -m_2L_1L_2\cos\delta\\
-m_2L_1L_2\cos\delta & (m_1+m_2)L_1^2
\end{pmatrix}.
$$

## 7. Система первого порядка

Вектор состояния $x = (\theta_1, \omega_1, \theta_2, \omega_2)^T$:

$$
\dot x = F(x, t) = \begin{pmatrix}
\omega_1\\
a_1\\
\omega_2\\
a_2
\end{pmatrix},
\qquad
\begin{pmatrix}a_1\\ a_2\end{pmatrix} = M^{-1}(\theta)\begin{pmatrix}f_1\\ f_2\end{pmatrix}.
$$

Эта функция реализована как `derivatives(state, t, params)` в `physics.py`
(вместо явной формулы с обратной матрицей используется `np.linalg.solve`).
Интегратору она передаётся как $f(y, t)$ через `make_rhs(params)`.

## 8. Баланс энергии

Умножим уравнение Эйлера–Лагранжа $i$ на $\dot q_i$ и сложим. Так как $\mathcal{L}$ не
зависит от времени явно, а $R$ — квадратичная форма по скоростям
($\sum_i \dot q_i\,\partial R/\partial\dot q_i = 2R$), получаем

$$
\frac{dE}{dt} = -2R + \tau\,\omega_1,
\qquad E = T + V .
$$

Мощность сил трения (при $\tau = 0$):

$$
P = \frac{dE}{dt} = -2R = -c_1\omega_1^2 - c_2(\omega_2-\omega_1)^2 \le 0 .
$$

Без трения и без момента энергия сохраняется.

## 9. Линеаризация около нижнего положения равновесия

Нижнее положение $\theta_1=\theta_2=0$, $\omega_1=\omega_2=0$. При малых углах
$\sin\theta\approx\theta$, $\cos\delta\approx 1$, квадратичные по скоростям члены
отбрасываются:

$$
M_0\ddot\theta + C\dot\theta + K\theta = \begin{pmatrix}\tau\\0\end{pmatrix},
$$

$$
M_0 = \begin{pmatrix}(m_1+m_2)L_1^2 & m_2L_1L_2\\ m_2L_1L_2 & m_2L_2^2\end{pmatrix},\quad
K = \begin{pmatrix}(m_1+m_2)gL_1 & 0\\ 0 & m_2gL_2\end{pmatrix},\quad
C = \begin{pmatrix}c_1+c_2 & -c_2\\ -c_2 & c_2\end{pmatrix}.
$$

В форме пространства состояний $\dot x = Ax + Bu$ матрица $A$ — якобиан $F$ в точке
равновесия, $u = \tau$. Именно эта форма нужна для LQR. В коде: `linearized_system(p)`
(аналитически: $A = \begin{pmatrix}0&I\\-M_0^{-1}K&-M_0^{-1}C\end{pmatrix}$ в перестановке состояний,
$B$ — столбец $M_0^{-1}(1,0)^T$) и `numeric_jacobian(p, x0)` (центральные разности в любой
точке, пригодится для верхнего положения).

## 10. Нормальные моды

Без трения ($C=0$) ищем решение $\theta = v\cos(\omega t)$:

$$
(K - \omega^2 M_0)\,v = 0 \quad\Longrightarrow\quad \det(K - \omega^2 M_0) = 0 .
$$

Раскрывая определитель и обозначив $a = m_1 + m_2$, получаем квадратное уравнение
относительно $\omega^2$:

$$
m_1L_1L_2\,\omega^4 - a\,g\,(L_1+L_2)\,\omega^2 + a\,g^2 = 0,
$$

$$
\omega_{\pm}^2 = \frac{g\left[a(L_1+L_2) \pm \sqrt{a^2(L_1+L_2)^2 - 4\,a\,m_1L_1L_2}\right]}{2\,m_1L_1L_2}.
$$

Форма моды из первой строки системы:

$$
\frac{v_2}{v_1} = \frac{a\,(g - \omega^2 L_1)}{m_2 L_2\,\omega^2}.
$$

Для $m_1=m_2=m$, $L_1=L_2=L$:

$$
\omega_{\pm}^2 = (2\pm\sqrt2)\,\frac{g}{L}, \qquad \frac{\theta_2}{\theta_1} = \pm\sqrt2 ,
\qquad T_{\pm} = \frac{2\pi}{\omega_{\pm}} .
$$

В первой моде ($\omega_-$) стержни качаются в фазе, во второй ($\omega_+$) — в противофазе,
причём второй стержень отклоняется сильнее первого. Если задать начальные углы в
пропорции $v_2/v_1$ одной моды и нулевые скорости, маятник колеблется с одной частотой
(пока углы малы).

## 11. Где что проверяется

| Формула | Проверка в `run_checks.py` |
|---|---|
| разделы 2–5: $T$, $V$, $R$, уравнения Эйлера–Лагранжа | `check_euler_lagrange`: $T$, $V$, $R$ записаны независимо, через декартовы координаты, производные берутся численно; невязка ≈ $10^{-8}$ |
| раздел 6–7: `derivatives` | то же |
| раздел 8: $dE/dt = -2R$ | `check_friction`: мгновенно и через интеграл мощности |
| раздел 9: $A$, $B$ | `check_linearization`: численный якобиан против аналитических $A$, $B$ |
| раздел 10: моды | `check_normal_modes`: периоды и формы мод для двух наборов параметров |
