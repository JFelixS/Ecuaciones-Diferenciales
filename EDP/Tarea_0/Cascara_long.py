# -*- coding: utf-8 -*- 
## Este primer renglon nos da la codificacion con la cual
## vamos a estar trabajando.


## Curso: Ecuaciones Diferenciale Parciales
## Profesores:
## Oscar Reula & Manuel Tiglio
## Fecha: 27 de Agosto de 2020
## Tarea 0
## Solución correspondiente a dos hemisferios a potenciales opuestos


##  Aquí resolveremos numéricamente el problema de una cáscara con sus hemisferios a potenciales opuestos. La solución analítica se la puede encontrar, por ejemplo en: http://www.famaf.unc.edu.ar/~reula/Docencia/Electromagnetismo/electrodynamics.pdf, página 81.Debido a la simetría axial solución se puede expresar como suma de polinomios de Legendre P_n​(x), donde x=cos(θ) y potencias de r, la coordenada radial.
## \phi(x,r) = \\sum_{n=0}^{\\infty} A_n r^{-(n+1)}P_n(x).

## Supondremos la cáscara tiene radio a y reescalaremos la variable radia como r→a/r​ de esa manera queda adimensional y la cáscara a radio unidad.

## Debemos calcular ahora los coeficientes $A_n$. Notemos que $A_{2n}=0, \\forall n$ ya que,
## la cáscara tiene simetría impar ante cambios $x \to -x$ y por lo tanto ninguno de los polinomios pares puede contribuir. Los impares se calculan imponiendo la condición de contorno: 
## \phi(x,1) = V_0 \theta(x)
## donde $\theta$ es la función escalón.",
## Tomaremos $V_0=1$ en lo que sigue.,
 ## Siempre se puede reescalear la solución al final,


## Tenemos así que (recordando que la cáscara está a $r=1$),
## \theta(x) = \\sum_{n=0}^{\\infty} A_n P_n(x),
## multiplicando por $P_{m}(x)$ ambos lados e integrando en $[-1,1]$ obtenemos,
## \int_{-1}^{1} P_n(x) \\theta(x) dx = \\frac{2}{2n+1} A_n.

## Haciendo solo as integrales entre $[0,1]$ y multiplicando or $2$, (ya que las funciones son impares) obtenemos,
## A_n = (2n+1) \\int_0^1 P_n(x) dx \\;\\;\\; n \\; \\mbox{impar}.

## Haremos estos cálculos analíticamente y luego graficaremos la solución.
## Para ilustrar el caso donde es conveniente hacer las integrales en forma numérica consideraremos el caso,
## \phi(x,1) = V_0 \\theta(x) e^{-x^2}.
## Los coeficientes vendrán dados por:
## \theta(x) e^{-x^2} = \\sum_{n=0}^{\\infty} B_n P_n(x),
## B_n = (2n+1) \\int_0^1 P_n(x) e^{-x^2} dx \\;\\;\\; n \\; \\mbox{impar}.

## A estos los calcularemos numéricamente.


##  Primero definimos todo lo necesario de las librerías de Python.

#matplotlib inline
#matplotlib notebook
#from sympy import *
from sympy import simplify, diff, integrate, Integral, Sum, lambdify, legendre
import sympy.functions as sp 
import numpy as np
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import animation, rc
from IPython.display import HTML
#init_printing(use_unicode=True)
#x, y, z, u, v, r = symbols('x y z, u, v, r')
from sympy.abc import x, y, r, n, k, m
#k, m, n = symbols('k m n', integer=True)
#f, step, potential = symbols('f step potential', cls=Function)
#var('n m x')
#theta = Symbol("theta")
#phi = Symbol("phi")

## Luego probamos que las definiciones funcionan correctamente con algunos ejemplos:

sp.legendre(3,x)   #sp se refiere a simpy.functions
print(sp.legendre(3,x))

## Podemos por ejemplo integrarlos (o diferenciarlos, o hacer otras operaciones con ellos), por ejemplo:

integrate(sp.legendre(5,x),x) #integrate es una función que importamos desde sympy
print(integrate(sp.legendre(5,x),x))

##  Por ejemplo podemos probar una de las fórmulas más útiles para hacer integrales usando los polinomios de Legendre:dP_(l+1)/dx​​ − dP_(l−1)/dx​​ − l(l+1)P_l ​= 0.
## Solo se puede constatar la fórmula para valores particulares (no la sabe simplificar completamente). Poner cualquier valor entero positivo (>1) para k.


k=5
diff(sp.legendre(k+1,x),x) - diff(sp.legendre(k-1,x),x) - (2*k+1)*sp.legendre(k,x)
#diff es una función importada desde sympy
print(diff(sp.legendre(k+1,x),x) - diff(sp.legendre(k-1,x),x) - (2*k+1)*sp.legendre(k,x))


## Grafiquemos ahora algunos polinomios de Legendre. Para ello primero los transformamos en funciones propiamente dichas (y no símbolos) y lo hacemos adaptados para ser graficados con las facilidades numéricas the Python (numpy).

x_vals=np.linspace(-1,1,200)
f = lambdify(x, sp.legendre(8,x),"numpy") #convierte la expresión en una función llamable que acepta numpy arrays
l_vals=f(x_vals)
mpl.plot(x_vals,l_vals)
## mpl.show() Por el momento no la visualizamos, pero
# abajo en la linea 130 mostramos las graficas.


## Primero hacemos el caso analítico: 

## Las integrales se pueden hacer o bien analíticamente o bien numéricamente. En el caso que sigue tomamos una grilla de L puntos e integramos numéricamente. Esto siempre funciona, pero hay que tener cuidado que la aproximación sea buena, en particular para polinomios de orden alto hay que tomar intervalos de grilla bien pequeños.

L = 200 # números de puntos de grilla, ajustar al grado del polinomio, se necesitan al menos 10 puntos por cada l.
L = 100

Integral(sp.legendre(3,x),(x,0,1)).as_sum(L).evalf()
#as_sum(L) hace una integral por trapecios con L puntos, pero la expresa como una fracción
#evalf(N) evalúa la expresión con N dígitos significativos (8 por defecto)
print(Integral(sp.legendre(3,x),(x,0,1)).as_sum(L).evalf())


## También en algunos casos las integrales se pueden hacer en forma exacta.

#Integral(sp.legendre(3,x),(x,0,1)).doit().evalf() #para expresar como número decimal

Integral(sp.legendre(3,x),(x,0,1)).doit()         #para expresar como fracción
print(Integral(sp.legendre(3,x),(x,0,1)).doit())

##  Primero hacemos un cálculo algebraico, en el caso que el potencial viene dado por θ(x). Calculamos los coeficientes (solo para los n impares, n=2m+1).



A=lambda n: (2*(2*n+1)+1)*Integral(sp.legendre(2*n+1,x),(x,0,1)).doit()
A(3)
print(A(3))


A(15).evalf()
print(A(15).evalf())


## Contruimos el potencial en la frontera (r=1) para saber que las condiciones de contorno se satisfacen aproximadamente.

N=15 # número de términos en la suma
#step_approx=lambdify(x,Sum(A(m).doit()*sp.legendre(2*m+1,x),(m,0,N)).doit(),"numpy")
step_approx=lambdify(x,Sum(A(m).evalf()*sp.legendre(2*m+1,x),(m,0,N)).doit(),"numpy")


#Sum es una función de sympy. Recordar que m es un símbolo que también importamos desde simpy.
#Un ejemplo para entender la expresión de arriba es
from sympy.abc import p    #notar que ya importamos m antes
Sum(p, (p,0,3))
print(Sum(p, (p,0,3)))


## Graficamos la aproximación y comparamos con la solución exacta. Note que si ponemos muchos términos la solución no converge.

fig = mpl.figure(figsize = (10,10))
step = lambdify(x,sp.sign(x),"numpy")
f_vals = step(x_vals)
mpl.plot(x_vals, f_vals, color = "k", label = "Exact")
for N in [5,15,25]:
    x_vals = np.linspace(-1, 1, 200)
    step_approx=lambdify(x,Sum(A(m).evalf()*sp.legendre(2*m+1,x),(m,0,N)).doit(),"numpy")
    z_vals = step_approx(x_vals)
    mpl.plot(x_vals, z_vals, label = f"N = {N:d}")
    mpl.ylabel("step")

mpl.ylim(-2,2) #seteamos los límites del ploteo
mpl.legend() #para ver los labels
# mpl.show()


##  Generamos ahora el potencial completo:

N = 15
potential=lambdify((x,r),Sum(A(m).evalf()*sp.legendre(2*m+1,x)*r**(-2*(m+1)),(m,0,N)).doit(),"numpy")

## Graficamos la aproximación a la función escalón que calculamos:

y_vals = np.linspace(1,10,200)
X,Y = np.meshgrid(x_vals,y_vals)
z_vals = potential(X,Y)
fig = mpl.figure(figsize = (12,10))
ax = fig.gca(projection='3d')
ax.set_zlim(-1.2, 1.2)
surf = ax.plot_surface(X,Y,z_vals, rstride=2, alpha=0.9, cstride=2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
levels = np.arange(-1., 1., 0.2)
levels_y = np.arange(0.9, 2., 0.1)
cset = ax.contour(X, Y, z_vals, levels, zdir='z', offset=1.5, cmap=cm.coolwarm)
cset = ax.contour(X, Y, z_vals, levels, zdir='x', offset=1.4, cmap=cm.coolwarm)
cset = ax.contour(X, Y, z_vals, levels_y, zdir='y', offset=-2, cmap=cm.coolwarm)
fig.colorbar(surf, shrink=0.5, aspect=5)
# mpl.show()

## Ahora calcularemos el segundo caso, V(x) = \theta(x)*e^(−x^2). Primero lo hacemos con un número pequeño de puntos en la integración numérica (20).

B_20=lambda m: (2*(2*m+1)+1)*Integral(sp.legendre(2*m+1,x)*sp.exp(-x*x),(x,0,1)).as_sum(20)

## Verificamos que funciona

B_20(3).evalf()
print(B_20(3).evalf())


## Ahora ponemos los números en una tupla para poderlos graficar de forma simple.

E_20=[]
for j in range(10):
        E_20.append(B_20(j))

## Repetimos con más puntos para ver la diferencia en la integración numérica:

#L = 100
print(L)


B_40=lambda m: (2*(2*m+1)+1)*Integral(sp.legendre(2*m+1,x)*sp.exp(-x*x),(x,0,1)).as_sum(40)

E_40=[]
for j in range(10):
        E_40.append(B_40(j))


mpl.plot(range(10),E_40,'bo', label = '40 puntos')
mpl.plot(range(10),E_20, 'r^', label = '20 puntos')
mpl.legend()
mpl.xlabel('m')
mpl.ylabel('Integral numerica de B(m)')
## mpl.show()


## Vemos que se pierde precisión a medida que crece el orden del polinomio como era de esperarse. Ahora ya estamos en condiciones de proseguir. Chequeamos que la condición de contorno se satisfaga. Para ello sumamos los 20 primeros términos de la suma. (Vea que sucede si toma más términos).

N=20 # número de términos en la suma
step=lambdify(x,Sum(B_40(m).doit()*sp.legendre(2*m+1,x),(m,0,N)).doit(),"numpy")


## Graficamos la aproximación a la función escalón que calculamos:


x_vals = np.linspace(-1, 1, 200)
z_vals = step(x_vals)
f = lambdify(x,sp.sign(x)*sp.exp(-x*x),"numpy")
f_vals = f(x_vals)
mpl.plot(x_vals, z_vals, x_vals, f_vals)
mpl.ylabel("step")
mpl.xlabel("x")
# mpl.show()

## Ahora calculamos una aproximación a la solución completa, también tomando 20 términos de la serie.

potential=lambdify((x,r),Sum(B_40(m).doit()*sp.legendre(2*m+1,x)*r**(-2*(m+1)),(m,0,N)).doit(),"numpy") 

## Finalmente la graficamos:

y_vals = np.linspace(1,10,200)
X,Y = np.meshgrid(x_vals,y_vals)
z_vals = potential(X,Y)
fig = mpl.figure(figsize = (12,10))
ax = fig.gca(projection='3d')
ax.set_zlim(-1.2, 1.2)
surf = ax.plot_surface(X,Y,z_vals, rstride=2, alpha=0.9, cstride=2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
levels = np.arange(-1., 1., 0.2)
levels_y = np.arange(0.9, 2., 0.1)
cset = ax.contour(X, Y, z_vals, levels, zdir='z', offset=1.5, cmap=cm.coolwarm)
cset = ax.contour(X, Y, z_vals, levels, zdir='x', offset=1.4, cmap=cm.coolwarm)
cset = ax.contour(X, Y, z_vals, levels_y, zdir='y', offset=-2, cmap=cm.coolwarm)
fig.colorbar(surf, shrink=0.5, aspect=5)
mpl.show()





