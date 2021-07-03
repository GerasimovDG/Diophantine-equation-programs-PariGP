# Diophantine-equation-programs-PariGP
A set of programs for working with Diophantine equations. The programs are written in [Pari/GP](https://pari.math.u-bordeaux.fr/download.html).

Для работы всех программ необходим [Pari/GP](https://pari.math.u-bordeaux.fr/download.html) версии не ниже 2.13.1.

##### Содержание 
- [ABC=N.gp и ABC=N_only-answer.gp](#ABC=N)
- [cryptosystem.gp](#crypto)
- [generalizedPellDN.gp](#generalizedPellDN)
- [matrix_uniformity.gp](#matrix_uniformity)
- [meanNumberOfSolution.gp](#meanNumberOfSolution)

<a name="ABC=N"><h2>ABC=N.gp и ABC=N_only-answer.gp</h2></a>
Поиск примитивных и непримитивных решений уравнения 𝐴𝑥^2 + 𝐵𝑥𝑦 + 𝐶𝑦^2 = 𝑁, где Δ = 𝐵^2 − 4𝐴𝐶 > 0, 𝐴, 𝐵, 𝐶− целые числа, а 𝑁 - целое
число, неравное 0.

Все функции файла **ABC=N.gp** выводят подробности выполнения на экран, функции из **ABC=N_only-answer.gp** выводят на экран только результат вычислений.
### Список функций:
- [ABC_N(A,B,C,N)](#ABC-N)
- [validationInputData(A,B,C)](#validationInputData(A,B,C))
- [transformationABC(a,b,c,N)](#transformationABC(a,b,c,N))
### Описание функций:
1.  <a name="ABC-N"><h4><i>ABC_N(A,B,C,N)</i></h4></a> - находит примитивные и непримитивные решения уравнения 𝐴𝑥^2 + 𝐵𝑥𝑦 + 𝐶𝑦^2 = 𝑁, где Δ = 𝐵^2 − 4𝐴𝐶>0, 𝐴, 𝐵, 𝐶− целые числа, а 𝑁 - целое
число, неравное 0.
    + Аргументы функции:
      + 𝐴, 𝐵, 𝐶 - целые числа.
      + 𝑁 - целое число, 𝑁 ̸= 0.
    + Возвращаемое значение - массив [𝐿𝑖𝑠𝑡1, 𝐿𝑖𝑠𝑡2], где:
      + 𝐿𝑖𝑠𝑡1 - список примитивных решений.
      + 𝐿𝑖𝑠𝑡2 - список непримитивных решений.

2. <a name="validationInputData(A,B,C)"><h4><i>validationInputData(A,B,C)</i></h4></a> - внутренняя валидация входных данных для функции ***ABC_N(A,B,C,N)***. Возвращает 1 *(true)*, если валидация пройдена, 0 *(false)*, если нет.
4. <a name="transformationABC(a,b,c,N)"><h4><i>transformationABC(a,b,c,N)</i></h4></a> - производит унимодулярную замену переменных, которая преобразует уравнение 𝑎𝑥^2+𝑏𝑥𝑦+𝑐𝑦^2 = 𝑁 с НОД(𝑎, 𝑁) ̸= 1 в уравнение вида 𝐴𝑋^2 + 𝐵𝑋𝑌 + 𝐶𝑌^2 = 𝑁, где НОД(𝐴, 𝑁) = 1. Замена происходит с использованием матрицы
𝑇 = (︃ 𝛼, 𝛽; 𝛾, 𝛿 )︃.
    + Аргументы функции:
      + 𝑎, 𝑏, 𝑐 - целые числа.
      + 𝑁 - целое число, 𝑁 ̸= 0.
    + Возвращаемое значение: массив [𝐴, 𝐵, 𝐶, 𝛼, 𝛽 , 𝛾, 𝛿], где:
      + 𝐴, 𝐵, 𝐶 - аргументы нового уравнения 𝐴𝑋2 + 𝐵𝑋𝑌 + 𝐶𝑌 2 = 𝑁.
      + 𝛼, 𝛽 , 𝛾, 𝛿 - аргументы матрицы 𝑇.

<a name="crypto"><h2>cryptosystem.gp</h2></a>
...
<a name="generalizedPellDN"><h2>generalizedPellDN.gp</h2></a>
...
<a name="matrix_uniformity"><h2>matrix_uniformity.gp</h2></a>
...
<a name="meanNumberOfSolution"><h2>meanNumberOfSolution.gp</h2></a>
...

