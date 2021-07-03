# Diophantine-equation-programs-PariGP
A set of programs for working with Diophantine equations. The programs are written in [Pari/GP](https://pari.math.u-bordeaux.fr/download.html).

Для работы всех программ необходим [Pari/GP](https://pari.math.u-bordeaux.fr/download.html) версии не ниже 2.13.1.

## Содержание 
- [ABC=N.gp и ABC=N_only-answer.gp](#ABC=N)
- [cryptosystem.gp](#crypto)
- [generalizedPellDN.gp](#generalizedPellDN)
- [matrix_uniformity.gp](#matrix_uniformity)
- [meanNumberOfSolution.gp](#meanNumberOfSolution)
_____

<a name="ABC=N"><h2>ABC=N.gp и ABC=N_only-answer.gp</h2></a>
Поиск примитивных и непримитивных решений уравнения 𝐴𝑥^2 + 𝐵𝑥𝑦 + 𝐶𝑦^2 = 𝑁, где Δ = 𝐵^2 − 4𝐴𝐶 > 0, 𝐴, 𝐵, 𝐶− целые числа, а 𝑁 - целое
число, неравное 0.

Все функции файла **ABC=N.gp** выводят подробности выполнения на экран, функции из **ABC=N_only-answer.gp** выводят на экран только результат вычислений.
### Список функций:
- [ABC_N(A,B,C,N)](#ABC-N)
- [validationInputData(A,B,C)](#validationInputData(A,B,C))
- [searchTetta(A,B,C,N)](#searchTetta(A,B,C,N))
- [transformationABC(a,b,c,N)](#transformationABC(a,b,c,N))
### Описание функций:
1.  <a name="ABC-N"><h4><i>ABC_N(A,B,C,N)</i></h4></a> - находит примитивные и непримитивные решения уравнения 𝐴𝑥^2+𝐵𝑥𝑦+𝐶𝑦^2 = 𝑁, где Δ=𝐵^2−4𝐴𝐶>0, 𝐴, 𝐵, 𝐶− целые числа, а 𝑁 - целое
число, неравное 0.
    + Аргументы функции:
      + 𝐴, 𝐵, 𝐶 - целые числа.
      + 𝑁 - целое число, 𝑁 ≠ 0.
    + Возвращаемое значение - массив [𝐿𝑖𝑠𝑡1, 𝐿𝑖𝑠𝑡2], где:
      + 𝐿𝑖𝑠𝑡1 - список примитивных решений.
      + 𝐿𝑖𝑠𝑡2 - список непримитивных решений.

2. <a name="validationInputData(A,B,C)"><h4><i>validationInputData(A,B,C)</i></h4></a> - внутренняя валидация входных данных для функции ***ABC_N(A,B,C,N)***. Возвращает 1 *(true)*, если валидация пройдена, 0 *(false)*, если нет.
3. <a name="searchTetta(A,B,C,N)"><h4><i>searchTetta(A,B,C,N)</i></h4></a> - находит все 𝜃, удовлетворяющие сравнению 𝑎𝜃^2 + 𝑏𝜃 + 𝑐 ≡ 0 (mod |𝑁|).
    + Аргументы функции: 𝑎, 𝑏, 𝑐, 𝑁 - целые числа.
    + Возвращаемое значение: список подходящих 𝜃 или пустой список, если таких 𝜃 нет.
4. <a name="transformationABC(a,b,c,N)"><h4><i>transformationABC(a,b,c,N)</i></h4></a> - производит унимодулярную замену переменных, которая преобразует уравнение 𝑎𝑥^2+𝑏𝑥𝑦+𝑐𝑦^2 = 𝑁 с НОД(𝑎, 𝑁) ̸≠ 1 в уравнение вида 𝐴𝑋^2 + 𝐵𝑋𝑌 + 𝐶𝑌^2 = 𝑁, где НОД(𝐴, 𝑁) = 1. Замена происходит с использованием матрицы
𝑇 = (︃ 𝛼, 𝛽; 𝛾, 𝛿 )︃.
    + Аргументы функции:
      + 𝑎, 𝑏, 𝑐 - целые числа.
      + 𝑁 - целое число, 𝑁 ̸≠ 0.
    + Возвращаемое значение: массив [𝐴, 𝐵, 𝐶, 𝛼, 𝛽 , 𝛾, 𝛿], где:
      + 𝐴, 𝐵, 𝐶 - аргументы нового уравнения 𝐴𝑋2 + 𝐵𝑋𝑌 + 𝐶𝑌 2 = 𝑁.
      + 𝛼, 𝛽 , 𝛾, 𝛿 - аргументы матрицы 𝑇.
_____
<a name="crypto"><h2>cryptosystem.gp</h2></a> 
...
_____
<a name="generalizedPellDN"><h2>generalizedPellDN.gp</h2></a>
...
_____
<a name="matrix_uniformity"><h2>matrix_uniformity.gp</h2></a> - работа с матрицами второго порядка, имеющими неприводимый характеристический многочлен 𝜆^2 − 𝑑.
### Список функций:
- [is_matrix_uniformity(matrixA, matrixB)](#is_matrix_uniformity)
- [matrix_similarity_classes(d, flag = 0)](#matrix_similarity_classes)
- [matrix_similarity_classes_interval_d(dStart, dEnd)](#matrix_similarity_classes_interval_d)
- [mean_count_of_classes(file_name)](#mean_count_of_classes)
- [number_of_d_with_classes(file_name)](#number_of_d_with_classes)
### Описание функций:
1.  <a name="is_matrix_uniformity"><h4><i>is_matrix_uniformity(matrixA, matrixB)</i></h4></a> - Определяет подобны ли матрицы
𝑚𝑎𝑡𝑟𝑖𝑥𝐴 и 𝑚𝑎𝑡𝑟𝑖𝑥𝐵.
    + Аргументы функции: 𝑚𝑎𝑡𝑟𝑖𝑥𝐴, 𝑚𝑎𝑡𝑟𝑖𝑥𝐵 - матрицы второго порядка, имеющие неприводимый характеристический многочлен 𝜆^2 − 𝑑.
    + Возвращаемое значение: 1 - (True) если матрицы подобны, 0 (False) - иначе.
2.  <a name="matrix_similarity_classes"><h4><i>matrix_similarity_classes(d, flag = 0)</i></h4></a> - Находит количество классов подобия матриц с характеристическим многочленом 𝜆^2 − 𝑑. Опционально выводит представителей классов подобия.
    + Аргументы функции:
      + 𝑑
      + flag - если равен 1, выводит представителей классов подобия. По умолчанию равен 0.
    + Возвращаемое значение: массив [𝑐𝑜𝑢𝑛𝑡, {𝑙𝑖𝑠𝑡}], где:
      + count - количество классов подобия.
      + list - список представителей классов подобия (опционально. При 𝑓𝑙𝑎𝑔 = 1).
3. <a name="matrix_similarity_classes_interval_d"><h4><i>matrix_similarity_classes_interval_d(dStart, dEnd)</i></h4></a> - Вычисляет количество классов подобия матриц с характеристическим многочленом 𝜆^2 − 𝑑 для всех (либо только для простых) 𝑑 на отрезке [𝑑𝑆𝑡𝑎𝑟𝑡, 𝑑𝐸𝑛𝑑]. Результат вычислений записывается в файл ***"matrix_similarity_D=\<dStart\>-\<dEnd\>.txt"*** (и в файл ***"matrix_similarity_D=\<dStart\>-\<dEnd\>_with_classes.txt"***, если выбрано сохранение представителей классов подобия).
    + Аргументы функции:
      + 𝑑𝑆𝑡𝑎𝑟𝑡 - начальное значение для 𝑑.
      + 𝑑𝐸𝑛𝑑 - конечное значение для 𝑑.
4. <a name="mean_count_of_classes"><h4><i>mean_count_of_classes(file_name)</i></h4></a>  - Вычисляет среднее число классов подобия, используя таблицу с данными о количестве классов подобия.
    + Аргументы функции: 𝑓𝑖𝑙𝑒_𝑛𝑎𝑚𝑒 - имя файла* с таблицей вида: 

        | значение 𝑑  | количество классов подобия |
        | ------------- | ------------- |

    > \* - файл, полученный в результате выполнения функции ***"[matrix_similarity_classes_interval_d(dStart, dEnd)](#matrix_similarity_classes_interval_d)"***.
    
    Результат записывается в файл ***"\<file_name\>_average_count_of_classes.txt"***

5. <a name="number_of_d_with_classes"><h4><i>number_of_d_with_classes(file_name)</i></h4></a> - Вычисляет как часто встречаются опре-
деленное количество классов в файле с таблицей вида: значение 𝑑 - количество клас-
сов подобия.
    + Аргументы функции: 𝑓𝑖𝑙𝑒_𝑛𝑎𝑚𝑒 - имя файла* с таблицей вида: 
    
        | значение 𝑑  | количество классов подобия |
        | ------------- | ------------- |
    > \* - файл, полученный в результате выполнения функции ***"[matrix_similarity_classes_interval_d(dStart, dEnd)](#matrix_similarity_classes_interval_d)"***.
    
    Результат записывается в файл ***"<file_name>_number_classes_at_d.txt"***
_____
<a name="meanNumberOfSolution"><h2>meanNumberOfSolution.gp</h2></a>
...

