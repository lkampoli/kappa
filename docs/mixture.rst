Mixture
*******

Смесь в поуровневом приближении

Смесь в поуровневом приближении, состоящая из произвольного числа молекул и атомов, описывается объектом класса Mixture (наследует от класса Approximation). Объект класса Mixture хранит информацию о всех взаимодействиях частиц смеси, позволяет рассчитывать удельные теплоемкости, числовые и массовые плотности, транспортные коэффиценты.

Примечание: так как объекты класса Mixture хранят информацию о матрицах, используемых для расчета коэффициентов переноса (для оптимизации скорости расчетов), они не являются thread-safe. Класс имеет следующие конструкторы:

Класс имеет следующие методы:
