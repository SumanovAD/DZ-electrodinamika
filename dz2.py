import urllib.request
import math
import os
import toml
import scipy.special as scp
import matplotlib.pyplot as plt
import requests

def load_variant_parameters(url, target_variant):
    """Загрузка параметров из TOML-файла по ссылке"""
    response = requests.get(url)
    response.raise_for_status()
    data = toml.loads(response.text)
    for entry in data.get("data", []):
        if entry.get("variant") == target_variant:
            D = float(entry["D"])
            fmin = float(entry["fmin"])
            fmax = float(entry["fmax"])
            return D, fmin, fmax
    raise ValueError(f"Вариант {target_variant} не найден в файле.")

if __name__ == '__main__':
    class Rcs:
        def __init__(self, r, fmin, fmax):
            self.r = r
            self.fmin = fmin
            self.fmax = fmax
            self.mass_f = []
            self.mass_rcs = []

        def calc(self):
            f = self.fmin
            def hn(n, x):
                return complex(scp.spherical_jn(n, x), scp.spherical_yn(n, x))
            while f <= self.fmax:
                l = 3e8 / f
                k = 2 * math.pi / l
                s = 0
                for n in range(1, 20):
                    an = scp.spherical_jn(n, k * self.r) / hn(n, k * self.r)
                    bn = (k * self.r * scp.spherical_jn(n - 1, k * self.r) - n * scp.spherical_jn(n, k * self.r)) / \
                         (k * self.r * hn(n - 1, k * self.r) - n * hn(n, k * self.r))
                    s += ((-1) ** n) * (n + 0.5) * (bn - an)
                self.mass_f.append(f)
                self.mass_rcs.append((l ** 2) * (abs(s) ** 2) / math.pi)
                f += 1e6
            plt.plot(self.mass_f, self.mass_rcs)
            plt.xlabel("f, Гц")
            plt.ylabel("ЭПР, м^2")
            plt.grid()
            plt.show()

    class Output:
        def __init__(self, mass_f, mass_rcs):
            self.mass_f = mass_f
            self.mass_rcs = mass_rcs

        def output(self):
            if not os.path.isdir("results"):
                os.mkdir("results")
            os.chdir("results")
            
            # Запись в CSV файл
            with open('results.csv', 'w') as f:
                f.write('№,Частота (Гц),ЭПР (м^2)\n')  # Заголовок
                for i, (freq, rcs) in enumerate(zip(self.mass_f, self.mass_rcs), 1):
                    f.write(f'{i},{freq},{rcs}\n')

    toml_url = "https://jenyay.net/uploads/Student/Modelling/task_rcs_01.toml"
    variant = 6

    D, fmin, fmax = load_variant_parameters(toml_url, variant)
    r = D / 2

    sph = Rcs(r, fmin, fmax)
    sph.calc()

    o = Output(sph.mass_f, sph.mass_rcs)
    o.output()
