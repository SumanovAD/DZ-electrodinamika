import numpy as np
import math 
import tools


class Gaussian:
    
    #Источник, создающий гауссов импульс
   

    def __init__(self, magnitude, dg, wg):
        
        #magnitude - максимальное значение в источнике;
        #dg - коэффициент, задающий начальную задержку гауссова импульса;
        #wg - коэффициент, задающий ширину гауссова импульса.
       
        self.magnitude = magnitude
        self.dg = dg
        self.wg = wg

    def getE(self, time):
        return self.magnitude * np.exp(-((time - self.dg) / self.wg) ** 2)


#if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
W0 = 120.0 * np.pi

    # Число Куранта
Sc = 1.0

    #Скорость света
c=3e8

    # Параметры сигнала
eps1=7
fmin=0
fmax=2.5e9
f0=(fmax+fmin)/2
DeltaF=fmax-fmin
A_0=100
A_max=100
wg=2 * np.sqrt(np.log(A_max)) / (np.pi * fmax)
dg=wg * np.sqrt(np.log(A_0))

    #Шаг по пространству в метрах
dx=c/(fmax*20)/np.sqrt(eps1)
Nl=c/(f0*dx)

    #шаг во времени в секундах
dt=(Sc*dx)/c

    #Дискретизация параметров сигнала
wg=wg/dt
dg=dg/dt

    # Время расчета в секундах
maxTimesec = 3e-9
maxTime=int(np.ceil(maxTimesec/dt))

    # Размер области моделирования в метрах
maxSize_m = 3
maxSize=int(np.ceil(maxSize_m/dx))

    # Положение источника в метрах
sourcePos_m = 1
sourcePos=int(np.ceil(sourcePos_m/dx))

 # Датчики для регистрации поля
probesPos1 = [1.2]
probesPos=[int(np.ceil(probesPos1/dx))]
probes = [tools.Probe(int(np.ceil(pos/dx)), maxTime) for pos in probesPos1]

    # Где начинается поглощающий диэлектрик
layerloss = 0.1
layer_lossx = int(np.ceil(layerloss/dx))

    # Параметры среды
    # Диэлектрическая проницаемость
eps = np.ones(maxSize)
eps[:] = eps1

    # Магнитная проницаемость
mu = np.ones(maxSize - 1)

    # Потери в среде. loss = sigma * dt / (2 * eps * eps0)
loss = np.zeros(maxSize)
loss[layer_lossx:] = 0.01

    # Магнитные потери в среде. loss_m = sigma_m * dt / (2 * mu * mu0)
loss_m = np.zeros(maxSize - 1)
loss_m[:] = loss[:1]

    # Коэффициенты для расчета поля E
ceze = (1.0 - loss) / (1.0 + loss)
cezh = (Sc * W0) / (eps * (1.0 + loss))

    # Коэффициенты для расчета поля H
chyh = (1.0 - loss_m) / (1.0 + loss_m)
chye = Sc / (mu * W0 * (1.0 + loss_m))

Ez = np.zeros(maxSize)

Ez_left1 = 0  # Ez[1] на предыдущем шаге
Ez_left2 = 0  # Ez[1] на шаге до предыдущего
Ez0_old1 = 0  # Ez[0] на предыдущем шаге
Ez0_old2 = 0  # Ez[0] на шаге до предыдущего

Hy = np.zeros(maxSize - 1)
source = Gaussian(1,dg, wg)

    # Параметры отображения поля E
display_field = Ez
display_ylabel = 'Ez, V/m'
display_ymin = -1.1
display_ymax = 1.1

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
display = tools.AnimateFieldDisplay(maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel,dx,dt)

display.activate()
display.drawProbes(probesPos)
display.drawSources([sourcePos])

for q in range(1, maxTime):
        # Расчет компоненты поля H
    Hy[:] = chyh * Hy + chye * (Ez[1:] - Ez[:-1])

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
    Hy[sourcePos - 1] -= Sc / (W0 * mu[sourcePos - 1]) * source.getE( q)

        # Расчет компоненты поля E
    Ez[1:-1] = ceze[1: -1] * Ez[1:-1] + cezh[1: -1] * (Hy[1:] - Hy[:-1])
    # Левое граничное условие ABC второго порядка (ABC 2-го порядка)
    Ez[0] = Ez0_old2 + ((Sc - 1) / (Sc + 1)) * (Ez[1] - Ez_left2)

    # Обновление истории значений
    Ez_left2 = Ez_left1
    Ez_left1 = Ez[1]

    Ez0_old2 = Ez0_old1
    Ez0_old1 = Ez[0]


        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
    Ez[sourcePos] += (Sc / (np.sqrt(eps[sourcePos] * mu[sourcePos])) *
                          source.getE(q + 0.5))
    #PEC находится справа
    Ez[-1]=0

        # Регистрация поля в датчиках
    for probe in probes:
        probe.addData(Ez, Hy)

    if q % 10 == 0:
        display.updateData(display_field, q)

display.stop()

   # Отображение сигнала, сохраненного в датчиках
tools.showProbeSignals(probes, -1.1, 1.1,dx,dt,maxTime)
    #отображение спектра
tools.Spectrum(f0,DeltaF,wg,dg)
