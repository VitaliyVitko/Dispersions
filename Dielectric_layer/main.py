from sympy import *
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm


class Dispersion:
    # find root with Newton method
    def fr(self, variable, close_root):
        solve_root = I
        while (im(solve_root) > 0) or (im(solve_root) < 0):
            solve_root = nsolve(self, variable, close_root, verify=False)
            close_root = re(solve_root)
        return solve_root

    # rearrangement of massive if necessary
    def rearr(massive_2, massive_1):
        arr = []
        massive_2.reverse()
        massive_1.reverse()
        for p in range(len(massive_1)):
            arr.append([massive_2[p], massive_1[p]])
        return arr

    def calculation(self, variable, wavenumber_massive, parameters, close_root, name):
        result_massive = []  # result massive for save
        dc = []  # result massive for plot
        with tqdm(total=len(wavenumber_massive)) as progressbar:  # add active bar for total number of points in calculation
            for p in range(len(wavenumber_massive)):
                progressbar.update(1)
                if (len(result_massive) == 0):
                    r0 = close_root
                else:
                    if (len(result_massive)==1):
                        r0 = result_massive[p - 1][1]
                    else:
                        r0=2*result_massive[p - 1][1]-result_massive[p - 2][1]
                solved_root = Dispersion.fr(self.subs(k, wavenumber_massive[p]).subs(parameters).subs(paramb), variable, r0);
                result_massive.append([wavenumber_massive[p], solved_root])
                dc.append(solved_root / (2 * pi).evalf())
                np.savetxt('disp_' + str(name) + '.csv', result_massive, delimiter=" ,");  # save results for ORIGIN PRO
        return dc


def plotSetup():
    # labels and grid
    plt.ylabel('frequency, rad/Hz')
    plt.title('Dispersion characteristics')
    plt.grid()
    plt.xlabel('wavenumber, rad/m')
    return


# variables and parameters
var('omega, k1, k2, k,kf, m, a, epsilon1, epsilon2, epsilonf, c,af,H,M,gamma,omegaH,omegaM')
#base parameters for all structures
paramb = [(m, 1), (a, 400 * 10 ** -6), (epsilon1, 1), (epsilon2, 1500), (c, float(299792458))]
pf1 = [(af, 5 * 10 ** -6), (epsilonf, 14), (M, 1750), (H, 1500), (gamma, 2.8 * 10 ** 6)]

# dispersion law
k11 = (sqrt(k ** 2 - omega ** 2 * epsilon1 / c ** 2))
k22 = (sqrt(omega ** 2 * epsilon2 / c ** 2 - k ** 2))
# EMW in dielectric layer
DL = (k2 * a - atan(k1 / k2) - (m - 1) * (pi / 2))  # tan(k2*a)-k1/k2#
pfs = [(omegaH, 2 * pi * gamma * H), (omegaM, 2 * pi * gamma * M)]
mud = 1 + omegaM * omegaH / (omegaH ** 2 - omega ** 2)
mua = omega * omegaM / (omegaH ** 2 - omega ** 2)
mup = (mud ** 2 - mua ** 2) / mud
kf1 = (sqrt(omega ** 2 * epsilonf / c ** 2 * mup - k ** 2))
#ferrite in free space
DLF = 2 * kf * k1 * mud * (mud - mua) * (mud + mua) / (
        kf ** 2 * mud ** 2 + k ** 2 * mua ** 2 - k1 ** 2 * (mud ** 2 - mua ** 2) ** 2) - tan(kf * af)
#ferrite-ferroelectric dilayer
DLFF = 2 * kf * k1 * mup + (k1 ** 2 * mup ** 2 - kf ** 2 - k ** 2 * mua ** 2 / mud ** 2) * tan(af * kf) - tan(
    a * k2) * (kf / k2 * (k2 ** 2 - k1 ** 2) * mup + (
        kf ** 2 * k1 / k2 + (k * k1 / k2 * mua / mud + k2 * mup) * (k * mua / mud + k1 * mup)) * tan(af * kf))

print("dispersion equation:", DL.subs(k1, k11).subs(k2, k22))
k0 = 600  # (2 * pi / (1.55 * 10 ** -6)).evalf()
cc = (c / epsilon1 ** 0.5 * k0).subs(paramb).evalf()
cf = sqrt(omegaH * (omegaM + omegaH)).subs(pfs).subs(pf1).evalf()  # (omegaM/2+omegaH).subs(pfs).subs(pf1).evalf()

print(cf)
# dispersion sub
DDL = DL.subs(k1, k11).subs(k2, k22)
DDLF = DLF.subs(kf, kf1).subs(k1, k11).subs(pfs)
DDLFE = DLFF.subs(kf, kf1).subs(k1, k11).subs(k2, k22).subs(pfs)
# check the root
print("first root in THz:", (Dispersion.fr(DDL.subs(k, k0).subs(paramb),omega, cc) / (2 * pi) * 10 ** -12).evalf())
print("first root in THz:", (Dispersion.fr(DDLF.subs(k, k0).subs(pf1).subs(paramb), omega, cc) / (2 * pi) * 10 ** -12).evalf())
print("first root in THz:", (Dispersion.fr(DDLF.subs(k, k0).subs(pf1).subs(paramb), omega, cf) / (2 * pi) * 10 ** -12).evalf())

# wave number
kk = [float(kk) for kk in np.arange(k0, 500000, 5000)]
print("number of points", len(kk))

# light cone
omega1 = [(c / sqrt(epsilon1) * kk).subs(paramb) for kk in kk]
omega2 = [(c / sqrt(epsilon2) * kk).subs(paramb) for kk in kk]
# calculation dispersion EMW for base dielectric layer
c1 = (c / epsilon1 ** 0.5 * k).subs(k, kk[0]).subs(paramb).evalf()
res0 = Dispersion.calculation(DDL, omega, kk, paramb, c1, 'm1')
#plot results
plt.plot(kk,res0,'g-',kk, omega1, 'c-', kk, omega2, 'b-')
plotSetup()
plt.show()
# calculation of dispersion for dielectric with parameters pv1
pv1 = [(a, 4 * 10 ** -6)]
c2 = (c / epsilon1 ** 0.5 * k).subs(k, kk[0]).subs(pv1).subs(paramb).evalf()
res1 = Dispersion.calculation(DDL, omega, kk, pv1, c2, 'm2')
#plot results
plt.plot(kk, res1, 'g-', kk, omega2, 'b-', kk, omega1, 'c-')
plotSetup()
plt.show()
# calculation of dispersion for dielectric with parameters pv2
pv2 = [(epsilon2, 1100)]
c3 = (c / epsilon1 ** 0.5 * k).subs(k, kk[0]).subs(pv2).subs(paramb).evalf()
res2 = Dispersion.calculation(DDL, omega, kk, pv2, c3, 'm3')
omega3 = [(c / sqrt(epsilon2) * kk).subs(pv2).subs(paramb) for kk in kk]
#plot results
plt.plot(kk, res2, 'g--', kk, omega3, 'b--', kk, omega1, 'c-')
plotSetup()
plt.show()
# calculation of dispersion for ferrite
res3 = Dispersion.calculation(DDLF, omega, kk, pf1, cf, 'm4')
plt.plot(kk, res3, 'r-', kk, omega1, 'c-')
#plot results
plt.axis([int(kk[0]), int(kk[len(kk)-1]), int(res3[0]), int(res3[len(kk)-1])])
plotSetup()
plt.show()
# calculation of dispersion in ferrite-ferroelectric
res4 = Dispersion.calculation(DDLFE, omega, kk, pf1, cf - 0.1 * 10 ** 9, 'm5')
res5 = Dispersion.calculation(DDLFE, omega, kk, pf1, cf + 0.1 * 10 ** 9, 'm6')
#plot results
plt.plot(kk, res3, 'r-', kk, res4, 'b-', kk, res5, 'b-', kk, omega1, 'c-', kk, res2, 'g-')
plt.axis([k0, 6000, 5.5*10**9, 6.7*10**9])
plotSetup()
plt.show()
print("end of calculation")
# save results for ORIGIN PRO
np.savetxt('disp_cone.csv', Dispersion.rearr(kk, omega1), delimiter=" ,")  # save results
np.savetxt('disp_cone2.csv', Dispersion.rearr(kk, omega2), delimiter=" ,")  # save results
np.savetxt('disp_cone3.csv', Dispersion.rearr(kk, omega3), delimiter=" ,")  # save results