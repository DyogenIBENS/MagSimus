#!/usr/bin/python
#  -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss

# very useful http://scipy.github.io/devdocs/tutorial/stats.html
# also interesting : http://stackoverflow.com/questions/10678546/creating-new-distributions-in-scipy
# http://stackoverflow.com/questions/22854332/using-the-methods-of-scipys-rv-continuous-when-creating-a-cutom-continuous-dist

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TAKE CARE, since pdf are discrete functions (pmf), the X array must only contain integer values !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import sys

nbGenesInChrom = 1330
minInvLength = 1
maxInvLength = nbGenesInChrom - 1

XmarginWidth = 100
X = range(0-XmarginWidth, maxInvLength+1 + XmarginWidth)
ylim_pmf = 5.0/maxInvLength

# not a discrete distribution
# uniformDistrib = ss.uniform(loc=minInvLength, scale=maxInvLength)

figure = plt.figure()
axes = []
for i in range(1,6):
    axes.append(figure.add_subplot('32' + str(i)))
for ax in axes:
    ax.set_xlim((min(X),max(X))) # this fixes the x-axis for all plots
    ax.set_ylim([0,1])
    ax.set_xlabel('longueur de l\'inversion (genes)')
    ax.set_ylabel('pmf', color='b')
    for tl in ax.get_yticklabels():
        tl.set_color('b')
    ax.set_ylim(0, ylim_pmf)
    ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='y')
axesTwin = [ax.twinx() for ax in axes]
for ax in axesTwin:
    ax.set_xlim((min(X),max(X))) # this fixes the x-axis for all plots
    ax.set_ylim([0,1]) # this fixes the x-axis for all plots
    ax.set_ylabel('cdf', color='k')
    for tl in ax.get_yticklabels():
        tl.set_color('k')

# discrete uniform
idxAxis = 0
# sub-classing inspired from http://stackoverflow.com/questions/23276631/scipy-generating-custom-random-variable-from-pmf
class uniform_discrete(ss.rv_discrete):
    def _pmf(self, X, nbGenes):
        # listOfvaluesAboveBounds = (X < self.a)  | (self.b < X) # where values are not in the support
        # X[listOfvaluesAboveBounds] = 0 # set all the values not in the support to 0
        L = nbGenes
        assert np.all(L > 1), 'not enough gene for an inversion'
        return np.true_divide(np.ones(len(X)), L-1)
#freeze the distribution (automatically associate the arg: nbGenes=nbGenesInChrom
uniform_discreteDistrib = uniform_discrete(a=minInvLength, b=maxInvLength).__call__(nbGenesInChrom)
# Le support de la distribution = [minInvLength, maxInvLength], les bornes sont inclues
uniform_discreteDistrib.name = 'uniforme %s'
Ypmf = uniform_discreteDistrib.pmf(X)
axes[idxAxis].step(X, Ypmf, color='b', label=(uniform_discreteDistrib.name % 'pmf'))
axes[idxAxis].fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
axes[idxAxis].set_title(uniform_discreteDistrib.name % '')
axesTwin[idxAxis].step(X, uniform_discreteDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(uniform_discreteDistrib.name % 'cdf'))

# Triangular
idxAxis = 1
class rb_triangular_distribution(ss.rv_discrete):
    def _pmf(self, X, nbGenes):
        #listOfvaluesAboveBounds = (X < self.a)  | (X > self.b) # where values are not in the support
        L = nbGenes
        res = np.true_divide(2 * (L-X+1), np.multiply(L, L+1) - 2)
        # if len(listOfvaluesAboveBounds) > 0:
        #     res[listOfvaluesAboveBounds] = 0
        return res
#freeze the distribution (automatically associate the arg: nbGenes=nbGenesInChrom
rbTriangularDistrib = rb_triangular_distribution(a=minInvLength, b=maxInvLength).__call__(nbGenesInChrom)
rbTriangularDistrib.name = 'RB triangulaire %s'
# plt.figure()
Ypmf = rbTriangularDistrib.pmf(X)
axes[idxAxis].step(X, Ypmf, color='b', label=(rbTriangularDistrib.name % 'pmf'))
axes[idxAxis].fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
axes[idxAxis].set_title(rbTriangularDistrib.name % '')
axesTwin[idxAxis].step(X, rbTriangularDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(rbTriangularDistrib.name % 'cdf'))
# plt.show()

# 1/x
idxAxis = 2
class oneOverX(ss.rv_discrete):
    def _pmf(self, X, wholleIntegralOneOverX):
        #listOfvaluesAboveBounds = (X < self.a)  | (X > self.b) # where values are not in the support
        #X[listOfvaluesAboveBounds] = 1.0 # set all the values not in the support to 1
        res = np.true_divide(1.0, np.multiply(X, wholleIntegralOneOverX))
        # res = np.true_divide(res, wholleIntegralOneOverX)
        # if len(listOfvaluesAboveBounds) > 0:
        #     res[listOfvaluesAboveBounds] = 0
        return res
#binWidth = 1
#wholleIntegralOneOverX = np.sum(np.true_divide(1, X))
oneOverXDistrib = oneOverX(a=minInvLength, b=maxInvLength).__call__(1.0)
wholleIntegralOneOverX = np.sum(oneOverXDistrib.pmf(X))
# Take care to not use a X that passes through zero, otherwise negative values for the pmf
oneOverXDistrib = oneOverX(a=minInvLength, b=maxInvLength).__call__(wholleIntegralOneOverX)#*1.122)
oneOverXDistrib.name = '1/x %s'
# plt.figure()
Ypmf = oneOverXDistrib.pmf(X)
axes[idxAxis].step(X, Ypmf, color='b', label=(oneOverXDistrib.name % 'pmf'))
axes[idxAxis].fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
axes[idxAxis].set_title(oneOverXDistrib.name % 'distribution')
axesTwin[idxAxis].step(X, oneOverXDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(oneOverXDistrib.name % 'cdf'))
# plt.show()

# gamma
idxAxis = 3
shape = 0.1
scale = 800
# freeze distribution
#gammaDistribClassic = ss.gamma(0.2, scale=300.0)
# gamma truncated
# TODO be sure that the integration is equal to one
class gammaTruncatedDiscrete(ss.rv_discrete):
    # def _pmf(self, X, wholleIntegralOneOverX, shape, scale):
    #     res = np.true_divide(ss.gamma(shape, scale=scale).pdf(X), wholleIntegralOneOverX)
    #     res[X<self.a] = 0
    #     res[X>self.b] = 0
    #     return res
    def _cdf(self, X, shape, scale):
        res = ss.gamma(shape, scale=scale).cdf(X)
        # res[X<self.a] = 0
        # res[X>self.b] = 0
        return res
#wholleIntegralOneOverX = np.sum(gammaDistribClassic.pdf(X))
gammaTruncatedDistrib = gammaTruncatedDiscrete(a=minInvLength, b=maxInvLength).__call__(shape, scale)
print >> sys.stderr, "%.2f%% micro-rearrangements (1<=x<=1 anc. genes) with Gamma(shape=%s,scale=%s)" % \
                     ((gammaTruncatedDistrib.cdf(1) * 100), shape, scale)
print >> sys.stderr, "%.2f%% micro-rearrangements (1<=x<=2 anc. genes) with Gamma(shape=%s,scale=%s)" % \
                     ((gammaTruncatedDistrib.cdf(2) * 100), shape, scale)
print >> sys.stderr, "%.2f%% micro-rearrangements (1<=x<=3 anc. genes) with Gamma(shape=%s,scale=%s)" % \
                     ((gammaTruncatedDistrib.cdf(3) * 100), shape, scale)
print >> sys.stderr, "%.2f%% micro-rearrangements (1<=x<=4 anc. genes) with Gamma(shape=%s,scale=%s)" % \
                     ((gammaTruncatedDistrib.cdf(4) * 100), shape, scale)
print >> sys.stderr, "%.2f%% micro-rearrangements (1<=x<=5 anc. genes) with Gamma(shape=%s,scale=%s)" % \
                     ((gammaTruncatedDistrib.cdf(5) * 100), shape, scale)

# wholleIntegralOneOverX = np.sum(gammaTruncatedDistrib.pmf(X))
# gammaTruncatedDistrib = gammaTruncatedDiscrete(a=minInvLength, b=maxInvLength).__call__(shape, scale)
gammaTruncatedDistrib.name = 'Gamma(shape=%s, scale=%s) %s' % (shape, scale, '%s')
# plt.figure()
Ypmf = gammaTruncatedDistrib.pmf(X)
axes[idxAxis].step(X, Ypmf, color='b', label=(gammaTruncatedDistrib.name % 'pmf'))
axes[idxAxis].fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
axes[idxAxis].set_title(gammaTruncatedDistrib.name % '')
axesTwin[idxAxis].step(X, gammaTruncatedDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(gammaTruncatedDistrib.name % 'cdf'))
# plt.show()
#

# intermediary von Mises
mu=0.1
kappa = 2
# freeze distribution
X_vmi = np.linspace(-1 + 2*mu, 1, num=1000)
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.vonmises_line.html#scipy.stats.vonmises_line
loc = mu
scale = (1.0 - mu) / np.pi
# create a von Mises with a support [-1+2mu, 1]
vmi = ss.vonmises_line(kappa, loc=loc, scale=scale)
vmi.name = 'vonMises (mu=%s, kappa=%s) sur [%s, %s] %s' % (mu, kappa, -1 + 2*mu, 1, '%s')
# plt.figure()
# Y_pdf= vmi.pdf(X_vmi)
# plt.vlines(mu, 0, max(Y_pdf))
# plt.step(X_vmi, Y_pdf, color='b', label=(vmi.name % 'pmf'))
# plt.fill_between(X_vmi, Y_pdf, 0, color=(0.75,0.75,0.75), step='pre')
# plt.xlabel('x-axis')
# plt.ylabel('pdf', color='b')
# for tl in plt.gca().get_yticklabels():
#     tl.set_color('b')
# plt.twinx().step(X_vmi, vmi.cdf(X_vmi), color='k', linestyle='-', linewidth=2, label=(vmi.name % 'cdf'))
# plt.ylabel('cdf', color='k')
# plt.title(vmi.name % '')
# plt.show()

# modified von Mises
idxAxis = 4
mu=0.01
kappa = 2
# # create a von Mises with a support [-1+2mu, 1]
# loc = mu
# scale = (1.0 - mu) / np.pi
# vmi = ss.vonmises_line(kappa, loc=loc, scale=scale)
# The [0,mu] interval must at least contain 20 bins
# base change [1, maxInv] -> [1/maxInv,1]
def b0(x): return np.true_divide(x, maxInvLength)
# # base change [0,1] -> [-mu/(1-mu)pi, pi]
# def b1(x): return np.true_divide(x-mu, 1-mu) * np.pi
# s1 = np.pi/float(1-mu)
# # base change [0,mu] -> [-mu/(1-mu)pi, -pi]
# def b2(x): return b1(x) - np.true_divide(x, mu)* np.pi
# s2 = float(1-2*mu)/((1-mu)*mu) * np.pi
# TODO
# base change [0,1] -> [-1+2mu, 1]
#def b1(x): return 2 * (1 - mu) * x + 2 * mu - 1
# base change [0,mu] -> [0, -1+2mu]
def b2(x): return x * float(2 * mu - 1) / mu
# TODO be sure that the integration is equal to one



class modifiedVonMises(ss.rv_discrete):
    # def _pmf(self, X, wholleIntegralOneOverX, mu, kappa):
    #     # http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.piecewise.html
    #     X0 = b0(X)
    #     X1 = b1(X0)
    #     # probability mass 1
    #     pm1 = np.multiply(ss.vonmises(kappa, loc=0).pdf(X1), s1)
    #     # print 'pm1 =', pm1
    #     X2 = b2(X0)
    #     # probability mass 2
    #     pm2 = np.multiply(ss.vonmises(kappa, loc=0).pdf(X2), s2)
    #     # print 'pm2 =', pm2
    #     # print 'pm1 + pm2 =', pm1 + pm2
    #     res = np.where(X0>mu, pm1, pm1 + pm2)
    #     res = np.true_divide(res, maxInvLength)
    #     res = np.true_divide(res, wholleIntegralOneOverX)
    #     res[X < self.a] = 0
    #     res[X > self.b] = 0
    #     # Waring !! if some values are over 1 they are -> to 1
    #     return res
    def _cdf(self, X, mu, kappa):
        # create a von Mises with a support [-1+2mu, 1]
        loc = mu
        scale = (1.0 - mu) / np.pi
        vmi = ss.vonmises_line(kappa, loc=loc, scale=scale)

        X0 = b0(X)
        # cumulated probability 1 :  mu <= X0
        cp_over0 = vmi.cdf(X0)
        cp_between0AndMu = cp_over0 - vmi.cdf(0)
        # 0 <= X0 < mu: X2 = b2(X0), sinon 0
        X2 = np.where(((X0<mu) & (X0>=0)), b2(X0), np.zeros(len(X0)))
        # -1 + 2 * mu <= X2 < 0: cp3 = cdf(X2), sinon 0
        cp_i_betweenMinus1PlusTwoMuAnd0 = np.where(((X2>=(-1+2 * mu)) & (X2<0)), vmi.cdf(X2), np.zeros(len(X2)))
        # 0 <= X0 < mu: cp2 = cdf(0) - cp3, sinon 0
        # cumulated probability 2 : 0<= X0 < mu
        cp_additionalBetween0AndMu =  np.where(((X0<mu) & (X0>=0)), vmi.cdf(0) - cp_i_betweenMinus1PlusTwoMuAnd0 , np.zeros(len(X2)))
        # mu <= X0, res = cp1, sinon, cp1 + cp2
        res = np.where(X0>=mu, cp_over0, cp_between0AndMu + cp_additionalBetween0AndMu)
        #res[(X < self.a) | (X0 > self.b)] = 0
        return res

modifiedVonMisesDistrib = modifiedVonMises(a=minInvLength, b=maxInvLength).__call__(mu, kappa)
#print modifiedVonMisesDistrib.pmf(0)
# C'est là que ça cloche, ce résultat est beaucoup trop élevé
# wholleIntegralOneOverX = np.sum(modifiedVonMisesDistrib.pmf(X))
# modifiedVonMisesDistrib = modifiedVonMises(a=minInvLength, b=maxInvLength).__call__(wholleIntegralOneOverX, mu, kappa)
modifiedVonMisesDistrib.name = 'vonMises_modifiee(mu=%s, kappa=%s) %s' % (mu, kappa, '%s')
Ypmf = modifiedVonMisesDistrib.pmf(X)
axes[idxAxis].step(X, Ypmf, color='b', label=(modifiedVonMisesDistrib.name % 'pmf'))
axes[idxAxis].fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
axes[idxAxis].set_title(modifiedVonMisesDistrib.name % '')
#axes[idxAxis].vlines(mu*maxInvLength, 0, +1)
axesTwin[idxAxis].step(X, modifiedVonMisesDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(modifiedVonMisesDistrib.name % 'cdf'))
plt.show(block=False)

fig, axes = plt.subplots(nrows=2, ncols=2)
print axes
ax = axes[0][0]
mu=0.001
kappa=2
modifiedVonMisesDistrib = modifiedVonMises(a=minInvLength, b=maxInvLength).__call__(mu, kappa)
# wholleIntegralOneOverX = np.sum(modifiedVonMisesDistrib.pmf(X))
# modifiedVonMisesDistrib = modifiedVonMises(a=minInvLength, b=maxInvLength).__call__(mu, kappa)
modifiedVonMisesDistrib.name = 'vonMises_modifiee(mu=%s, kappa=%s) %s' % (mu, kappa, '%s')
Ypmf = modifiedVonMisesDistrib.pmf(X)
ax.step(X, Ypmf, color='b', label=(modifiedVonMisesDistrib.name % 'pmf'))
ax.set_ylim([0,ylim_pmf])
ax.fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
ax.set_title(modifiedVonMisesDistrib.name % 'distribution')
ax.twinx().step(X, modifiedVonMisesDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(modifiedVonMisesDistrib.name % 'cdf'))

ax = axes[0][1]
mu=0.001
kappa=3
modifiedVonMisesDistrib = modifiedVonMises(a=minInvLength, b=maxInvLength).__call__(mu, kappa)
# wholleIntegralOneOverX = np.sum(modifiedVonMisesDistrib.pmf(X))
# modifiedVonMisesDistrib = modifiedVonMises(a=minInvLength, b=maxInvLength).__call__(mu, kappa)
modifiedVonMisesDistrib.name = 'vonMises_modifiee(mu=%s, kappa=%s) %s' % (mu, kappa, '%s')
Ypmf = modifiedVonMisesDistrib.pmf(X)
ax.step(X, Ypmf, color='b', label=(modifiedVonMisesDistrib.name % 'pmf'))
ax.set_ylim([0,ylim_pmf])
ax.fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
ax.set_title(modifiedVonMisesDistrib.name % 'distribution')
ax.twinx().step(X, modifiedVonMisesDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(modifiedVonMisesDistrib.name % 'cdf'))


#GAMMA
ax = axes[1][0]
# shape = 0.2
# scale = 300
shape = 0.05
scale = 5000
gammaTruncatedDistrib = gammaTruncatedDiscrete(a=minInvLength, b=maxInvLength).__call__(shape, scale)
# wholleIntegralOneOverX = np.sum(gammaTruncatedDistrib.pmf(X))
# gammaTruncatedDistrib = gammaTruncatedDiscrete(a=minInvLength, b=maxInvLength).__call__(wholleIntegralOneOverX, shape, scale)
gammaTruncatedDistrib.name = 'Gamma(shape=%s, scale=%s) %s' % (shape, scale, '%s')
Ypmf = gammaTruncatedDistrib.pmf(X)
ax.step(X, Ypmf, color='b', label=(gammaTruncatedDistrib.name % 'pmf'))
ax.set_ylim([0,ylim_pmf])
ax.fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
ax.set_title(gammaTruncatedDistrib.name % 'distribution')
ax.twinx().step(X, gammaTruncatedDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(gammaTruncatedDistrib.name % 'cdf'))

ax = axes[1][1]
shape = 0.1
scale = 800
gammaTruncatedDistrib = gammaTruncatedDiscrete(a=minInvLength, b=maxInvLength).__call__(shape, scale)
# wholleIntegralOneOverX = np.sum(gammaTruncatedDistrib.pmf(X))
# gammaTruncatedDistrib = gammaTruncatedDiscrete(a=minInvLength, b=maxInvLength).__call__(wholleIntegralOneOverX, shape, scale)
gammaTruncatedDistrib.name = 'Gamma(shape=%s, scale=%s) %s' % (shape, scale, '%s')
Ypmf = gammaTruncatedDistrib.pmf(X)
ax.step(X, Ypmf, color='b', label=(gammaTruncatedDistrib.name % 'pmf'))
ax.set_ylim([0,ylim_pmf])
ax.fill_between(X, Ypmf, 0, color=(0.75,0.75,0.75), step='pre')
ax.set_title(gammaTruncatedDistrib.name % 'distribution')
ax.twinx().step(X, gammaTruncatedDistrib.cdf(X), color='k', linestyle='-', linewidth=2, label=(gammaTruncatedDistrib.name % 'cdf'))

plt.show()











# # deprecated
#
# # modified von Mises (mu=0.1, kappa=3), between [0,1] donne un beau résultat
# mu=0.1
# kappa = 3
# vmc = ss.vonmises(kappa, loc=0)
# # The [0,mu] interval must at least contain 20 bins
# numDivs = 20 * (1.0/mu)
# print numDivs
# X = np.linspace(0,1,num=numDivs)
# b1 = lambda x: np.true_divide(x-mu, 1-mu) * np.pi
# s1 = np.pi/float(1-mu)
# b2 = lambda x: b1(x) - np.true_divide(x, mu)* np.pi
# s2 = float(1-2*mu)/((1-mu)*mu) * np.pi
# class gammaTruncatedDiscrete(ss.rv_continuous):
#     def _pdf(self, X):
#         listOfvaluesAboveBounds = (X < self.a)  | (self.b < X)
#         indicesOverMu = mu <= X
#         indicesBelowMu = X < mu
#         # http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.piecewise.html
#         res = np.piecewise(X, [indicesOverMu, indicesBelowMu],
#                            [lambda x: vmc.pdf(b1(x))*s1, lambda x: vmc.pdf(b1(x))*s1 + vmc.pdf(b2(x))*s2])
#         res[listOfvaluesAboveBounds] = 0
#         return res
# #wholleIntegralOneOverX = np.sum(gammaDistribClassic.pdf(X))
# gammaTruncatedDistrib = gammaTruncatedDiscrete(a=0, b=1)#.__call__(wholleIntegralOneOverX)
# gammaTruncatedDistrib.name = 'vonMises_modifiee(mu=%s, kappa=%s) over [%s, %s] %s' % (mu, kappa, 0, 1, '%s')
# plt.figure()
# plt.step(X, gammaTruncatedDistrib.pdf(X), color='k', label=(gammaTruncatedDistrib.name % 'pdf'))
# plt.step(X, gammaTruncatedDistrib.cdf(X), color='k', label=(gammaTruncatedDistrib.name % 'cdf'))
# plt.legend()
# plt.show()

