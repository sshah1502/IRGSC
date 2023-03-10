import os
import sys
import numpy as np
from matplotlib import pyplot as plt

p1 = np.genfromtxt('tf1_oid_stack_sshah1502.csv', delimiter=',', skip_header=1)
pis=p1[:,13]
indfis = np.where(pis!=-999)[0]

p0 = np.genfromtxt('tf1_oid_sshah1502.csv', delimiter=',', skip_header=1)
pim=p0[:,13]
indfim = np.where(pim!=-999)[0]
print(pis[indfis], pim[indfim])

bins=50
bins = np.histogram(np.hstack((np.int64(pis[indfis]), np.int64(pim[indfim]))), bins=bins)[1]

plt.clf()
plt.figure(figsize=(10,10))
y0,x0,_=plt.hist(pis[indfis], bins, histtype="step", edgecolor='black', label='PS1 i from Stack photometry')
y1,x1,_=plt.hist(pim[indfim], bins, histtype="step", edgecolor='red', label='PS1 i from mean photometry')
plt.xlabel('$i_{psf}$', fontsize=18)
plt.ylabel('Counts', fontsize=18)
plt.legend(loc='best', fontsize=18)
plt.grid(linestyle='--')
percent_completeness=90
completeness_mags = np.percentile(pis[indfis],90)
completeness_magm = np.percentile(pim[indfim],90)
#completeness_mag = x1[np.where(y1==np.max(y1))[0]]
#ercent_completeness = 100 * np.max(y1)/y0[np.where(x0==completeness_mag  )[0]]
plt.title('The catalog is %0.2f' % (percent_completeness) + '%% complete in $i_{psf, stack}$ at %0.2f'%(completeness_mags) \
+' '+ 'and' + ' ' + '$i_{psf, mean}$ at  %0.2f' % (completeness_magm))
plt.savefig('hist_i.png')
plt.clf()

sys.exit(0)

bins=50
bins = np.histogram(np.hstack((np.int64(cj), np.int64(app_j))), bins=bins)[1]

plt.clf()
plt.figure(figsize=(10,10))
y0,x0,_=plt.hist(app_j, bins, facecolor='g', edgecolor='k', alpha=0.3, label='Besancon J')
y1,x1,_=plt.hist(cj, bins, facecolor='b', edgecolor='k', alpha=0.3, label='Computed J')
plt.xlabel('J magnitude', fontsize=18)
plt.ylabel('Counts', fontsize=18)
plt.legend(loc='best', fontsize=15)
plt.grid(linestyle='--')
percent_completeness=90
completeness_mag = np.percentile(cj,90)
#completeness_mag = x1[np.where(y1==np.max(y1))[0]]
#ercent_completeness = 100 * np.max(y1)/y0[np.where(x0==completeness_mag  )[0]]
plt.title('The catalog is %0.2f' % (percent_completeness) + '%% complete in $J_{computed}$ is at %0.2f'%(completeness_mag))
plt.savefig('hist_bj.png')
plt.clf()
