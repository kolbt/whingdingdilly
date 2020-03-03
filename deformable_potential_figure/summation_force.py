'''
Let's compute the average force from a reference on a neighbor
'''

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

thetaL = np.arange(0., 2*np.pi, 0.01)
thetaR = np.arange(0., 2*np.pi, 0.01)

## Let's write a 3-body script
#fs = []
#img = []
#for i in thetaL:
#    img.append([])
#    x = len(img) - 1
#    for j in thetaR:
#        # Left particle points left
#        if (np.pi / 2.) < i < (3. * np.pi / 2.):
#            # Skip, if ref also points left
#            if (np.pi / 2.) < j < (3. * np.pi / 2.):
#                fs.append(0.)
#                img[x].append(fs[-1])
#                continue
#            if np.cos(j) >= 0:
#                fs.append(np.cos(j))
#                img[x].append(fs[-1])
#                continue
#        # Left particle has rightward component
#        else:
#            # If the sum points left
#            if np.cos(i) + np.cos(j) <= 0:
#                fs.append(0.)
#                img[x].append(fs[-1])
#                continue
#            # If the sum is positive
#            else:
#                fs.append(np.cos(i) + np.cos(j))
#                img[x].append(fs[-1])
#                continue
  
# Let's write a 2-body script
fs = []
img = []
for i in thetaL:
    img.append([])
    x = len(img) - 1
    for j in thetaR:
        # Does the left particle point right
        if i <= (np.pi / 2.) or i >= (3 * np.pi / 2.):
            # Right particle points left
            if (np.pi / 2.) <= j <= (3 * np.pi / 2.):
                fs.append(max(np.cos(i), -np.cos(j)))
                img[x].append(fs[-1])
                continue
            # Right particle points right
            else:
                fs.append(np.cos(i))
                img[x].append(fs[-1])
                continue
        # Ignore left particle's orientation
        else:
            # Right particle points left
            if (np.pi / 2.) <= j <= (3 * np.pi / 2.):
                fs.append(-np.cos(j))
                img[x].append(fs[-1])
                continue
            # They point away (in their x-component) from each other
            else:
                fs.append(0.)
                img[x].append(fs[-1])
                continue
                
avgF = sum(fs) / len(fs)
print(avgF)

plt.imshow(img, extent=[0,2*np.pi,0,2*np.pi])
plt.colorbar()
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)
ax = plt.gca()
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
plt.show()
            
