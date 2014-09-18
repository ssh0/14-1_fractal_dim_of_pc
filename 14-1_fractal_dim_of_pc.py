#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, June 2014.

from Tkinter import *
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as optimize


class Percolation:

    def __init__(self):
        self.sub = None
        self.L = 61  # lattice size
        self.p = 0.5
        self.lattice = np.zeros([self.L, self.L], dtype=bool)

    def percolate(self, p=0.5):
        self.p = p
        if self.sub is None or not self.sub.winfo_exists():
            lattice = self.lattice
            rn = np.random.random([self.L, self.L])
            lattice[rn < p] = True
            lattice[rn >= p] = False
            self.lattice = lattice

    def labeling(self):
        label = np.zeros([self.L + 2, self.L + 2], dtype=int)
        n = 1
        r = range(1, self.L + 1)
        for i in r:
            for j in r:
                if self.lattice[i - 1][j - 1]:
                    nn = []
                    if label[i - 1][j] > 0:
                        nn.append(label[i - 1][j])
                    if label[i][j - 1] > 0:
                        nn.append(label[i][j - 1])
                    if len(nn) > 0:
                        label[i][j] = min(nn)
                    else:
                        label[i][j] = n
                        n += 1
        tag = range(1, n + 1)

        for i in reversed(r):
            for j in reversed(r):
                if label[i][j] > 0:
                    nn = []
                    if label[i + 1][j] > 0:
                        nn.append(label[i + 1][j])
                    if label[i][j + 1] > 0:
                        nn.append(label[i][j + 1])
                    nn.append(label[i][j])
                    min_tag = min(nn)
                    nn = set([x for x in nn if x != min_tag])
                    for t in nn:
                        tag[t - 1] = min_tag
                        label[label == t] = tag[t - 1]

        self.lattice = label[1:-1, 1:-1]
        left = set(self.lattice[0])
        right = set(self.lattice[self.L - 1])
        top = set([self.lattice[t][0] for t in range(self.L)])
        bottom = set([self.lattice[t][self.L - 1] for t in range(self.L)])
        self.ptag = (left.intersection(right)
                     | top.intersection(bottom)) - set([0])

    def view_expansion(self):
        lattice = np.zeros([self.L, self.L])
        lattice[self.lattice == list(self.ptag)[0]] = 1
        M_b = []
        s = np.sum
        ave = np.average
        append = M_b.append
        for k in range(1, int(self.L) / 2):
            nonzero = np.nonzero(lattice[k:-k, k:-k])
            tmp = np.array([0])
            for i, j in zip(nonzero[0] + k, nonzero[1] + k):
                tmp = np.append(
                    tmp, s(lattice[i - k:i + k + 1, j - k:j + k + 1]))
            append(ave(tmp))

        b = np.array([2. * k + 1 for k in range(1, int(self.L) / 2)])
        M_b = np.array(M_b)

        def fit_func(parameter0, b, M_b):
            log = np.log
            c1 = parameter0[0]
            c2 = parameter0[1]
            residual = log(M_b) - c1 - c2 * log(b)
            return residual

        parameter0 = [0.1, 2.0]
        result = optimize.leastsq(
            fit_func, parameter0, args=(b[:-1], M_b[:-1]))
        c1 = result[0][0]
        D = result[0][1]

        def fitted(b, c1, D):
            return np.exp(c1) * (b ** D)

        fig = plt.figure("Fractal Dimesion")
        ax = fig.add_subplot(111)
        ax.plot(b, M_b, '-o', label="p = %f" % self.p)
        ax.plot(b, fitted(b, c1, D), label="fit func: D = %f" % D)
        ax.set_xlabel(r'$b$', fontsize=16)
        ax.set_ylabel(r'$M(b)$', fontsize=16)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ymargin(0.05)
        fig.tight_layout()
        plt.legend(loc='best')
        print "D = %f" % D
        plt.show()

    def draw_canvas(self):
        default_size = 640  # default size of canvas
        r = int(default_size / (2 * self.L))
        fig_size = 2 * r * self.L
        margin = 10
        sub = Toplevel()

        sub.title('figure  ' + '(p=%s)' % str(self.p))
        self.canvas = Canvas(sub, width=fig_size + 2 * margin,
                             height=fig_size + 2 * margin)
        self.canvas.create_rectangle(margin, margin,
                                     fig_size + margin, fig_size + margin,
                                     outline='black', fill='white')
        self.canvas.pack()

        c = self.canvas.create_rectangle
        rect = self.lattice
        colors = ['blue', 'green', 'red', 'purple']
        colordict = dict(zip(list(self.ptag),
                             colors * (int(len(self.ptag) / len(colors)) + 1)
                             )
                         )

        nonzero_rect = np.nonzero(rect)
        for m, n in zip(nonzero_rect[0], nonzero_rect[1]):
            if rect[m][n] in self.ptag:
                c(2 * m * r + margin, 2 * n * r + margin,
                  2 * (m + 1) * r + margin, 2 * (n + 1) * r + margin,
                  outline=colordict[rect[m][n]], fill=colordict[rect[m][n]])
            else:
                c(2 * m * r + margin, 2 * n * r + margin,
                  2 * (m + 1) * r + margin, 2 * (n + 1) * r + margin,
                  outline='black', fill='black')


class TopWindow:

    def quit(self):
        self.root.destroy()
        sys.exit()

    def show_window(self, pr, pushed, b4_pushed):
        self.root = Tk()
        self.root.title('Percolation')
        f = Frame(self.root)
        self.label = Label(f, text='p =')
        self.label.pack(side='left')
        self.entry = Entry(f, width=20)
        self.entry.pack(side='left')
        self.entry.delete(0, END)
        self.entry.insert(0, 0.5927)
        self.entry.focus_set()

        b1 = Button(f, text='run', command=pushed)
        b1.pack(side='left', expand=YES, fill='x')

        b4 = Button(f, text='count', command=b4_pushed)
        b4.pack(side='left', expand=YES, fill='x')

        b2 = Button(f, text='write canvas to sample.eps', command=pr)
        b2.pack(side='left', expand=YES, fill='x')

        b3 = Button(f, text='quit', command=self.quit)
        b3.pack(side='right', expand=YES, fill='x')

        f.pack(fill='x')

        self.root.mainloop()

if __name__ == '__main__':
    top = TopWindow()
    per = Percolation()
    count = 1

    def pr():
        global count
        p = float(top.entry.get())
        d = per.canvas.postscript(file="figure_%d(p=%s).eps" % (count, str(p)))
        print "saved the figure to a eps file"
        count += 1

    def pushed():
        p = float(top.entry.get())
        per.percolate(p)
        per.labeling()
        per.draw_canvas()

    def b4_pushed():
        if not per.ptag:
            print "Can't count. (There is no percolation cluster.)"
        else:
            per.view_expansion()

    top.show_window(pr, pushed, b4_pushed)
