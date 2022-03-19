import matplotlib.pyplot as plt


def scatter(xs, ys, title = None, sizes = None, show = True, xlim = None, ylim = None, cols = None):
    fig, ax = plt.subplots()
    if xlim:
        ax.set_xlim(xlim) # [xmin, xmax]
    if ylim:
        ax.set_ylim(ylim) # [ymin, ymax]
    ax.scatter(xs, ys, s = sizes, c = cols)
    if title:
        plt.title(title)
    if show:
        plt.show()
    return fig
