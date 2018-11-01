import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def add_arrow3D(X, Y, Z, ax=None, arrowstyle="-|>", lw=1, **kwargs):
    if ax is None:
        ax = plt.gca(projection='3d')
    arrow = Arrow3D(X, Y, Z, mutation_scale=10, lw=lw, arrowstyle=arrowstyle, **kwargs)
    ax.add_artist(arrow)


if __name__ == '__main__':
    add_arrow3D([0, 1], [0, 1], [0, 1])