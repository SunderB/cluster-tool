import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animate_2d(snapshots, x, y, c, xlim=[-30, 30], ylim=[-30, 30], s=10, cmap="viridis", xlabel="x", ylabel="y", fig_size=[8,8]):    
    # Plot data in animation
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_size)

    # Get first snapshot
    snapshot_mask = (snapshots == 0)

    # Plot the stars
    x_snap = x[snapshot_mask]
    y_snap = y[snapshot_mask]
    c_snap = c[snapshot_mask]
    normalise = mpl.colors.Normalize(vmin=min(c), vmax=max(c))
    scat1 = ax.scatter(x_snap, y_snap, s=s, c=c_snap, cmap=cmap, norm=normalise)

    # Set axis limits and label
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    fig.tight_layout()

    def update(frame):
        # Get current snapshot
        snapshot_mask = (snapshots == frame)

        # Plot the stars
        x_snap = x[snapshot_mask]
        y_snap = y[snapshot_mask]
        c_snap = c[snapshot_mask]
        
        data = np.stack([x_snap, y_snap]).T
        scat1.set_offsets(data)
        scat1.set_array(c_snap)
        scat1.set_cmap(cmap)
        scat1.set_norm(normalise)

        # Return the updated scatter plot
        return scat1

    # Return the animation
    return animation.FuncAnimation(fig=fig, func=update, frames=max(snapshots)+1, interval=1)

def animate_3d(snapshots, x, y, z, c, xlim=[-30, 30], ylim=[-30, 30], zlim=[-30, 30], s=10, cmap="viridis", xlabel="x", ylabel="y", zlabel="z", fig_size=[8,8]) -> animation.FuncAnimation:    
    # Plot data in animation
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_size)

    # Get first snapshot
    snapshot_mask = (snapshots == 0)

    # Plot the stars
    x_snap = x[snapshot_mask]
    y_snap = y[snapshot_mask]
    z_snap = z[snapshot_mask]
    c_snap = c[snapshot_mask]
    normalise = mpl.colors.Normalize(vmin=min(c), vmax=max(c))
    scat1 = ax.scatter(x_snap, y_snap, z=z_snap, s=s, c=c_snap, cmap=cmap, norm=normalise)

    # Set axis limits and label
    ax.set(xlim=xlim, ylim=ylim, zlim=zlim, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
    fig.tight_layout()

    def update(frame):
        # Get current snapshot
        snapshot_mask = (snapshots == frame)

        # Plot the stars
        x_snap = x[snapshot_mask]
        y_snap = y[snapshot_mask]
        z_snap = z[snapshot_mask]
        c_snap = c[snapshot_mask]

        scat1._offsets3d = (x_snap, y_snap, z_snap)
        scat1.set_array(c_snap)
        scat1.set_cmap(cmap)
        scat1.set_norm(normalise)

        # Return the updated scatter plot
        return scat1

    # Return the animation
    return animation.FuncAnimation(fig=fig, func=update, frames=max(snapshots)+1, interval=1)