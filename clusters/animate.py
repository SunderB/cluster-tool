import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animate_2d(snapshots, x, y, c, xlim=[-40, 40], ylim=[-40, 40], s=10, highlight_mask=None, cmap="viridis", xlabel="x", ylabel="y", fig_size=[9,8]):    
    # Plot data in animation
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_size)

    using_highlight_mask = str(type(highlight_mask)) != "<class 'NoneType'>"

    # Get first snapshot
    snapshot_mask = (snapshots == 0)
    if (using_highlight_mask):
        # Exclude highlighted stars from main plot
        snapshot_mask = np.logical_and(snapshot_mask, np.logical_not(highlight_mask))
        snapshot_highlight_mask = np.logical_and(snapshot_mask, highlight_mask)

    # Plot the stars
    x_snap = x[snapshot_mask]
    y_snap = y[snapshot_mask]
    c_snap = c[snapshot_mask]
    normalise = mpl.colors.Normalize(vmin=min(c), vmax=max(c))
    scat1 = ax.scatter(x_snap, y_snap, s=s, c=c_snap, cmap=cmap, norm=normalise)

    # Plot highlighted stars
    if (using_highlight_mask):
        x_snap = x[snapshot_highlight_mask]
        y_snap = y[snapshot_highlight_mask]
        c_snap = c[snapshot_highlight_mask]
        scat2 = ax.scatter(x_snap, y_snap, s=s*5, c=c_snap, marker="^", cmap=cmap, norm=normalise)

    # Set axis limits and label
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)
    fig.tight_layout()

    def update(frame):
        # Get current snapshot
        snapshot_mask = (snapshots == frame)
        if (using_highlight_mask):
            # Exclude highlighted stars from main plot
            snapshot_mask = np.logical_and(snapshot_mask, np.logical_not(highlight_mask))
            snapshot_highlight_mask = np.logical_and(snapshot_mask, highlight_mask)

        # Plot the stars
        x_snap = x[snapshot_mask]
        y_snap = y[snapshot_mask]
        c_snap = c[snapshot_mask]
        
        data = np.stack([x_snap, y_snap]).T
        scat1.set_offsets(data)
        scat1.set_array(c_snap)
        scat1.set_cmap(cmap)
        scat1.set_norm(normalise)
        scat1.set_sizes(np.ones(len(c_snap))*s)

         # Plot highlighted stars
        if (using_highlight_mask):
            x_snap = x[snapshot_highlight_mask]
            y_snap = y[snapshot_highlight_mask]
            c_snap = c[snapshot_highlight_mask]
            data = np.stack([x_snap, y_snap]).T
            scat2.set_offsets(data)
            scat2.set_array(c_snap)
            scat2.set_cmap(cmap)
            scat2.set_norm(normalise)
            scat2.set_sizes(np.ones(len(c_snap))*s*5)

            # Return the updated scatter plots
            return scat1, scat2

        # Return the updated scatter plot
        return scat1

    # Return the animation
    return animation.FuncAnimation(fig=fig, func=update, frames=max(snapshots)+1, interval=1)

def animate_3d(snapshots, x, y, z, c, xlim=[-40, 40], ylim=[-40, 40], zlim=[-40, 40], s=10, cmap="viridis", xlabel="x", ylabel="y", zlabel="z", fig_size=[9,8]) -> animation.FuncAnimation:    
    # Plot data in animation
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    fig.set_size_inches(tuple(fig_size))

    # Get first snapshot
    snapshot_mask = (snapshots == 0)

    # Plot the stars
    x_snap = x[snapshot_mask]
    y_snap = y[snapshot_mask]
    z_snap = z[snapshot_mask]
    c_snap = c[snapshot_mask]
    normalise = mpl.colors.Normalize(vmin=min(c), vmax=max(c))
    scat1 = ax.scatter(x_snap, y_snap, z_snap, s=s, c=c_snap, cmap=cmap, norm=normalise)

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