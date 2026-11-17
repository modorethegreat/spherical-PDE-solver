from mode_benchmark import numerical_step
from sph_harm_transform import ISHT, SHT, ls, phi, phi_mesh, theta_mesh
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import animation

# parameters
H = 2
R = 1
g = 1
b = 0.5
dt = 0.001


def simulation_step(h_lm, v_lm):
    coupling = (H / R) * (ls * (ls + 1))[:, None]
    return numerical_step(h_lm, v_lm, coupling, g / R, b, dt)


def simulate():
    h_lm = SHT(np.full(theta_mesh.shape, H, dtype=float))
    v_lm = np.zeros_like(h_lm)
    h_list = []

    print('Solving PDEs...')
    for i in tqdm(range(10000)):
        if i % 1500 == 0:
            height = ISHT(h_lm)
            for _ in range(2):
                lon = np.random.uniform(0.3, 2 * np.pi - 0.3)
                colat = np.random.uniform(0.3, np.pi - 0.3)
                height += np.exp(-((phi_mesh - lon)**2 + (theta_mesh - colat)**2) / (2 * 0.1**2))
            h_lm = SHT(height)
        h_lm, v_lm = simulation_step(h_lm, v_lm)
        if i % 20 == 0:
            h_list.append(ISHT(h_lm))
    return h_list


def time_plot(h_list):
    print('Plotting 3D movie...')
    phi_mesh_ext = np.hstack([phi_mesh, phi_mesh[:, [0]]])[:, :len(phi)//2]
    theta_mesh_ext = np.hstack([theta_mesh, theta_mesh[:, [0]]])[:, :len(phi)//2]
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection='3d')
    ax.relim()
    ax.autoscale_view()
    
    def init():
        return
    
    def animate(i):
        ax.clear()
        h = h_list[i]
        h_ext = np.hstack([h, h[:, [0]]])[:, :len(phi)//2]
        x = h_ext * np.cos(phi_mesh_ext) * np.sin(theta_mesh_ext)
        y = h_ext * np.sin(phi_mesh_ext) * np.sin(theta_mesh_ext)
        z = h_ext * np.cos(theta_mesh_ext)

        ax.plot_wireframe(x, y, z, rstride=3, cstride=3)
        ax.view_init(elev=0, azim=90)
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        ax.set_zlim([-1.5, 1.5])
        ax.set_box_aspect([1, 1, 1])
        ax.set_axis_off()
        plt.title(str(i))
        if i % 100 == 0:
            print(i)
        return
    
    anim = animation.FuncAnimation(fig, animate, init_func = init, frames = 500, interval = 0.03, blit = False)
    anim.save('./test.mp4', fps = 1/0.04, writer = 'ffmpeg', dpi = 200)


if __name__ == "__main__":
    time_plot(simulate())
