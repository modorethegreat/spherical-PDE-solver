from sph_harm_transform import *
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
    v_lm_new = (v_lm - dt * g/R * h_lm) / (np.ones((lmax+1, lmax+1)) + dt * b)
    h_lm_new = h_lm + dt * H/R * (ls * (ls + 1))[:, None] * v_lm_new
    return h_lm_new, v_lm_new


h0 = np.ones(theta_mesh.shape) * H
v0 = np.zeros((len(theta), len(phi), 2))
h_lm0 = SHT(h0)
v_lm0 = np.zeros(h_lm0.shape)

h_list = []

h_lm = h_lm0
v_lm = v_lm0

print('Solving PDEs...')

for i in tqdm(range(10000)):
    if i % 1500 == 0:
        h = ISHT(h_lm)
        h += np.exp(-((phi_mesh - np.random.uniform(0.3, np.pi*2 - 0.3))**2 + (theta_mesh - np.random.uniform(0.3, np.pi - 0.3))**2)/0.1**2/2)
        h += np.exp(-((phi_mesh - np.random.uniform(0.3, np.pi*2 - 0.3))**2 + (theta_mesh - np.random.uniform(0.3, np.pi - 0.3))**2)/0.1**2/2)
        # h += np.exp(-((phi_mesh - np.pi * 0.8)**2 + (theta_mesh - np.pi/2)**2)/0.1**2/2)
        h_lm = SHT(h)
    h_lm, v_lm = simulation_step(h_lm, v_lm)
    if i % 20 == 0:
        h_list.append(ISHT(h_lm))


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


time_plot(h_list)