import json
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import os

#Path to json
#Getting folder name
script_dir = os.path.dirname(os.path.abspath(__file__))
build_dir = os.path.join(script_dir, "..", "build")
filename = os.path.join(build_dir, "simulation_log.jsonl")

#check for file existance
if not os.path.exists(filename):
    raise FileNotFoundError(f"Cannot find {filename}")

with open(filename, 'r') as f:
    for line in f:
        print(line)
#list to data storing
time = []
pressure_data = []
Sw_data = []
Sg_data = []
So_data = []
positions_data = []

with open(filename, 'r') as f:
    for line in f:
        entry = json.loads(line)
        t = entry["time"]
        time.append(t)

        
        positions = np.array(entry["positions"])        #array is n*3
        positions_data.append(positions)

        #pressure & sats
        pressure = np.array(entry["fields"]["pressure"])
        Sw = np.array(entry["fields"]["Sw"])
        Sg = np.array(entry["fields"]["Sg"])
        So = 1 - Sw - Sg

        pressure_data.append(pressure)
        Sw_data.append(Sw)
        Sg_data.append(Sg)
        So_data.append(So)

#nampy array
time = np.array(time)

#average value parameter plot
avg_pressure = np.array([np.mean(p) for p in pressure_data])
avg_Sw = np.array([np.mean(s) for s in Sw_data])
avg_Sg = np.array([np.mean(s) for s in Sg_data])
avg_So = np.array([np.mean(s) for s in So_data])

plt.figure(figsize=(10,5))
plt.plot(time, avg_Sw, label="Sw (Water)")
plt.plot(time, avg_Sg, label="Sg (Gas)")
plt.plot(time, avg_So, label="So (Oil)")
plt.xlabel("Time [s]")
plt.ylabel("Average Saturation")
plt.title("Saturations vs Time")
plt.legend()
plt.grid(True)

plt.figure(figsize=(10,5))
plt.plot(time, avg_pressure, label="Pressure")
plt.xlabel("Time [s]")
plt.ylabel("Average Pressure [Pa]")
plt.title("Pressure vs Time")
plt.grid(True)

# scatter 3d plots for last timestep
last_idx = -1
positions = positions_data[last_idx]
Sw = Sw_data[last_idx]
Sg = Sg_data[last_idx]
So = So_data[last_idx]
pressure = pressure_data[last_idx]

fig = plt.figure(figsize=(15,10))

#subplot axes
ax1 = fig.add_subplot(221, projection='3d')
ax2 = fig.add_subplot(222, projection='3d')
ax3 = fig.add_subplot(223, projection='3d')
ax4 = fig.add_subplot(224, projection='3d')

#initial timstep index
t_idx = 0
positions = positions_data[t_idx]
z_plane = positions[:,2]

#scatter plots
sc1 = ax1.scatter(positions[:,0], positions[:,1], z_plane, c=Sw_data[t_idx], cmap='viridis', s=2, alpha=0.4)
ax1.set_title("Sw (Water Saturation)"); ax1.set_xlabel("X"); ax1.set_ylabel("Y"); ax1.set_zlabel("Z")

sc2 = ax2.scatter(positions[:,0], positions[:,1], z_plane, c=Sg_data[t_idx], cmap='viridis', s=2, alpha=0.4)
ax2.set_title("Sg (Gas Saturation)"); ax2.set_xlabel("X"); ax2.set_ylabel("Y"); ax2.set_zlabel("Z")

sc3 = ax3.scatter(positions[:,0], positions[:,1], z_plane, c=So_data[t_idx], cmap='viridis', s=2, alpha=0.4)
ax3.set_title("So (Oil Saturation)"); ax3.set_xlabel("X"); ax3.set_ylabel("Y"); ax3.set_zlabel("Z")

sc4 = ax4.scatter(positions[:,0], positions[:,1], z_plane, c=pressure_data[t_idx], cmap='viridis', s=2, alpha=0.4)
ax4.set_title("Pressure"); ax4.set_xlabel("X"); ax4.set_ylabel("Y"); ax4.set_zlabel("Z")

plt.subplots_adjust(bottom=0.15)  #lspace for slidder

#slidder axis
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Timestep', 0, len(time)-1, valinit=t_idx, valstep=1)

#update func
def update(val):
    t = int(slider.val)
    positions = positions_data[t]
    z_plane = positions[:,2]

    sc1._offsets3d = (positions[:,0], positions[:,1], z_plane)
    sc1.set_array(Sw_data[t])
    
    sc2._offsets3d = (positions[:,0], positions[:,1], z_plane)
    sc2.set_array(Sg_data[t])
    
    sc3._offsets3d = (positions[:,0], positions[:,1], z_plane)
    sc3.set_array(So_data[t])
    
    sc4._offsets3d = (positions[:,0], positions[:,1], z_plane)
    sc4.set_array(pressure_data[t])
    
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.tight_layout()

#2d scatter plot for single layer
plt.figure(figsize=(12,10))

plt.subplot(221)
plt.scatter(positions[:,0], positions[:,1], c=Sw, cmap='Blues', s=2)
plt.title("Sw (Water Saturation)")
plt.xlabel("X"); plt.ylabel("Y")
plt.colorbar(label="Sw")

plt.subplot(222)
plt.scatter(positions[:,0], positions[:,1], c=Sg, cmap='Greens', s=2)
plt.title("Sg (Gas Saturation)")
plt.xlabel("X"); plt.ylabel("Y")
plt.colorbar(label="Sg")

plt.subplot(223)
plt.scatter(positions[:,0], positions[:,1], c=So, cmap='Oranges', s=2)
plt.title("So (Oil Saturation)")
plt.xlabel("X"); plt.ylabel("Y")
plt.colorbar(label="So")

plt.subplot(224)
plt.scatter(positions[:,0], positions[:,1], c=pressure, cmap='Reds', s=2)
plt.title("Pressure")
plt.xlabel("X"); plt.ylabel("Y")
plt.colorbar(label="Pressure")

plt.tight_layout()
plt.show()