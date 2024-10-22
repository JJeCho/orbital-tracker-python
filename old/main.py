import numpy as np
from sgp4.api import Satrec, jday
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import urllib.request
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import simpledialog

from planetary_constants import earth
earth_radius = earth[ 'radius' ]
tle_url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle"

class SatelliteGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Satellite Collision Simulation")
        self.geometry("500x600")
        
        self.status_label = tk.Label(self, text="Fetching TLE data from CelesTrak...")
        self.status_label.pack(pady=10)
        
        self.tle_data = self.fetch_tle_data(tle_url)
        if not self.tle_data:
            messagebox.showerror("Error", "Failed to fetch TLE data.")
            return

        self.sat_listbox = tk.Listbox(self, selectmode='multiple', height=20, width=50)
        self.sat_listbox.pack(pady=10)

        self.populate_sat_listbox()

        self.start_button = tk.Button(self, text="Start Simulation", command=self.start_simulation)
        self.start_button.pack(pady=10)
        
    def fetch_tle_data(self, url):
        try:
            response = urllib.request.urlopen(url)
            data = response.read().decode('utf-8').splitlines()
            return [(data[i], data[i+1], data[i+2]) for i in range(0, len(data), 3)]
        except Exception as e:
            print(f"Error fetching TLE data: {e}")
            return None
    
    def populate_sat_listbox(self):
        for sat_info in self.tle_data:
            self.sat_listbox.insert(tk.END, sat_info[0].strip())  # Add satellite name to the listbox
    
    def start_simulation(self):
        selected_indices = self.sat_listbox.curselection()
        if len(selected_indices) < 2:
            messagebox.showwarning("Warning", "Please select at least 2 satellites for collision prediction.")
            return
        
        selected_tles = [(self.tle_data[i][1], self.tle_data[i][2]) for i in selected_indices]

        satellites = initialize_satellites(selected_tles)

        start_time = datetime.utcnow()
        delta_t = timedelta(minutes=10)
        end_time = start_time + timedelta(hours=24)

        collision_threshold = 50.0

        sat_positions, collision_times = simulate_and_detect_collisions(satellites, delta_t, end_time, collision_threshold)

        if all(len(positions) > 0 for positions in sat_positions):
            if collision_times:
                for collision in collision_times:
                    self.status_label.config(text=f"Collision detected at {collision[0]} between satellite {collision[1]} and {collision[2]}.")
            else:
                self.status_label.config(text="No collision detected.")
            
            animate_orbits(sat_positions, collision_time=len(sat_positions[0]) if collision_times else None)
        else:
            messagebox.showerror("Error", "No satellite positions were calculated. Please check your TLE data.")

def initialize_satellites(tle_list):
    satellites = []
    for tle1, tle2 in tle_list:
        sat = Satrec.twoline2rv(tle1, tle2)
        satellites.append(sat)
    return satellites

def simulate_and_detect_collisions(satellites, delta_t, end_time, collision_threshold):
    start_time = datetime.utcnow()
    sat_positions = [[] for _ in range(len(satellites))]
    collision_times = []
    
    time = start_time
    collision_detected = False

    while time < end_time and not collision_detected:
        jd, fr = jday(time.year, time.month, time.day, time.hour, time.minute, time.second)

        for i, sat in enumerate(satellites):
            e, r, v = sat.sgp4(jd, fr)
            if e == 0:
                sat_positions[i].append(r)

        for i in range(len(satellites)):
            for j in range(i + 1, len(satellites)):
                if len(sat_positions[i]) > 0 and len(sat_positions[j]) > 0:
                    dist = np.linalg.norm(np.array(sat_positions[i][-1]) - np.array(sat_positions[j][-1]))
                    if dist < collision_threshold:
                        collision_times.append((time, i, j, dist))
                        collision_detected = True
                        break
            if collision_detected:
                break
        
        time += delta_t
    
    return sat_positions, collision_times

def animate_orbits(sat_positions, collision_time=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)
    ax.plot_surface(x, y, z, color='b', alpha=0.6)

    scatters = [ax.scatter([], [], [], label=f'Satellite {i+1}') for i in range(len(sat_positions))]

    ax.set_xlim([-2 * earth_radius, 2 * earth_radius])
    ax.set_ylim([-2 * earth_radius, 2 * earth_radius])
    ax.set_zlim([-2 * earth_radius, 2 * earth_radius])
    ax.legend()

    def update(frame):
        for i, scatter in enumerate(scatters):
            if frame < len(sat_positions[i]):
                pos = sat_positions[i][frame]
                x_pos, y_pos, z_pos = [pos[0]], [pos[1]], [pos[2]]
                scatter._offsets3d = (x_pos, y_pos, z_pos)
        return scatters

    total_frames = min([len(p) for p in sat_positions])
    if collision_time:
        total_frames = collision_time

    ani = FuncAnimation(fig, update, frames=total_frames, interval=100)
    plt.show()

if __name__ == "__main__":
    app = SatelliteGUI()
    app.mainloop()