import numpy as np
import colorsys
from sgp4.api import Satrec, jday
from datetime import datetime, timedelta
import urllib.request
import tkinter as tk
from tkinter import messagebox
from vpython import sphere, vector, rate, color, scene, label, textures, vec, box

from planetary_constants import earth
earth_radius = earth['radius']

class SatelliteGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Satellite Collision Simulation")

        self.status_label = tk.Label(self, text="Fetching TLE data from multiple sources...")
        self.status_label.grid(row=0, column=0, columnspan=2, pady=5)

        self.tle_sources = [
            ("Active Satellites", "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle"),
            ("Space Stations", "https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle"),
            ("Cosmos 1408 Debris", "https://celestrak.org/NORAD/elements/gp.php?GROUP=cosmos-1408-debris&FORMAT=tle"),
            ("Fengyun-1C Debris", "https://celestrak.org/NORAD/elements/gp.php?GROUP=fengyun-1c-debris&FORMAT=tle"),
            ("Iridium 33 Debris", "https://celestrak.org/NORAD/elements/gp.php?GROUP=iridium-33-debris&FORMAT=tle"),
            ("Cosmos 2251 Debris", "https://celestrak.org/NORAD/elements/gp.php?GROUP=cosmos-2251-debris&FORMAT=tle")
        ]

        self.selected_satellites = set()

        self.tle_data = self.fetch_tle_data()
        if not self.tle_data:
            messagebox.showerror("Error", "Failed to fetch TLE data.")
            return

        self.filter_vars = {}
        self.search_var = tk.StringVar()
        self.filtered_tle_data = self.tle_data.copy()

        self.create_filter_frame()

        self.create_search_bar()

        self.create_satellite_listbox()

        self.create_user_input_fields()

        self.start_button = tk.Button(self, text="Start Simulation", command=self.start_simulation)
        self.start_button.grid(row=self.row_counter, column=0, columnspan=2, pady=10)
        self.row_counter += 1

    def fetch_tle_data(self):
        tle_data = []
        self.grouped_tle_data = {}
        for name, url in self.tle_sources:
            try:
                self.status_label.config(text=f"Fetching TLE data for {name}...")
                self.update_idletasks()
                response = urllib.request.urlopen(url)
                data = response.read().decode('utf-8').splitlines()
                tle_entries = [(f"{name}: {data[i].strip()}", data[i+1], data[i+2], name) for i in range(0, len(data), 3)]
                tle_data.extend(tle_entries)
                self.grouped_tle_data[name] = tle_entries
            except Exception as e:
                print(f"Error fetching TLE data from {name}: {e}")
                continue
        self.status_label.config(text="TLE data fetched successfully.")
        return tle_data

    def create_filter_frame(self):
        self.row_counter = 1
        self.filter_frame = tk.LabelFrame(self, text="Filter by TLE Group")
        self.filter_frame.grid(row=self.row_counter, column=0, columnspan=2, pady=5, sticky='w')
        self.row_counter += 1
        for name, _ in self.tle_sources:
            var = tk.IntVar(value=1)
            cb = tk.Checkbutton(self.filter_frame, text=name, variable=var, command=self.update_filters)
            cb.pack(anchor='w')
            self.filter_vars[name] = var

    def create_search_bar(self):
        search_frame = tk.Frame(self)
        search_frame.grid(row=self.row_counter, column=0, columnspan=2, pady=5, sticky='w')
        self.row_counter += 1
        tk.Label(search_frame, text="Search Satellites:").pack(side='left')
        search_entry = tk.Entry(search_frame, textvariable=self.search_var)
        search_entry.pack(side='left', padx=5)
        search_entry.bind('<KeyRelease>', lambda event: self.update_filters())

    def create_satellite_listbox(self):
        listbox_frame = tk.Frame(self)
        listbox_frame.grid(row=self.row_counter, column=0, columnspan=2, pady=5)
        self.row_counter += 1

        scrollbar = tk.Scrollbar(listbox_frame, orient='vertical')
        scrollbar.pack(side='right', fill='y')

        self.sat_listbox = tk.Listbox(listbox_frame, selectmode='multiple', height=15, width=70, yscrollcommand=scrollbar.set)
        self.sat_listbox.pack(side='left', fill='both')

        scrollbar.config(command=self.sat_listbox.yview)

        self.sat_listbox.bind('<<ListboxSelect>>', self.on_satellite_select)

        self.populate_sat_listbox()

    def create_user_input_fields(self):
        self.start_time_label = tk.Label(self, text="Start Date and Time (YYYY-MM-DD HH:MM:SS UTC):")
        self.start_time_label.grid(row=self.row_counter, column=0, sticky='w', pady=5, padx=5)
        default_start_time = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        self.start_time_entry = tk.Entry(self)
        self.start_time_entry.insert(0, default_start_time)
        self.start_time_entry.grid(row=self.row_counter, column=1, sticky='w', pady=5, padx=5)
        self.row_counter += 1

        self.end_time_label = tk.Label(self, text="End Date and Time (YYYY-MM-DD HH:MM:SS UTC):")
        self.end_time_label.grid(row=self.row_counter, column=0, sticky='w', pady=5, padx=5)
        default_end_time = (datetime.utcnow() + timedelta(hours=24)).strftime("%Y-%m-%d %H:%M:%S")
        self.end_time_entry = tk.Entry(self)
        self.end_time_entry.insert(0, default_end_time)
        self.end_time_entry.grid(row=self.row_counter, column=1, sticky='w', pady=5, padx=5)
        self.row_counter += 1

        self.playback_speed_label = tk.Label(self, text="Playback Speed (e.g., 1x, 2x, 0.5x):")
        self.playback_speed_label.grid(row=self.row_counter, column=0, sticky='w', pady=5, padx=5)
        self.playback_speed_entry = tk.Entry(self)
        self.playback_speed_entry.insert(0, "1x")
        self.playback_speed_entry.grid(row=self.row_counter, column=1, sticky='w', pady=5, padx=5)
        self.row_counter += 1

        self.collision_threshold_label = tk.Label(self, text="Collision Threshold Distance (km):")
        self.collision_threshold_label.grid(row=self.row_counter, column=0, sticky='w', pady=5, padx=5)
        self.collision_threshold_entry = tk.Entry(self)
        self.collision_threshold_entry.insert(0, "50.0")
        self.collision_threshold_entry.grid(row=self.row_counter, column=1, sticky='w', pady=5, padx=5)
        self.row_counter += 1

    def update_filters(self):
        selected_groups = [name for name, var in self.filter_vars.items() if var.get() == 1]
        search_term = self.search_var.get().lower()
        self.filtered_tle_data = [
            (name, tle1, tle2, group) for name, tle1, tle2, group in self.tle_data
            if group in selected_groups and search_term in name.lower()
        ]
        self.populate_sat_listbox()

    def populate_sat_listbox(self):
        self.sat_listbox.delete(0, tk.END)
        for idx, sat_info in enumerate(self.filtered_tle_data):
            sat_name = sat_info[0].strip()
            self.sat_listbox.insert(tk.END, sat_name)
            if sat_name in self.selected_satellites:
                self.sat_listbox.selection_set(idx)

    def on_satellite_select(self, event):
        selection = self.sat_listbox.curselection()
        selected_names_in_view = set(self.sat_listbox.get(i) for i in selection)
        current_view_satellites = set(self.sat_listbox.get(0, tk.END))
        self.selected_satellites -= current_view_satellites
        self.selected_satellites |= selected_names_in_view

    def start_simulation(self):
        if len(self.selected_satellites) < 2:
            messagebox.showwarning("Warning", "Please select at least 2 satellites for collision prediction.")
            return

        selected_satellites = [sat_info for sat_info in self.tle_data if sat_info[0].strip() in self.selected_satellites]
        self.satellite_names = [name.strip() for name, _, _, _ in selected_satellites]
        selected_tles = [(tle1, tle2) for _, tle1, tle2, _ in selected_satellites]

        self.satellites = initialize_satellites(selected_tles)

        start_time_str = self.start_time_entry.get()
        end_time_str = self.end_time_entry.get()
        playback_speed_str = self.playback_speed_entry.get()
        collision_threshold_str = self.collision_threshold_entry.get()

        try:
            self.start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S")
            self.end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            messagebox.showerror("Error", "Invalid date/time format. Please use YYYY-MM-DD HH:MM:SS")
            return

        if self.start_time >= self.end_time:
            messagebox.showerror("Error", "Start time must be before end time.")
            return

        self.delta_t = timedelta(minutes=10)

        try:
            self.collision_threshold = float(collision_threshold_str)
        except ValueError:
            messagebox.showerror("Error", "Invalid collision threshold distance. Please enter a number.")
            return

        try:
            if playback_speed_str.endswith('x'):
                self.playback_speed = float(playback_speed_str[:-1])
            else:
                self.playback_speed = float(playback_speed_str)
        except ValueError:
            messagebox.showerror("Error", "Invalid playback speed. Please enter a number (e.g., 1x, 2x, 0.5x).")
            return

        self.destroy()

        run_simulation(self.satellites, self.satellite_names, self.start_time, self.end_time, self.delta_t, self.collision_threshold, self.playback_speed)

def initialize_satellites(tle_list):
    satellites = []
    for tle1, tle2 in tle_list:
        sat = Satrec.twoline2rv(tle1, tle2)
        satellites.append(sat)
    return satellites

def simulate_and_detect_collisions(satellites, start_time, delta_t, end_time, collision_threshold, satellite_names):
    sat_positions = [[] for _ in range(len(satellites))]
    collision_times = []
    time = start_time

    while time <= end_time:
        jd, fr = jday(time.year, time.month, time.day, time.hour, time.minute, time.second)

        for i, sat in enumerate(satellites):
            e, r, v = sat.sgp4(jd, fr)
            if e == 0:
                sat_positions[i].append(r)
            else:
                sat_positions[i].append([np.nan, np.nan, np.nan])

        for i in range(len(satellites)):
            for j in range(i + 1, len(satellites)):
                if not np.isnan(sat_positions[i][-1][0]) and not np.isnan(sat_positions[j][-1][0]):
                    dist = np.linalg.norm(np.array(sat_positions[i][-1]) - np.array(sat_positions[j][-1]))
                    if dist < collision_threshold:
                        collision_times.append((time, satellite_names[i], satellite_names[j], dist))
        time += delta_t

    return sat_positions, collision_times

def run_simulation(satellites, satellite_names, start_time, end_time, delta_t, collision_threshold, playback_speed):
    sat_positions, collision_times = simulate_and_detect_collisions(
        satellites, start_time, delta_t, end_time, collision_threshold, satellite_names)

    if not all(len(positions) > 0 for positions in sat_positions):
        print("No satellite positions were calculated. Please check your TLE data.")
        return

    if collision_times:
        for collision in collision_times:
            print(f"Collision detected at {collision[0]} between satellite {collision[1]} and {collision[2]} at distance {collision[3]:.2f} km.")
    else:
        print("No collision detected.")

    animate_orbits_vpython(sat_positions, collision_times, start_time, delta_t, playback_speed, satellite_names)

def get_color(index, total):
    hue = index / total
    rgb = colorsys.hsv_to_rgb(hue, 1.0, 1.0)
    return vec(*rgb)

def animate_orbits_vpython(sat_positions, collision_times, start_time, delta_t, playback_speed, satellite_names):
    root = tk.Tk()
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.destroy()

    scene.title = "Orbital Collision Simulation"
    scene.width = screen_width
    scene.height = screen_height
    scene.background = color.black
    scene.autoscale = False
    scene.range = 2 * earth_radius
    scene.center = vector(0, 0, 0)

    earth_obj = sphere(pos=vector(0, 0, 0), radius=earth_radius, texture=textures.earth, shininess=0)

    satellite_objects = []
    colors_list = [get_color(i, len(sat_positions)) for i in range(len(sat_positions))]

    legend_items = []
    legend_y_start = 1.5 * earth_radius
    legend_x = -2.5 * earth_radius
    legend_z = 0 
    legend_spacing = 0.15 * earth_radius

    for i, (positions, name) in enumerate(zip(sat_positions, satellite_names)):
        sat_color = colors_list[i % len(colors_list)]
        initial_pos = positions[0] if not np.isnan(positions[0][0]) else [0, 0, 0]
        sat = sphere(pos=vector(initial_pos[0], initial_pos[1], initial_pos[2]),
                     radius=earth_radius * 0.02, color=sat_color,
                     make_trail=True, trail_type="curve", retain=50)
        satellite_objects.append(sat)

        legend_pos = vector(legend_x, legend_y_start - i * legend_spacing, legend_z)
        color_box = box(pos=legend_pos, size=vector(earth_radius * 0.05, earth_radius * 0.05, earth_radius * 0.05),
                        color=sat_color)
        legend_label = label(pos=legend_pos + vector(earth_radius * 0.08, 0, 0),
                             text=name, height=12, xoffset=0, yoffset=0,
                             border=0, box=False, color=color.white)
        legend_items.append((color_box, legend_label))

    num_frames = len(sat_positions[0])
    times = [start_time + i * delta_t for i in range(num_frames)]

    time_label = label(pos=vector(0, -1.5 * earth_radius, 0), text='', height=16,
                       xoffset=0, yoffset=0, border=0, box=False, color=color.white)

    for frame in range(num_frames):
        rate(int(playback_speed * 10))
        current_time = times[frame]
        time_label.text = f'Time: {current_time.strftime("%Y-%m-%d %H:%M:%S UTC")}'
        for i, sat_obj in enumerate(satellite_objects):
            pos = sat_positions[i][frame]
            if not np.isnan(pos[0]):
                sat_obj.pos = vector(pos[0], pos[1], pos[2])

    while True:
        rate(10)


if __name__ == "__main__":
    app = SatelliteGUI()
    app.mainloop()
