import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from datetime import datetime, timedelta

def quaternion_rotate_Sphere(points, q): #function to rotate celestial objects
    q = np.array(q, dtype=float)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)
    rotated = []
    for v in points:
        v_q = np.array([v[0], v[1], v[2], 0.0], dtype=float)
        temp = quaternion_multiply(q, quaternion_multiply(v_q, q_conj))
        rotated.append(temp[:3])
    return np.array(rotated, dtype=float)

def quaternion_rotate_plane(x, y, q): #function to rotate orbits
    q = np.array(q, dtype=float)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)
    v_q = np.array([x, y, 0.0, 0.0], dtype=float)
    rotated = quaternion_multiply(q, quaternion_multiply(v_q, q_conj))
    return rotated[0], rotated[1], rotated[2]

def quaternion_rotate_point(x,y,z,q):
    q = np.array(q, dtype=float)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)
    v_q = np.array([x, y, z, 0.0], dtype=float)
    rotated = quaternion_multiply(q, quaternion_multiply(v_q, q_conj))
    return rotated[0], rotated[1], rotated[2]

def quaternion_multiply(q1, q2): #quaternion multiplication 
    x1, y1, z1, w1 = q1
    x2, y2, z2, w2 = q2
    return np.array([
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
        w1*w2 - x1*x2 - y1*y2 - z1*z2
    ], dtype=float)

def jd_to_datetime(JD):
    """Convert Julian Date to a Python datetime in UTC."""
    days = int(math.floor(JD - 2451545.0))
    secs = (JD - 2451545.0 - days) * 86400.0
    return datetime(2000,1,1,12,0,0) + timedelta(days=days, seconds=secs)

def Calculate_visual_orbit_coordinates(t):
    
    phi = W_Earth * t
    theta = W_Moon * t
    
    #Earth cords with log scale (base 10)
    x_earth_visual = np.log10(Semi_major_earth) * math.cos(phi)
    y_earth_visual = np.log10(Semi_major_earth * math.sqrt(1.0 - E_earth**2)) * math.sin(phi)
    
    #Moon cords with log scale (base 10)
    x_moon_visual = x_earth_visual + np.log10(Semi_major_moon) * math.cos(theta)
    y_moon_visual = y_earth_visual + np.log10(Semi_major_moon * math.sqrt(1.0 - E_moon**2)) * math.sin(theta)

    x_moon_visual, y_moon_visual, z_moon_visual = quaternion_rotate_plane(x_moon_visual, y_moon_visual, Quaternion_Moon) #Rotating moon cords for realistic orbit
    
    return x_earth_visual, y_earth_visual, x_moon_visual, y_moon_visual, z_moon_visual

def Calculate_real_orbits(t):
    
    phi = W_Earth * t
    theta = W_Moon * t
    
    x_earth_real = Semi_major_earth * math.cos(phi)
    y_earth_real = Semi_major_earth * math.sqrt(1.0 - E_earth**2) * math.sin(phi)
      
    x_moon_rel = Semi_major_moon * math.cos(theta)
    y_moon_rel = Semi_major_moon * math.sqrt(1.0 - E_moon**2) * math.sin(theta)
    z_moon_rel = 0.0
    
    #apply inclination of moon
    x_moon_rel_rot, y_moon_rel_rot, z_moon_rel_rot = quaternion_rotate_point(
        x_moon_rel, y_moon_rel, z_moon_rel, Quaternion_Moon
    )
    
    #moon real position
    x_moon_real = x_earth_real + x_moon_rel_rot
    y_moon_real = y_earth_real + y_moon_rel_rot
    z_moon_real = z_moon_rel_rot  # z_earth_real is 0
    
    return x_earth_real, y_earth_real, x_moon_real, y_moon_real, z_moon_real
    
    
half_angle_Earth = 23.5/2.0 #angle of Earth relative to ecliptic plane
half_angle_Moon = 5.145/2.0 #angle of orbit of moon relative to ecliptic plane

Quaternion_Earth = [math.sin(math.radians(half_angle_Earth)), 0.0, 0.0, math.cos(math.radians(half_angle_Earth))] #quaternion for Earth rotation
Quaternion_Moon = [math.sin(math.radians(half_angle_Moon)), 0.0, 0.0, math.cos(math.radians(half_angle_Moon))] #quaternion for moon orbit rotation

R_moon = 1738.0
R_earth = 6378.137
R_sun = 695700.0
E_earth = 0.0167
E_moon = 0.0549
Semi_major_moon = 384400
W_Earth = 1.72021536E-2
W_Moon = 2.30153184E-1
Semi_major_earth = 149.598E6

dt = 0.01
JD_end = 2462586.0
JD_start = 2460939.5000000
t = JD_start

time_values = []

latitude = math.radians(float(input("Enter observer latitude in degrees: ")))
longitude = math.radians(float(input("Enter observer longitude in degrees: ")))

#lists for earth coordinates
Earth_X_cords = []
Earth_Y_cords = []

#lists for moon coordinates
Moon_X_cords = []
Moon_Y_cords = []
Moon_Z_cords = []

#angle 2D arrays for sphere point-by-point generation
theta = np.linspace(0.0, np.pi, 10, dtype=float)
phi = np.linspace(0.0, 2.0 * np.pi, 10, dtype=float)
theta, phi = np.meshgrid(theta, phi)

#earth sphere generation with log scale (base 4186) for radius
log_R_earth = np.log(R_earth) / np.log(4186.0)
Earth_sphere_x = log_R_earth * np.sin(theta) * np.cos(phi)
Earth_sphere_y = log_R_earth * np.sin(theta) * np.sin(phi)
Earth_sphere_z = log_R_earth * np.cos(theta)

#transform array to look like [[x,y,z],[x,y,z]] 
earth_points = np.stack([
    Earth_sphere_x.flatten(),
    Earth_sphere_y.flatten(),
    Earth_sphere_z.flatten()
], axis=1, dtype=float)

#rotating every point by angle 23.5
rotated_points = quaternion_rotate_Sphere(earth_points, Quaternion_Earth)

#transform to previous shape for surface plotting
Earth_sphere_x = rotated_points[:, 0].reshape(Earth_sphere_x.shape)
Earth_sphere_y = rotated_points[:, 1].reshape(Earth_sphere_y.shape)
Earth_sphere_z = rotated_points[:, 2].reshape(Earth_sphere_z.shape)

#Generate sun sphere with more realistic sun location relative to earth orbit 
log_R_sun = np.log(R_sun) / np.log(4186.0)
Sun_sphere_x = log_R_sun * np.sin(theta) * np.cos(phi)
Sun_sphere_y = log_R_sun * np.sin(theta) * np.sin(phi)
Sun_sphere_z = log_R_sun * np.cos(theta)

#Generate moon sphere with log scale (base 4186) for radius
log_R_moon = np.log(R_moon) / np.log(4186.0)
Moon_sphere_x = log_R_moon * np.sin(theta) * np.cos(phi)
Moon_sphere_y = log_R_moon * np.sin(theta) * np.sin(phi)
Moon_sphere_z = log_R_moon * np.cos(theta)

prev_align = False
events = []
#Calculating orbits, checking for eclipses
while t < JD_end:
    
    x_ev, y_ev, x_mv, y_mv, z_mv = Calculate_visual_orbit_coordinates(t)
    
    x_er, y_er, x_mr, y_mr, z_mr = Calculate_real_orbits(t)
    Earth_X_cords.append(x_ev)
    Earth_Y_cords.append(y_ev)
    
    Moon_X_cords.append(x_mv)
    Moon_Y_cords.append(y_mv)
    Moon_Z_cords.append(z_mv)
    
    #Eclipse section
    x_sr = E_earth * Semi_major_earth #x_sun_real
    y_sr = E_earth * (1.0 - E_earth) #y_sun_real
    z_sr = 0.0
    
    #taking into account rotation of the Earth
    longitude_rotated = longitude + (2.0 * math.pi) * t
    
    x_ow = R_earth * math.cos(latitude) * math.cos(longitude_rotated)
    y_ow = R_earth * math.cos(latitude) * math.sin(longitude_rotated)
    z_ow = R_earth * math.sin(latitude)
    
    #rotate observer location with quaternion by 23.5 deg
    x_ow, y_ow, z_ow = quaternion_rotate_point(x_ow, y_ow, z_ow, Quaternion_Earth)
    
    #Observer cords
    x_o =  x_er + x_ow
    y_o =  y_er + y_ow
    z_o =  z_ow
    
    #Observer to objects vector
    o_s = [x_sr - x_o, y_sr - y_o, -z_o] #observer_sun
    o_m = [x_mr - x_o, y_mr - y_o, z_mr - z_o] #observer_moon
    #Distance of vectors between observer, sun and moon
    dos = math.sqrt(o_s[0]**2 + o_s[1]**2 + o_s[2]**2)
    dom = math.sqrt(o_m[0]**2 + o_m[1]**2 + o_m[2]**2)
    
    dot = o_s[0]*o_m[0] + o_s[1]*o_m[1] + o_s[2]*o_m[2]
    sep = math.degrees(math.acos(max(min(dot/(dos*dom), 1.0), -1.0)))
    
    #Angular radius of Sun and Moon
    alpha_s = math.degrees(math.asin(R_sun / dos))
    alpha_m = math.degrees(math.asin(R_moon / dom))
    
    #Correction for Earth radius
    align = (sep <= (alpha_s + alpha_m))
    if align and not prev_align:
        
        #Refine around eclipse
        t0 = max(t - dt, JD_start)
        t1 = min(t + dt, JD_end)
        min_sep = sep
        min_t = t
        
        for i in range(1001):
            
            t_f = t0 + (t1 - t0) * i / 1000.0
            
            x_erf, y_erf, x_mrf, y_mrf, z_mrf = Calculate_real_orbits(t_f)
            
            longitude_rotated_f = longitude + (2.0 * math.pi) * t_f
            x_owf = R_earth * math.cos(latitude) * math.cos(longitude_rotated_f)
            y_owf = R_earth * math.cos(latitude) * math.sin(longitude_rotated_f)
            z_owf = R_earth * math.sin(latitude)
            
            x_owf, y_owf, z_owf = quaternion_rotate_point(x_owf, y_owf, z_owf, Quaternion_Earth)
            
            x_of =  x_erf + x_owf
            y_of =  y_erf + y_owf
            z_of =  z_owf
            
            o_sf = [x_sr - x_of, y_sr - y_of, -z_of]
            o_mf = [x_mrf - x_of, y_mrf - y_of, z_mrf - z_of]
            
            dos_f = math.sqrt(o_sf[0]**2 + o_sf[1]**2 + o_sf[2]**2)
            dom_f = math.sqrt(o_mf[0]**2 + o_mf[1]**2 + o_mf[2]**2)
            
            dot_f = o_sf[0]*o_mf[0] + o_sf[1]*o_mf[1] + o_sf[2]*o_mf[2]
            sep_f = math.degrees(math.acos(max(min(dot_f/(dos_f*dom_f), 1.0), -1.0)))
            
            if sep_f < min_sep:
                min_sep = sep_f
                min_t = t_f
                
        # classify eclipse at min_t
        
        x_erf, y_erf, x_mrf, y_mrf, z_mrf = Calculate_real_orbits(min_t)
        
        longitude_rotated_f = longitude + (2.0 * math.pi) * min_t
        
        x_owf = R_earth * math.cos(latitude) * math.cos(longitude_rotated_f)
        y_owf = R_earth * math.cos(latitude) * math.sin(longitude_rotated_f)
        z_owf = R_earth * math.sin(latitude)
        
        x_owf, y_owf, z_owf = quaternion_rotate_point(x_owf, y_owf, z_owf, Quaternion_Earth)
        
        x_of =  x_erf + x_owf
        y_of =  y_erf + y_owf
        z_of =  z_owf
        
        o_sf = [x_sr - x_of, y_sr - y_of, -z_of]
        o_mf = [x_mrf - x_of, y_mrf - y_of, z_mrf - z_of]
        
        dos_f = math.sqrt(o_sf[0]**2 + o_sf[1]**2 + o_sf[2]**2)
        dom_f = math.sqrt(o_mf[0]**2 + o_mf[1]**2 + o_mf[2]**2)
        
        dot_f = o_sf[0]*o_mf[0] + o_sf[1]*o_mf[1] + o_sf[2]*o_mf[2]
        sep_f = math.degrees(math.acos(max(min(dot_f/(dos_f*dom_f), 1.0), -1.0)))
       
        alpha_s_f = math.degrees(math.asin(R_sun / dos_f))
        alpha_m_f = math.degrees(math.asin(R_moon / dom_f))
        # distance from Earth's center to line (Sun-Moon)
        cross_x = o_sf[1]*o_mf[2] - o_sf[2]*o_mf[1]
        cross_y = o_sf[2]*o_mf[0] - o_sf[0]*o_mf[2]
        cross_z = o_sf[0]*o_mf[1] - o_sf[1]*o_mf[0]
        dist_line = math.sqrt(cross_x**2 + cross_y**2 + cross_z**2) / \
                    math.sqrt((o_mf[0]-o_sf[0])**2 + (o_mf[1]-o_sf[1])**2 + (o_mf[2]-o_sf[2])**2)
                    
        if dist_line <= R_earth and alpha_m_f >= alpha_s_f:
            typ = "Total"
        else:
            typ = "Partial"
            
        events.append((min_t-22, typ))
    prev_align = align
    t += dt
    time_values.append(t)
    
 
for (jd_event, etype) in events:
    t = jd_to_datetime(jd_event)
    print(t.strftime("%Y-%m-%d %H:%M UTC"), etype, "Eclipse")

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

#Orbit lines
line, = ax.plot([], [], [], 'b', alpha=0.5, linewidth=2)
line2, = ax.plot([], [], [], 'r-', alpha=0.5, linewidth=1)

#Drawing Spheres in t = 0
Earth_surface = ax.plot_surface(Earth_sphere_x, Earth_sphere_y, Earth_sphere_z, color='blue', alpha=1)
Sun_surface = ax.plot_surface(Sun_sphere_x, Sun_sphere_y, Sun_sphere_z, color='yellow', alpha=1)
Moon_surface = ax.plot_surface(Moon_sphere_x, Moon_sphere_y, Moon_sphere_z, color='gray', alpha=1)

#Set logarithmic axis limits
ax.set(xlim=[-10, 10], ylim=[-10, 10], zlim=[-8, 8], xlabel='log10(x)', ylabel='log10(y)', zlabel='log10(z)')
time_text = ax.text(10, 10, 5, '', fontsize=15)
ax.grid(True)
ax.legend()

#Function to update frames in animation
def update(frame):
    global Earth_surface, Moon_surface
    current_x = Earth_X_cords[frame]
    current_y = Earth_Y_cords[frame]
    
    current_x_moon = Moon_X_cords[frame] 
    current_y_moon = Moon_Y_cords[frame]
    current_z_moon = Moon_Z_cords[frame]
    
    
    time_text.set_text(f"t = {time_values[frame]} days")
    
    #Orbit line behind Earth
    line.set_data(Earth_X_cords[:frame+1], Earth_Y_cords[:frame+1])
    line.set_3d_properties([0.0] * (frame+1))
    
    #Orbit line behind Moon
    line2.set_data(Moon_X_cords[:frame+1], Moon_Y_cords[:frame+1])
    line2.set_3d_properties(Moon_Z_cords[:frame+1])
    
    #Deleting Earth sphere in old coordinates and drawing new
    Earth_surface.remove()
    Earth_surface = ax.plot_surface(Earth_sphere_x + current_x, Earth_sphere_y + current_y, Earth_sphere_z, color='green', alpha=1)
    
    #Same with Moon
    Moon_surface.remove()
    Moon_surface = ax.plot_surface(Moon_sphere_x + current_x_moon, Moon_sphere_y + current_y_moon, Moon_sphere_z + current_z_moon, color='gray', alpha=1)
    
    return line, line2, Earth_surface, Moon_surface, time_text

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(time_values), interval=10, blit=False)
plt.show()
