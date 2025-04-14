import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def quaternion_rotate_Sphere(points, q): #function to rotate celestial objects
    q = np.array(q)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]])
    rotated = []
    for v in points:
        v_q = np.array([v[0], v[1], v[2], 0])
        temp = quaternion_multiply(q, quaternion_multiply(v_q, q_conj))
        rotated.append(temp[:3])
    return np.array(rotated)

def quaternion_rotate_plane(x , y , q): #function to rotate orbits
    q = np.array(q)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]])
    v_q = np.array([x, y, 0, 0])
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
    ])

def is_eclipse(x_earth, y_earth, x_moon, y_moon, z_moon, time, latitude, longtitude):
    
    #Calculating sun coordinates
    x_sun = E_earth * Semi_major_earth
    y_sun = E_earth * (1 - E_earth)
    z_sun = 0
    
    longtitude_rotated = longtitude + (2 * math.pi / 86164) * time #taking into a count rotation of the Earth
    
    #Calculating observer coordinates
    x_observer =  R_earth * math.sin(latitude) * math.cos(longtitude_rotated)
    y_observer =  R_earth * math.sin(latitude) * math.sin(longtitude_rotated)
    z_observer =  R_earth * math.cos(latitude)
    
    rx, ry, rz = x_earth + x_observer, y_earth + y_observer, z_observer
    S = np.array([x_sun - rx, y_sun - ry, z_sun - rz]) #Sun to observer distance
    M = np.array([x_moon - rx, y_moon - ry, z_moon - rz]) #Moon to observer distance
    
    cos_theta = np.dot(S, M) / (np.linalg.norm(S) * np.linalg.norm(M)) 
    theta = math.degrees(math.acos(np.clip(cos_theta, -1.0, 1.0))) #angular distance between moon and sun center 

    if abs(z_moon) < 0.5 and theta < 2.0: #Checking only when moon near sun and ecliptic plane
        
        d_se = math.hypot(x_earth - x_sun, y_earth - y_sun) #distance from sun to earth
        d_me = math.sqrt((x_moon - x_earth)**2 + (y_moon - y_earth)**2 + z_moon**2) #distance from moon to earth

        theta_sun = math.degrees(math.atan(R_sol / d_se)) #angular sun radius
        theta_moon = math.degrees(math.atan(R_moon / d_me)) #angular moon radius
        
        if theta_moon >= theta_sun and theta < (theta_moon - theta_sun) * 1.5: #If moon seems larger than sun or equally and difference of radius bigger than distance between them, Total solar eclipse
            return "Total"
        elif theta < (theta_moon + theta_sun) * 1.5: #If distance smaller than the sum of radius of sun and moon (they intersecting), particular solar eclipse
            return "Partial"

    return None

def time_to_date(t):
    base_jd = 2461120.5
    jd = base_jd + t / 86400
    F, I = math.modf(jd)
    I = int(I)
    A = math.floor((I - 1867216.25) / 36524.25)
    if I > 2299160:
        B = I + 1 + A - math.floor(A / 4)
    else:
        B = I
    C = B + 1524
    D = math.floor((C - 122.1) / 365.25)
    E = math.floor(365.25 * D)
    G = math.floor((C - E) / 30.6001)
    day = C - E + F - math.floor(30.6001 * G)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    day_frac = day % 1
    day = int(day)
    hours = int(day_frac * 24)
    minutes = int((day_frac * 24 - hours) * 60)
    return f"{year}-{month:02d}-{day:02d} {hours:02d}:{minutes:02d}"

half_angle_Earth = 23.5/2 #angle of Earth relative to ecliptic plane
half_angle_Moon = 5.145/2 #angle of orbit of moon relative to ecliptic plane

Quaternion_Earth = [math.sin(math.radians(half_angle_Earth)),0,0,math.cos(math.radians(half_angle_Earth))] #quaternion for Earth rotation
Quaternion_Moon = [math.sin(math.radians(half_angle_Moon)),0,0,math.cos(math.radians(half_angle_Moon))] #quaternion for moon orbit rotation

R_moon = math.log(1738, 4186)
R_earth = math.log(6378.137, 4186)
R_sol = math.log(695700, 4186)
E_earth = 0.0167
E_moon = 0.0549
Semi_major_moon = math.log(3.844E5, 10)
W_Earth = 1.99099E-7
W_Moon = 2.66381E-6
Semi_major_earth = math.log(149.598E6, 10)

dt = 86400
t = 0

time_values = []

lat = 24
lon = 74

eclipses = []

#lists for earth coordinates
Earth_X_cords = []
Earth_Y_cords = []

#lists for moon coordinates
Moon_X_cords = []
Moon_Y_cords = []
Moon_Z_cords = []

#angle 2D arrays for sphere point-by-point generation
theta = np.linspace(0, np.pi, 10)
phi = np.linspace(0, 2 * np.pi, 10)
theta, phi = np.meshgrid(theta, phi)

#earth sphere generation
Earth_sphere_x = R_earth * np.sin(theta) * np.cos(phi)
Earth_sphere_y = R_earth * np.sin(theta) * np.sin(phi)
Earth_sphere_z = R_earth * np.cos(theta)

#transform array to look like [[x,y,z],[x,y,z]] 
earth_points = np.stack([
    Earth_sphere_x.flatten(),
    Earth_sphere_y.flatten(),
    Earth_sphere_z.flatten()
], axis=1)

#rotating every point by angle 23.5
rotated_points = quaternion_rotate_Sphere(earth_points, Quaternion_Earth)

#transform to previous shape for surface plotting
Earth_sphere_x = rotated_points[:, 0].reshape(Earth_sphere_x.shape)
Earth_sphere_y = rotated_points[:, 1].reshape(Earth_sphere_y.shape)
Earth_sphere_z = rotated_points[:, 2].reshape(Earth_sphere_z.shape)

#Generate sun sphere with more realistic sun location relative to earth orbit 
Sun_sphere_x = E_earth * Semi_major_earth +  R_sol * np.sin(theta) * np.cos(phi)
Sun_sphere_y = E_earth * (1-E_earth) + R_sol * np.sin(theta) * np.sin(phi)
Sun_sphere_z = R_sol * np.cos(theta)

#Generate moon sphere
Moon_sphere_x = R_moon * np.sin(theta) * np.cos(phi)
Moon_sphere_y = R_moon * np.sin(theta) * np.sin(phi)
Moon_sphere_z = R_moon * np.cos(theta)

#Calculating orbits for Earth and Moon relative to time passed 
while t<(dt*365*5):
    
    t += dt
    time_values.append(t)
    
    phi = W_Earth * t
    theta = W_Moon * t
    
    #Earth cords
    x = Semi_major_earth*math.cos(phi)
    y = Semi_major_earth*math.sqrt(1-E_earth**2)*math.sin(phi)
    
    #Moon cords
    x_moon = x + Semi_major_moon*math.cos(theta)
    y_moon = y + Semi_major_moon*math.sqrt(1-E_moon**2)*math.sin(theta)
    
    
    x_moon, y_moon, z_moon = quaternion_rotate_plane(x_moon, y_moon, Quaternion_Moon) #Rotating moon cords for realistic orbit
    
    Earth_X_cords.append(x)
    Earth_Y_cords.append(y)
    
    Moon_X_cords.append(x_moon)
    Moon_Y_cords.append(y_moon)
    Moon_Z_cords.append(z_moon)
    
    #When near ascending or descending node, cheking more carefully
    if abs(math.sin(theta)) < 0.1:
        
        for t_fine in np.arange(t - 43200, t + 43200, 3600):
            
            phi_f = W_Earth * t_fine
            theta_f = W_Moon * t_fine
            
            x_f = Semi_major_earth * math.cos(phi_f)
            y_f = Semi_major_earth * math.sqrt(1 - E_earth**2) * math.sin(phi_f)
            
            x_mf = x_f + Semi_major_moon * math.cos(theta_f)
            y_mf = y_f + Semi_major_moon * math.sqrt(1 - E_moon**2) * math.sin(theta_f)
            x_mf, y_mf, z_mf = quaternion_rotate_plane(x_mf, y_mf, Quaternion_Moon)
            
            eclipse_type = is_eclipse(x_f, y_f, x_mf, y_mf, z_mf, t_fine, lat, lon)
            
            if eclipse_type:
                
                eclipses.append((t_fine, eclipse_type))
                t += 86400
                break

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

#Orbit lines
line, = ax.plot([], [], [], 'b', alpha=0.5, linewidth=2)
line2, = ax.plot([], [], [], 'r-', alpha=0.5, linewidth=1)

#Drawing Spheres in t = 0
Earth_surface = ax.plot_surface(Earth_sphere_x, Earth_sphere_y, Earth_sphere_z, color='blue', alpha=1)
Sun_surface = ax.plot_surface(Sun_sphere_x, Sun_sphere_y, Sun_sphere_z, color='yellow', alpha=1)
Moon_surface = ax.plot_surface(Moon_sphere_x, Moon_sphere_y, Moon_sphere_z, color='gray', alpha=1)

ax.set(xlim = [-10,10], ylim = [-10,10], zlim = [-8,8], xlabel='x', ylabel='y', zlabel='z')
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
    
    
    time_text.set_text(f"t = {time_values[frame]/86400} days")
    
    #Orbit line behind Earth
    line.set_data(Earth_X_cords[:frame+1], Earth_Y_cords[:frame+1])
    line.set_3d_properties([0] * (frame+1))
    
    #Orbit line behind Moon
    line2.set_data(Moon_X_cords[:frame+1], Moon_Y_cords[:frame+1])
    line2.set_3d_properties(Moon_Z_cords[:frame+1])
    
    #Deleting Earth sphere in old coordinates and drawing new
    Earth_surface.remove()
    Earth_surface = ax.plot_surface(Earth_sphere_x + current_x,Earth_sphere_y + current_y,Earth_sphere_z, color='green', alpha=1)
    
    #Same with Moon
    Moon_surface.remove()
    Moon_surface = ax.plot_surface(Moon_sphere_x + current_x_moon,Moon_sphere_y + current_y_moon,Moon_sphere_z + current_z_moon, color='gray', alpha=1)
    
    return line, line2, Earth_surface, Moon_surface, time_text

print("Predicted Solar Eclipses:")
for t_e, et in eclipses:
    print(f"{time_to_date(t_e)}: {et} Eclipse")

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(time_values), interval=10, blit=False)
plt.show()






