import tkinter as tk
from tkinter import messagebox
import numpy as np
import os
import matplotlib.pyplot as plt 
import tkinter.colorchooser as cc
from tkinter import ttk



image_file = os.path.join(os.path.dirname(__file__), 'example.png')

# create the main window
root = tk.Tk()
root.title('Creation Screen')
root.attributes('-fullscreen', True)
# create the input fields and labels
#style = ttk.Style()
#style.theme_use('classic')
#style.configure("BW.TLabel", foreground="black", background="blue")


blurb_label = tk.Label(root, text='''This prompt is designed to ask you (the user) how many spheres you would like to create to 
be used within the ray tracing algorithm. The range is recommended between 1 and 3 for your 
own sake, more than 3 will be a bit of a mess on the screen, but I have not disallowed you 
to do so, so you can actually enter an infinite number of spheres and face the consequences.
To keep adding spheres, simply click the 'new sphere' button. Once you are done adding spheres,
AND inputting your dimensions for the output, click the 'Done' button. 
                                     
Below I have incorporated an example image with 3 spheres to get you an idea of how that looks.
''')
blurb_label.pack()


image = tk.PhotoImage(file=image_file)
image_label = tk.Label(root, image=image)
image_label.pack()


image_filenames = ['stripes.png', 'checkerboard.png', 'solid.png', 'disco.png']

# Load the images as PhotoImage objects
images = [tk.PhotoImage(file=filename) for filename in image_filenames]


sphere_centers = []
sphere_radii = []
sphere_colors = []
plane_color_map = []



# define a function to retrieve the input values and store them as variables
    
    
def save_inputs2(center_entry_x, center_entry_y, center_entry_z, radius_entry, win):
    try:
        center_array = np.array([center_entry_x.get(), center_entry_y.get(), center_entry_z.get()], dtype=float)
        for i in range(len(center_array)):
                if not float(center_array[i]) and center_array[i]!=0:
                    raise ValueError
        sphere_centers.append(center_array)
        radius = float(radius_entry.get())
        if radius <= 0:
            raise ValueError
        sphere_radii.append(radius)
        choose_color()
        win.destroy()
    except ValueError:
        messagebox.showerror('Error', 'Invalid input, try entering a radius greater than 0!')
        return


def choose_color():
    color = str(cc.askcolor()[0]).replace('(',', ').replace(')',', ' ).split(', ')
    col_array = np.array(color)
    num_col = np.array(col_array[1:-1], dtype=float)/255
    sphere_colors.append(num_col)

def choose_plane_color():
    if chosen_plane_index == 4:
        plane_color_map.append(np.array([0,0,0]))
        plane_win.destroy()
        return
    elif chosen_plane_index == 3:
        color1 = str(cc.askcolor()[0]).replace('(',', ').replace(')',', ' ).split(', ')
        col1_array = np.array(color1)
        num_col1 = np.array(col1_array[1:-1], dtype=float)/255
        plane_color_map.append(num_col1)
        plane_win.destroy()
        return
    else:
        color1 = str(cc.askcolor()[0]).replace('(',', ').replace(')',', ' ).split(', ')
        color2 = str(cc.askcolor()[0]).replace('(',', ').replace(')',', ' ).split(', ')
        col1_array = np.array(color1)
        col2_array = np.array(color2)
        num_col1 = np.array(col1_array[1:-1], dtype=float)/255
        num_col2 = np.array(col2_array[1:-1], dtype=float)/255
        plane_color_map.append(num_col1)
        plane_color_map.append(num_col2)
        plane_win.destroy()

def sphere_params():
    win = tk.Toplevel(root)
    win.title('Sphere Descriptors')
    descriptor_label = tk.Label(win, text='''
    Here you are required to enter values for the specifics of the spheres you want in the created image.
    
    First, you need to place the sphere. Do this by choosing the center coordinates of each sphere.

    Second, you need to specify the radius of the sphere. This is simpler, just choose a radius between 0 and 1 depending on what size of sphere you want.

    Once you have input parameters, a color selection wheel will pop up- this will be the color of the sphere.

    For a reference, the example image seen below has 3 spheres. Sphere 1 has coordinates: (1, .6, 1) and radius: .6. Sphere 2 has coordinates: (-.5, .3, 1.5)
    and radius: .6. Sphere 3 has coordinates (-2, .3, 2) and radius .6. As you can see, the radius is not the only thing the size is dependent on. The camera 
    is located as coordinates (0., 0.35, -1.), so the closer you get to this value, the closer the sphere will appear as well.
    ''')
    descriptor_label.pack()
    ex_image = tk.PhotoImage(file=image_file)
    ex_image_label = tk.Label(win, image=ex_image)
    ex_image_label.pack()

    
    center_label_x = tk.Label(win, text=f'X Coordinate of Sphere: ', fg='brown')
    center_label_x.pack()
    center_entry_x = tk.Scale(win, from_=-4, to=4, orient=tk.HORIZONTAL, resolution=0.1, length=500, sliderlength=20, width=20, troughcolor="lightgray", activebackground="lightblue", highlightthickness=0)
    center_entry_x.pack()

    center_label_y = tk.Label(win, text=f'Y Coordinate of Sphere: ', fg='red')
    center_label_y.pack()
    center_entry_y = tk.Scale(win, from_=0, to=4, orient=tk.HORIZONTAL, resolution=0.1, length=300, sliderlength=20, width=20, troughcolor="lightgray", activebackground="lightblue", highlightthickness=0)
    center_entry_y.pack()

    center_label_z = tk.Label(win, text=f'Z Coordinate of Sphere: ', fg='blue')
    center_label_z.pack()
    center_entry_z = tk.Scale(win, from_=0, to=4, orient=tk.HORIZONTAL, resolution=0.1, length=300, sliderlength=20, width=20, troughcolor="lightgray", activebackground="lightblue", highlightthickness=0)
    center_entry_z.pack()

    radius_label = tk.Label(win, text=f'Radius of Sphere: ')
    radius_label.pack()
    radius_entry = tk.Scale(win, from_=0, to=4, orient=tk.HORIZONTAL, resolution=0.1, length=200, sliderlength=20, width=20, troughcolor="lightgray", activebackground="lightblue", highlightthickness=0)
    radius_entry.pack()
    save_button = tk.Button(win, text='Save Inputs', command=lambda: save_inputs2(center_entry_x, center_entry_y, center_entry_z, radius_entry, win))
    save_button.pack()

    win.mainloop()


# create a button to save the inputs
def save_button():
    save_but = tk.Button(root, text="Add Sphere", command=sphere_params)
    save_but.pack()

#def done_button():
#    done_but = tk.Button(root, text='Done!', command=quit)
#    done_but.pack()

def quit():
    root.destroy()

def height_width_button():
    x_label = tk.Label(root, text="Enter a width for your output image in pixels (Recommended 500):")
    x_label.pack()
    x_entry = tk.Entry(root)
    x_entry.pack()

    y_label = tk.Label(root, text="Enter a height for your output image in pixels (Recommended 400):")
    y_label.pack()
    y_entry = tk.Entry(root)
    y_entry.pack()

    height_width_button = tk.Button(root, text='DONE', command =lambda: save_inputs(x_entry, y_entry))
    height_width_button.pack()

def choose_plane_button():
    plane_button = tk.Button(root, text='Choose your style of plane', command=choose_plane)
    plane_button.pack()

def choose_plane():
    global image_label_2
    global plane_win
    global chosen_plane
    plane_color_map.clear()
    plane_win = tk.Toplevel(root)
    plane_win.title('Customize Your Plane')
    
    info_label = tk.Label(plane_win, text='''
    Here is a list of different images. To view the images, click the one you are interested
    in, then click update image to see what it looks like. The "floor" of each image is different
    and you are free to choose which one looks most interesting to you. Once you are done selecting
    which one you like, click the "Confirm your selection" button and you will be prompted with
    two different color selector tools. The colors will replace the black and white, with
    the first one replacing the white and the second one replacing the black. Then the window will 
    close and your selections will be saved.
    ''')
    info_label.pack()
    # Create a Listbox widget with the image filenames as options
    listbox = tk.Listbox(plane_win)
    for filename in image_filenames:
        listbox.insert(tk.END, filename)
    listbox.pack()

    # Create a Label widget to display the selected image
    image_label_2 = tk.Label(plane_win)
    image_label_2.pack()

    listbox.bind("<<ListboxSelect>>", lambda event: update_image(listbox, images))

 #   update_button = tk.Button(plane_win, text='Update Image', command=lambda: update_image(listbox, images))
 #   update_button.pack()

    confirm_button = tk.Button(plane_win, text='Confirm your selection', command=choose_plane_color)
    confirm_button.pack()

    plane_win.mainloop()

# Define a function to update the image displayed in the Label widget
def update_image(listbox, images):
    global chosen_plane_index
    # Get the index of the selected item in the Listbox

    index = listbox.curselection()[0]
    # Update the Label widget with the selected image
    image_label_2.config(image=images[index])
    if index == 0:
        chosen_plane_index = 1
    elif index == 1:
        chosen_plane_index = 2
    elif index == 2:
        chosen_plane_index = 3
    elif index == 3:
        chosen_plane_index = 4

def save_inputs(x_entry, y_entry):
    global X 
    global Y 
    try:
        X = int(x_entry.get())
    except ValueError:
        messagebox.showerror('Error', 'Invalid input (try entering integer for X and Y)')
        return

    try: 
        Y = int(y_entry.get())
    except ValueError:
        messagebox.showerror('Error', 'Invalid input (try entering integer for X and Y)')
        return

    try:
        if not sphere_centers:
            raise ValueError
    except ValueError:
        messagebox.showerror('Error', 'Try adding a sphere!')
    try:
        if not plane_color_map:
            raise ValueError
    except ValueError:
        messagebox.showerror('Error', 'Try customizing the plane!')

    if sphere_centers and plane_color_map:
        root.destroy()
    



save_button()
choose_plane_button()
#done_button()
height_width_button()

note = tk.Label(root, text='Note: Upon clicking \'Done\', the program will be unresponsive as the image renders. Simply let some time go by!')
note.pack() 
# run the main event loop
root.mainloop()

# easy, just define normalization of vector by using numpy
# takes a vector as the argument
def normalize_vector(v):
    return v / np.linalg.norm(v)



# define camera (viewing center)
camera = np.array([0, 0.35, -1.])

# define origin
origin = np.array([0., 0., 0.])

# define aspect ratio
asp_rat = float(X) / Y

# define screen center
screen = (-1., -1. / asp_rat + .25, 1., 1. / asp_rat + .25)

# define light
light = np.array([5., 5., -10.])
light_color = np.array([1, 1, 1])



# create black/blank canvas
img = np.zeros((Y, X, 3))

depth_max = 5

# define default lighting/material colors
ambient = 0.05
diffuse_c = 1.
specular_c = 1.
specular_k = 50

# if I add more shapes later, it's good to have a general intersection function to support this
def intersection(origin, direction, object):
    if object['type'] == 'plane':
        return plane_intersection(origin, direction, object['position'], object['normal'])
    elif object['type'] == 'sphere':
        return sphere_intersection(origin, direction, object['center'], object['radius'])
    
# ray equation = Origin + tDirection
# intersection must pass through sphere at 2 points- at some points, ||O + tD - center||^2 = radius^2
# ||D||^2*t^2 + 2t(D np.dot O - center) + ||O - center||^2 - r^2 = 0
# can just solve for t by finding coefficients
# ray equation = Origin + tDirection
# intersection must pass through sphere at 2 points- at some points, ||O + tD - center||^2 = radius^2
# ||D||^2*t^2 + 2t(D np.dot O - center) + ||O - center||^2 - r^2 = 0
# can just solve for t by finding coefficients
def sphere_intersection(origin, direction, center, radius):
    a = np.dot(direction, direction)
    b = 2 * np.dot(direction, origin - center)
    c = np.dot(origin - center, origin - center) - radius ** 2
    discriminant = b ** 2 - 4 * a * c
    if discriminant > 0:
        if b > 0:
            t0 = (-b + np.sqrt(discriminant)) / (2 * a)
        else:
            t0 = (-b - np.sqrt(discriminant)) / (2 * a)
        t1 = t0 / a
        t2 = c / t0
        tmin, tmax = min(t1, t2), max(t1, t2)
        # if both exist, we want the minimum UNLESS the minimum is negative because we don't care about the ray moving backwards
        # check tmax first because if tmax < 0, then tmin < 0 as well
        if tmax >= 0:
            if tmin < 0:
                return tmax
            else:
                return tmin
    # if neither times are positive and real- we return np.inf to show that the light never reaches the object, there is no intersection
    return np.inf

def plane_intersection(origin, direction, position, normal):
    # Return the distance from origin to the intersection of the ray (origin, direction) with the 
    # plane (position, normal), or +inf if there is no intersection.
    # origin and plane are 3D points, direction and normal are normalized vectors.
    denom = np.dot(direction, normal)
    if np.abs(denom) < 1e-6:
        return np.inf
    distance = np.dot(position - origin, normal) / denom
    if distance < 0:
        return np.inf
    return distance

def get_color(obj, point_of_intersection):
    color = obj['color']
    if not hasattr(color, '__len__'):
        color = color(point_of_intersection)
    return color

def get_normal(obj, point_of_intersection):
    # Find normal.
    if obj['type'] == 'sphere':
        normal = normalize_vector(point_of_intersection - obj['center'])
    elif obj['type'] == 'plane':
        normal = obj['normal']
    return normal

def ray_trace(origin, direction):
    # first check if there are any intersections at any time
    t = np.inf
    for i, obj in enumerate(objects):
        t_int = intersection(origin, direction, obj)
        if t_int < t:
            t, index = t_int, i
    
    # if there is no intersection at any time then return nothing
    if t == np.inf:
        return
    
    # now we have the object being intersected and the time of intersection we can continue:
    intersected_object = objects[index]
    point_of_intersection = origin + direction * t

    normal_intersection_to_center = get_normal(intersected_object, point_of_intersection)
    color = get_color(intersected_object, point_of_intersection)

    normal_intersection_to_light = normalize_vector(light - point_of_intersection)
    normal_intersection_to_camera = normalize_vector(camera - point_of_intersection)


    
    # now we need to see if the object is behind another one or in its shadow at all
    # do this by checking if the light also intersects another object or not
    # need to check for shadows using this point otherwise object will think it is in its own shadow
    # runs through all the other objects and checks if theres an intersection from the surface of our first object
    # back through the direction of the light
    # checks if the list contains any values, then checks if any of them are less than infinity (they exist)
    # if there is an intersection point, this means there is an object between the ray and the inital object
    # this means that the color should stay black, so we return before computing the color of the object at this point
    other_intersection_points = [intersection(point_of_intersection + normal_intersection_to_center*0.0001, normal_intersection_to_light, diff_obj) for k, diff_obj in enumerate(objects) if k != index]
    if other_intersection_points and min(other_intersection_points) < np.inf:
        return
    
    # now we start finding color of object
    col_ray = ambient
    col_ray += intersected_object.get('diffuse_c', diffuse_c) * max(np.dot(normal_intersection_to_center, normal_intersection_to_light), 0) * color
    col_ray += intersected_object.get('specular_c', specular_c) * max(np.dot(normal_intersection_to_center, normalize_vector(normal_intersection_to_light + normal_intersection_to_camera)), 0) ** specular_k * light_color
    return intersected_object, point_of_intersection, normal_intersection_to_center, col_ray


# starting color (black)



def create_sphere(center, radius, color):
    objects.append({'type': 'sphere', 'center':np.array(center), 
        'radius':np.array(radius), 'color':np.array(color), 'reflection':.5})


objects = []

def get_plane_color(M):
    if chosen_plane_index == 2:
        return plane_color_map[0] if (int(M[0] * 2) % 2) == (int(M[2] * 2) % 2) else plane_color_map[1]
    elif chosen_plane_index == 1:
        return plane_color_map[0] if (int(M[0] * 2) % 2) else plane_color_map[1]
    elif chosen_plane_index == 3:
        return plane_color_map[0]
    elif chosen_plane_index == 4:
        return (np.array([np.sin(M[0]*10) * np.cos(M[2]*10), np.sin(M[0]*10) * np.sin(M[2]*10), np.cos(M[0]*10)]) * 0.5 + 0.5)

objects.append({'type': 'plane',
                'position': np.array([0, -.5, 0.]),
                'normal': np.array([0., 1., 0.]),
                'color': lambda M: get_plane_color(M),
                'diffuse_c': .75,
                'specular_c': 1.,
                'reflection': .25})

    
# List of objects.


for i in range(len(sphere_centers)):
    create_sphere(sphere_centers[i], sphere_radii[i], sphere_colors[i])



def create_image():
    color = np.zeros(3)
    for i, x in enumerate(np.linspace(screen[0], screen[2], X)):
        for j, y in enumerate(np.linspace(screen[1], screen[3], Y)):
            color[:] = 0
            origin[:2] = (x, y)
            depth = 0
            ray_origin = camera
            ray_direction = normalize_vector(origin - camera)
            reflection = 1.
            while depth < depth_max:
                traced = ray_trace(ray_origin, ray_direction)
                if not traced:
                    break
                intersected_object, point_of_intersection, normal_intersection_to_center, col_ray = traced
                ray_origin, ray_direction = point_of_intersection + normal_intersection_to_center*1e-4, normalize_vector(ray_direction - 2 * np.dot(ray_direction, normal_intersection_to_center) * normal_intersection_to_center)
                depth += 1
                color += reflection * col_ray
                        
                reflection *= intersected_object.get('reflection', 1.)
            img[Y - j - 1, i, :] = np.clip(color, 0, 1)



    plt.imsave('fig.png', img)

create_image()
