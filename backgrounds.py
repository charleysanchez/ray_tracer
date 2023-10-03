# STRIPES:
def create_plane(position, normal):
    objects.append({'type':'plane', 'position':np.array(position), 
        'normal':np.array(normal),
        'color':lambda M: (color_plane0 
            if (int(M[0] * 2) % 2) else color_plane1),
        'diffuse_c':.75, 'specular_c':1., 'reflection':.25})
    
# CHECKERBOARD:
def create_plane(position, normal):
    objects.append({'type':'plane', 'position':np.array(position), 
        'normal':np.array(normal),
        'color':lambda M: (color_plane0 
            if (int(M[0] * 2) % 2) == (int(M[2] * 2) % 2) else color_plane1),
        'diffuse_c':.75, 'specular_c':1., 'reflection':.25})
    
# SOLID COLOR
def create_plane(position, normal):
    objects.append({'type':'plane', 'position':np.array(position), 
        'normal':np.array(normal),
        'color': color_plane0,
        'diffuse_c':.75, 'specular_c':1., 'reflection':.25})