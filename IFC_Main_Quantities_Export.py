import ifcopenshell.util.element
from ifcopenshell import geom
import numpy as np
from mathutils import Vector

model = ifcopenshell.open('C:/path/to/IFC_file.ifc')
exportedCsvFilePathAndName = 'C:/path/to/my_csv.csv'


def get_volume(geometry):
	"""_summary_: Returns the volume of given geometry. Using the signed volume of a tetrahedron technic
	"""
    def signed_triangle_volume(p1, p2, p3):
        v321 = p3[0] * p2[1] * p1[2]
        v231 = p2[0] * p3[1] * p1[2]
        v312 = p3[0] * p1[1] * p2[2]
        v132 = p1[0] * p3[1] * p2[2]
        v213 = p2[0] * p1[1] * p3[2]
        v123 = p1[0] * p2[1] * p3[2]
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)

    verts = geometry.verts
    faces = geometry.faces
    grouped_verts = [[verts[i], verts[i + 1], verts[i + 2]] for i in range(0, len(verts), 3)]
    volumes = [
        signed_triangle_volume(grouped_verts[faces[i]], grouped_verts[faces[i + 1]], grouped_verts[faces[i + 2]])
        for i in range(0, len(faces), 3)
    ]
    return abs(sum(volumes))


def get_area_vf(vertices, faces):
	"""_summary_: Returns the area of given vertices and faces. Used by other function to get specific areas by giving filtered faces as parameter.
	"""
    try:
        # Calculate the triangle normal vectors
        v1 = vertices[faces[:, 1]] - vertices[faces[:, 0]]
        v2 = vertices[faces[:, 2]] - vertices[faces[:, 0]]
        triangle_normals = np.cross(v1, v2)

        # Normalize the normal vectors to get their length (i.e., triangle area)
        triangle_areas = np.linalg.norm(triangle_normals, axis=1) / 2

        # Sum up the areas to get the total area of the mesh
        mesh_area = np.sum(triangle_areas)
    except:
        return 0

    return mesh_area


def get_side_area(geometry, axis, bbox):
	"""_summary_: Returns the side area of given geometry.
	:param geometry: ifcopenshell shape.geometry
	:param axis: main axis of geometry
	:param bbox: oriented bounded box
	:return str: main axis x or y
	"""
    verts = geometry.verts
    faces = geometry.faces
    vertices = np.array([[verts[i], verts[i + 1], verts[i + 2]] for i in range(0, len(verts), 3)])
    faces = np.array([[faces[i], faces[i + 1], faces[i + 2]] for i in range(0, len(faces), 3)])

    # Calculate the triangle normal vectors
    v1 = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    v2 = vertices[faces[:, 2]] - vertices[faces[:, 0]]
    triangle_normals = np.cross(v1, v2)

    # Normalize the normal vectors
    triangle_normals = triangle_normals / np.linalg.norm(triangle_normals, axis=1)[:, np.newaxis]

    # Get side normal of bbox
    side_normal = get_side_normal(bbox)
    
    # Find the faces with a normal vector aligned with bbox side normal
    filtered_faces = []
    for i in range(len(faces)):
        if np.dot(side_normal, triangle_normals[i]) >.98:
            filtered_faces.append(faces[i])
    
    filtered_faces = np.array(filtered_faces)
    return get_area_vf(vertices, filtered_faces)


def get_footprint_area(geometry):
	"""_summary_: Returns the footprint area of given geometry.
	"""
    verts = geometry.verts
    faces = geometry.faces
    vertices = np.array([[verts[i], verts[i + 1], verts[i + 2]] for i in range(0, len(verts), 3)])
    faces = np.array([[faces[i], faces[i + 1], faces[i + 2]] for i in range(0, len(faces), 3)])

    # Calculate the triangle normal vectors
    v1 = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    v2 = vertices[faces[:, 2]] - vertices[faces[:, 0]]
    triangle_normals = np.cross(v1, v2)

    # Normalize the normal vectors
    triangle_normals = triangle_normals / np.linalg.norm(triangle_normals, axis=1)[:, np.newaxis]

    # Find the faces with a normal vector pointing in the desired +Z normal direction
    filtered_face_indices = np.where(triangle_normals[:, 2] > tol)[0]
    filtered_faces = faces[filtered_face_indices]
    return get_area_vf(vertices, filtered_faces)


def get_oriented_bound_box(element):
	"""_summary_: Returns the oriented bounded box of an IFC element. Rotation only on Z axis.
	"""
    settings = geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    shape = geom.create_shape(settings , element)
    geometry = shape.geometry
    verts = geometry.verts
	# Convert vertices to 2D
    verts = [[verts[i], verts[i + 1]] for i in range(0, len(verts), 3)]

    points = np.asarray(verts)
    means = np.mean(points, axis=1)

    cov = np.cov(points, y = None,rowvar = 0,bias = 1)

    v, vect = np.linalg.eig(cov)

    tvect = np.transpose(vect)
    points_r = np.dot(points, np.linalg.inv(tvect))

    co_min = np.min(points_r, axis=0)
    co_max = np.max(points_r, axis=0)

    xmin, xmax = co_min[0], co_max[0]
    ymin, ymax = co_min[1], co_max[1]
	# Get Z min and max for later
    zmin = min([geometry.verts[i + 2] for i in range(0, len(geometry.verts), 3)])
    zmax = max([geometry.verts[i + 2] for i in range(0, len(geometry.verts), 3)])

    xdif = (xmax - xmin) * 0.5
    ydif = (ymax - ymin) * 0.5

    cx = xmin + xdif
    cy = ymin + ydif
    
    corners = np.array([
        [cx - xdif, cy - ydif],
        [cx - xdif, cy + ydif],
        [cx + xdif, cy - ydif],
        [cx + xdif, cy + ydif],
    ])
    
    corners = np.dot(corners, tvect)
    center = np.dot([cx, cy], tvect)

	# 2D oriented bounded box
    corners = [(el[0], el[1]) for el in corners]
    
	# Convert 2D to 3D by adding Z min / max.
    corners3d = [(corners[0][0], corners[0][1], zmin),
        (corners[0][0], corners[0][1], zmax),
        (corners[1][0], corners[1][1], zmin),
        (corners[1][0], corners[1][1], zmax),
        (corners[2][0], corners[2][1], zmin),
        (corners[2][0], corners[2][1], zmax),
        (corners[3][0], corners[3][1], zmin),
        (corners[3][0], corners[3][1], zmax)]
    
    return corners3d








def get_object_main_axis(bbox):
	"""_summary_: Returns the main object axis. Useful for profile-defined objects. Adapted so it ignores Z axis.
	:param blender-object o: Blender Object
	:return str: main axis x or y
	"""
	x = (Vector(bbox[4]) - Vector(bbox[0])).length
	y = (Vector(bbox[2]) - Vector(bbox[0])).length

	if x >= y:
		return "x"
	else:
		return "y"


def get_linear_length(bbox):
	"""_summary_: Returns the linear length of a bbox. Used for wall and beam length.
	"""
	x = (Vector(bbox[4]) - Vector(bbox[0])).length
	y = (Vector(bbox[2]) - Vector(bbox[0])).length
	
	if get_object_main_axis(bbox) == "x":
		return x
	elif get_object_main_axis(bbox) == "y":
		return y
	else:
		return 0


def get_height((bbox):
	"""_summary_: Returns the height of a bbox. Could be used for wall height if needed.
	"""
	return (Vector(bbox[1]) - Vector(bbox[0])).length


def get_side_normal(bbox):
	"""_summary_: Returns the side normal vector of given bbox. Used to get correct quantities for non axis-aligned geometry.
	"""
    z = Vector((0,0,1))
    if get_object_main_axis(bbox) == "x":
        x = Vector(bbox[4]) - Vector(bbox[0])
        normal = np.cross(x, z)
    elif get_object_main_axis(bbox) == "y":
        y = Vector(bbox[2]) - Vector(bbox[0])
        normal = np.cross(y, z)
    
    normal = normal / np.linalg.norm(normal)
    
    return normal


def get_parent_building(entity):
	"""_summary_: Retrieve whatever Building contains this entity, or None
	"""
    if entity.is_a("IfcElement"):
        parents = entity.ContainedInStructure
        decomposes = entity.Decomposes
        if not parents:
            if not decomposes:
                return None
            parent = decomposes[0].RelatingObject
        else:
            parent = parents[0].RelatingStructure
    elif entity.is_a("IfcSpatialElement"):
        decomposes = entity.Decomposes
        if not decomposes:
            return None
        parent = decomposes[0].RelatingObject
    elif entity.is_a("IfcStructuralItem"):
        assignments = entity.HasAssignments
        if not assignments:
            return None
        parent = assignments[0].RelatingGroup
    elif entity.is_a("IfcSystem"):
        services = entity.ServicesBuildings
        if not services:
            return None
        parent = services[0].RelatedBuildings[0]
    else:
        return None
    if parent.is_a("IfcBuilding"):
        return parent
    return get_parent_building(parent)


def get_my_linear(element):
	"""_summary_: Returns the linear length value that will be stored into CSV file. Returns empty string if not concerned.
	"""
    ml = ""
    if element.is_a("IfcWall") or element.is_a("IfcBeam"):
        bbox = get_oriented_bound_box(element)
        ml = str(get_linear_length(bbox))
	# Replace dots with commas for Excel
    return ml.replace(".",",")


def get_my_net_area(element):
	"""_summary_: Returns the net area value that will be stored into CSV file. Returns empty string if not concerned.
	"""
    m2 = ""
    try:
        if element.is_a("IfcWall"):
            settings = geom.settings()
            settings.set(settings.USE_WORLD_COORDS, True)
            shape = geom.create_shape(settings , element)
            m2 = str(get_side_area(shape.geometry, get_object_main_axis(get_oriented_bound_box(element)),get_oriented_bound_box(element)))
        if element.is_a("IfcSlab"):
            settings = geom.settings()
            settings.set(settings.USE_WORLD_COORDS, True)
            shape = geom.create_shape(settings , element)
            m2 = str(get_footprint_area(shape.geometry))
    except:
        print(element.GlobalId)
	# Replace dots with commas for Excel
    return m2.replace(".",",")


def get_my_gross_area(element):
	"""_summary_: Returns the gross area value that will be stored into CSV file. Returns empty string if not concerned.
	"""
    m2 = ""
    try:
        if element.is_a("IfcWall"):
            settings = geom.settings()
            settings.set(settings.USE_WORLD_COORDS, True)
            settings.set(settings.DISABLE_OPENING_SUBTRACTIONS, True)
            shape = geom.create_shape(settings , element)
            m2 = str(get_side_area(shape.geometry, get_object_main_axis(get_oriented_bound_box(element)),get_oriented_bound_box(element)))
        if element.is_a("IfcSlab") or element.is_a("IfcCovering"):
            settings = geom.settings()
            settings.set(settings.USE_WORLD_COORDS, True)
            settings.set(settings.DISABLE_OPENING_SUBTRACTIONS, True)
            shape = geom.create_shape(settings , element)
            m2 = str(get_footprint_area(shape.geometry))
    except:
        print(element.GlobalId)
	# Replace dots with commas for Excel
    return m2.replace(".",",")


def get_my_volume(element):
	"""_summary_: Returns the volume value that will be stored into CSV file. Returns empty string if not concerned.
	"""
    m3 = ""
    try:
        if element.is_a("IfcWall") or element.is_a("IfcSlab") or element.is_a("IfcColumn") or element.is_a("IfcBeam"):
            settings = geom.settings()
            settings.set(settings.USE_WORLD_COORDS, True)
            shape = geom.create_shape(settings , element )
            m3 = str(get_volume(shape.geometry))
    except:
        print(element.GlobalId)
	# Replace dots with commas for Excel
    return m3.replace(".",",")


def get_element_line_data(element):
	"""_summary_: Returns the line that will be added to CSV file for a given element. This is where you can tweak to adjust the data you want to export.
	"""
    newLine = ""
    newLine += str(get_parent_building(element)) + ";"
    newLine += str(ifcopenshell.util.element.get_container(element).Name) + ";"
    newLine += str(element.is_a()) + ";"
    newLine += str(element.ObjectType) + ";"
    
    newLine += str(1) + ";"
    newLine += get_my_linear(element) + ";"
    newLine += get_my_net_area(element) + ";"
    newLine += get_my_gross_area(element) + ";"
    newLine += get_my_volume(element) + ";"
    
    newLine += str(element.GlobalId) + ";"
    
    newLine += "\n"
    
    return newLine











allLines = []
tol = 1e-6

# Retrieve all elements important for quantification
walls = model.by_type("IfcWall")
slabs = model.by_type("IfcSlab")
columns = model.by_type("IfcColumn")
beams = model.by_type("IfcBeam")
windows = model.by_type("IfcWindow")
doors = model.by_type("IfcDoor")
ceilings = model.by_type("IfcCovering")

# Add titles to table. Change here too if you tweak exported data lines
allLines.append("Building;Level;IFCType;Type;U;Linear;Net Area;Gross Area;Volume;guid;\n")


# For each element, get specific data and add it to table
for wall in walls:
    allLines.append(get_element_line_data(wall))

for slab in slabs:
    allLines.append(get_element_line_data(slab))

for column in columns:
    allLines.append(get_element_line_data(column))

for beam in beams:
    allLines.append(get_element_line_data(beam))

for window in windows:
    allLines.append(get_element_line_data(window))

for door in doors:
    allLines.append(get_element_line_data(door))

for ceiling in ceilings:
    allLines.append(get_element_line_data(ceiling))



# write data to CSV file
with open( exportedCsvFilePathAndName, 'w' ) as file:
    file.writelines(allLines)