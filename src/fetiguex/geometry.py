import gmsh
import yaml
import math
from pathlib import Path

def load_config(config_path="config.yaml"):
    """Loads parameters from a YAML configuration file."""
    config_file = Path(config_path)
    if not config_file.is_file():
        raise FileNotFoundError(f"Configuration file {config_path} not found.")
    
    with open(config_file, "r") as file:
        return yaml.safe_load(file)

def generate_ct_mesh(config):
    """
    Generates a Gmsh mesh for a CT specimen using parameters from a yaml config file.
    Includes configurable local mesh refinement for the phase-field crack path.
    """
    geo_params = config["geometry"]
    mesh_params = config["mesh"]

    gmsh.initialize()
    gmsh.model.add("CT_Specimen")

    # Extract geometry dimensions
    L = geo_params["L"]
    width = geo_params["width_factor"] * L
    height = geo_params["height_factor"] * L
    hole_radius = geo_params["hole_radius"] * L
    
    # Notch parameters
    notch_tip_x = geo_params["notch_tip_x"] * L
    notch_half_w = geo_params["notch_half_width"] * L
    notch_angle_vert_rad = math.radians(geo_params["notch_angle_vertical"])
    
    x_hole = (geo_params["width_factor"] - geo_params["x_hole_from_right"]) * L
    y_mid = height / 2.0
    y_hole_top = y_mid + geo_params["hole_offset"] * L
    y_hole_bot = y_mid - geo_params["hole_offset"] * L

    # Create geometry
    # 1. Create Base Geometry (OpenCASCADE)
    rect = gmsh.model.occ.addRectangle(0, 0, 0, width, height)

    # 2. Construct Top Hole (Split into upper and lower arcs)
    pt_c_top = gmsh.model.occ.addPoint(x_hole, y_hole_top, 0)
    pt_r_top = gmsh.model.occ.addPoint(x_hole + hole_radius, y_hole_top, 0)
    pt_l_top = gmsh.model.occ.addPoint(x_hole - hole_radius, y_hole_top, 0)
    arc_top_upper = gmsh.model.occ.addCircleArc(pt_r_top, pt_c_top, pt_l_top)
    arc_top_lower = gmsh.model.occ.addCircleArc(pt_l_top, pt_c_top, pt_r_top)
    loop_top = gmsh.model.occ.addCurveLoop([arc_top_upper, arc_top_lower])
    hole_top = gmsh.model.occ.addPlaneSurface([loop_top])

    # 3. Construct Bottom Hole (Split into upper and lower arcs)
    pt_c_bot = gmsh.model.occ.addPoint(x_hole, y_hole_bot, 0)
    pt_r_bot = gmsh.model.occ.addPoint(x_hole + hole_radius, y_hole_bot, 0)
    pt_l_bot = gmsh.model.occ.addPoint(x_hole - hole_radius, y_hole_bot, 0)
    arc_bot_upper = gmsh.model.occ.addCircleArc(pt_r_bot, pt_c_bot, pt_l_bot)
    arc_bot_lower = gmsh.model.occ.addCircleArc(pt_l_bot, pt_c_bot, pt_r_bot)
    loop_bot = gmsh.model.occ.addCurveLoop([arc_bot_upper, arc_bot_lower])
    hole_bot = gmsh.model.occ.addPlaneSurface([loop_bot])
    
    # 4. Construct the Exact V-Notch Polygon
    # Calculate horizontal distance from the start of the angle to the tip
    # tan(theta_v) = dx / dy  => dx = dy * tan(theta_v)
    delta_x = notch_half_w * math.tan(notch_angle_vert_rad)
    x_angle_start = notch_tip_x - delta_x

    # Define the 5 perimeter points of the notch
    p1 = gmsh.model.occ.addPoint(0, y_mid + notch_half_w, 0)
    p2 = gmsh.model.occ.addPoint(x_angle_start, y_mid + notch_half_w, 0)
    p3 = gmsh.model.occ.addPoint(notch_tip_x, y_mid, 0) # The sharp tip
    p4 = gmsh.model.occ.addPoint(x_angle_start, y_mid - notch_half_w, 0)
    p5 = gmsh.model.occ.addPoint(0, y_mid - notch_half_w, 0)

    # Connect points into lines
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)
    l4 = gmsh.model.occ.addLine(p4, p5)
    l5 = gmsh.model.occ.addLine(p5, p1)

    # Form the surface to subtract
    loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4, l5])
    notch = gmsh.model.occ.addPlaneSurface([loop])

    # 5. Cut holes and the custom V-notch out of the main rectangle
    gmsh.model.occ.cut([(2, rect)], [(2, hole_top), (2, hole_bot), (2, notch)])
    gmsh.model.occ.synchronize()

    # 6. Define Physical Groups for Boundaries
    curves = gmsh.model.getEntities(1)
    
    top_hole_upper_curves = []
    top_hole_lower_curves = []
    bot_hole_upper_curves = []
    bot_hole_lower_curves = []
    
    for dim, tag in curves:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        y_center = (ymin + ymax) / 2.0  # Find the vertical center of the curve
        
        # Check Top Hole Region
        if (xmin >= x_hole - hole_radius - 1e-3 and xmax <= x_hole + hole_radius + 1e-3 and
            ymin >= y_hole_top - hole_radius - 1e-3 and ymax <= y_hole_top + hole_radius + 1e-3):
            
            if y_center > y_hole_top + 1e-4:
                top_hole_upper_curves.append(tag)
            else:
                top_hole_lower_curves.append(tag)
                
        # Check Bottom Hole Region
        elif (xmin >= x_hole - hole_radius - 1e-3 and xmax <= x_hole + hole_radius + 1e-3 and
              ymin >= y_hole_bot - hole_radius - 1e-3 and ymax <= y_hole_bot + hole_radius + 1e-3):
            
            if y_center > y_hole_bot + 1e-4:
                bot_hole_upper_curves.append(tag)
            else:
                bot_hole_lower_curves.append(tag)

    # Apply strictly separated physical tags
    gmsh.model.addPhysicalGroup(1, top_hole_upper_curves, tag=1, name="Top_Hole_Upper")
    gmsh.model.addPhysicalGroup(1, top_hole_lower_curves, tag=2, name="Top_Hole_Lower")
    gmsh.model.addPhysicalGroup(1, bot_hole_upper_curves, tag=3, name="Bottom_Hole_Upper")
    gmsh.model.addPhysicalGroup(1, bot_hole_lower_curves, tag=4, name="Bottom_Hole_Lower")
    
    surfaces = gmsh.model.getEntities(2)
    surface_tags = [tag for dim, tag in surfaces]
    gmsh.model.addPhysicalGroup(2, surface_tags, tag=5, name="Domain")

    # 7. Define Configurable Mesh Refinement (Box Field)
    ref_box = mesh_params["refinement_box"]
    
    gmsh.model.mesh.field.add("Box", 1)
    gmsh.model.mesh.field.setNumber(1, "VIn", mesh_params["h_fine"])
    gmsh.model.mesh.field.setNumber(1, "VOut", mesh_params["h_bulk"])
    gmsh.model.mesh.field.setNumber(1, "XMin", ref_box["x_min"])
    gmsh.model.mesh.field.setNumber(1, "XMax", ref_box["x_max"])
    gmsh.model.mesh.field.setNumber(1, "YMin", ref_box["y_min"])
    gmsh.model.mesh.field.setNumber(1, "YMax", ref_box["y_max"])
    gmsh.model.mesh.field.setNumber(1, "Thickness", ref_box["thickness"])

    gmsh.model.mesh.field.setAsBackgroundMesh(1)

    # 8. Generate and Save Mesh
    gmsh.model.mesh.generate(2)
    output_filename = mesh_params["output_file"]
    gmsh.write(output_filename)
    gmsh.finalize()
    print(f"Successfully generated mesh file: {output_filename}")

if __name__ == "__main__":
    configuration = load_config("config.yaml")
    generate_ct_mesh(configuration)