import gmsh
import yaml
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
    Generates a Gmsh mesh for a CT specimen using parameters from a dictionary.
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
    notch_length = geo_params["notch_length"] * L
    notch_width = geo_params["notch_width"] * L
    
    x_hole = 0.25 * L
    y_mid = height / 2.0
    y_hole_top = y_mid + geo_params["hole_offset"] * L
    y_hole_bot = y_mid - geo_params["hole_offset"] * L

    # 1. Create Base Geometry using OpenCASCADE (occ)
    rect = gmsh.model.occ.addRectangle(0, 0, 0, width, height)
    hole_top = gmsh.model.occ.addDisk(x_hole, y_hole_top, 0, hole_radius, hole_radius)
    hole_bot = gmsh.model.occ.addDisk(x_hole, y_hole_bot, 0, hole_radius, hole_radius)
    notch = gmsh.model.occ.addRectangle(0, y_mid - notch_width/2.0, 0, notch_length, notch_width)

    # Cut holes and notch
    gmsh.model.occ.cut([(2, rect)], [(2, hole_top), (2, hole_bot), (2, notch)])
    gmsh.model.occ.synchronize()

    # 2. Define Physical Groups for Boundaries
    curves = gmsh.model.getEntities(1)
    top_hole_curves, bot_hole_curves = [], []
    
    for dim, tag in curves:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        if (xmin >= x_hole - hole_radius - 1e-3 and xmax <= x_hole + hole_radius + 1e-3 and
            ymin >= y_hole_top - hole_radius - 1e-3 and ymax <= y_hole_top + hole_radius + 1e-3):
            top_hole_curves.append(tag)
        elif (xmin >= x_hole - hole_radius - 1e-3 and xmax <= x_hole + hole_radius + 1e-3 and
              ymin >= y_hole_bot - hole_radius - 1e-3 and ymax <= y_hole_bot + hole_radius + 1e-3):
            bot_hole_curves.append(tag)

    gmsh.model.addPhysicalGroup(1, top_hole_curves, tag=1, name="Top_Hole")
    gmsh.model.addPhysicalGroup(1, bot_hole_curves, tag=2, name="Bottom_Hole")
    
    surfaces = gmsh.model.getEntities(2)
    surface_tags = [tag for dim, tag in surfaces]
    gmsh.model.addPhysicalGroup(2, surface_tags, tag=3, name="Domain")

    # 3. Define Configurable Mesh Refinement (Box Field)
    # This forces h_fine inside the box and smoothly transitions to h_bulk outside
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

    # 4. Generate and Save
    gmsh.model.mesh.generate(2)
    output_filename = mesh_params["output_file"]
    gmsh.write(output_filename)
    gmsh.finalize()
    print(f"Successfully generated {output_filename} with configured refinement.")

if __name__ == "__main__":
    # Example usage
    configuration = load_config("config.yaml")
    generate_ct_mesh(configuration)