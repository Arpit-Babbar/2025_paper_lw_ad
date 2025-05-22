using Trixi2Vtk
using ReadVTK
using TenkaiAD
using Plots

tenkaiad_data_dir = joinpath(dirname(pathof(TenkaiAD)), "..", "paper_data")
m2000_outdir = joinpath(tenkaiad_data_dir, "output_m2000")
trixi2vtk(joinpath(m2000_outdir, "solution_018086.h5"), output_directory = m2000_outdir)
vtk_file = joinpath(m2000_outdir, "solution_018086.vtu")
vtk = VTKFile(vtk_file)
data = get_point_data(vtk)

points = get_points(vtk)
points_2d = @view points[1:2, :]
cells = get_cells(vtk)

# Contains things beyond points
point_data = ReadVTK.find_element(ReadVTK.piece(vtk), "PointData")
xml_data_array = ReadVTK.find_element(point_data, "DataArray") # Returns density for some reason, maybe because it is the first data_array
data_array = VTKDataArray(xml_data_array, vtk)
density = get_data(data_array)
ReadVTK.find_element(point_data, "Density")

x = points_2d[1, :]
y = points_2d[2, :]
z = density
