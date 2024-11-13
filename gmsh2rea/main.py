from read_msh import read_msh
from tet2hex import tet2hex
from export_rea import export_rea
from plot_element import plot_hex_elements, plot_tet_elements, plot_element

msh_file = "./nreltet.msh"
elements, boundaries = read_msh(msh_file)
plot_tet_elements(elements, boundaries)
elements, boundaries = tet2hex(elements, boundaries)
plot_hex_elements(elements, boundaries)
# export_rea('nreltet.rea', elements, boundaries)
