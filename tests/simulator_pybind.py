import cppimport.import_hook
import simulator

# a = simulator.add(10, 20)
c = simulator.Clim()
print(c)
# print(a)
t = simulator.TestEnvironment()
print(t.use_ppa)
print(t.n_layers)
print(t.z_star)
print(t.fapar_tot)
print(t.canopy_openness)

t.print()

print(c.co2)
