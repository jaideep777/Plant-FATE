import cppimport.import_hook
import plantFATE 

# a = simulator.add(10, 20)
c = plantFATE.Clim()
print(c)
# print(a)
t = plantFATE.LightEnvironment()
print(t.use_ppa)
print(t.n_layers)
print(t.z_star)
print(t.fapar_tot)
print(t.canopy_openness)

t.print()

print(c.co2)

print("TEST the Patch")

s = plantFATE.Patch("tests/params/p.ini")
s.init(1001, 1050)
s.simulate()
s.close
