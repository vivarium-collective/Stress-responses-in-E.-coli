import roadrunner
rr = roadrunner.RoadRunner("BIOMD0000000507.xml")
results = rr.simulate(0, 2000, 200)
rr.plot()
print(results)