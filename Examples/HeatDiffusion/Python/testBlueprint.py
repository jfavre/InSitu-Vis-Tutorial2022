import conduit
import conduit.blueprint
import ascent

mesh = conduit.Node()
conduit.blueprint.mesh.examples.braid("uniform", 8, 4, 0, mesh)
print(mesh)
a = ascent.Ascent()
a.open()

actions = conduit.Node()
add_act = actions.append()
add_act["action"] = "add_scenes"
scenes = add_act["scenes"]
scenes["s1/plots/p1/type"] = "pseudocolor"
scenes["s1/plots/p1/field"] = "braid"

a.publish(mesh)
a.execute(actions)

