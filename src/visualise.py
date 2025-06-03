import py3Dmol

with open("./data/output/5awl_solv_ions.gro", "r") as f:
    gro_data = f.read()

view = py3Dmol.view(width=800, height=600)
view.addModel(gro_data, "gro")
view.setStyle({"cartoon": {}})
view.addStyle({"resn": "NA"}, {"sphere": {"color": "yellow"}})
view.addStyle({"resn": "CL"}, {"sphere": {"color": "cyan"}})
view.addStyle({"resn": "SOL"}, {"stick": {"colorscheme": "greenCarbon"}})
view.setBackgroundColor("white")
view.zoomTo()

with open("data/output/structure.html", "w") as out:
    out.write(view._make_html())
    print("3D structure saved to data/output/structure.html")