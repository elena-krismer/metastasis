import libsbml
import cobra

"""
model = cobra.io.load_matlab_model("../../obj/models/brain_metastasis.mat")
cobra.io.write_sbml_model(model, "../../obj/models/brain_metastasis_sbml.xml")
"""
model = cobra.io.load_matlab_model("../../obj/models/MDA_MB_231.mat")
cobra.io.write_sbml_model(model, "../../obj/models/MDA_MB_231_sbml.xml")

model = cobra.io.load_matlab_model("../../obj/models/lung_metastasis.mat")
cobra.io.write_sbml_model(model, "../../obj/models/lung_metastasis_sbml.xml")